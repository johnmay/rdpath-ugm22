/* ==================================== *
 * Copyright (C) 2022 NextMove Software *
 * ==================================== */
// rdpath_dropin.h
// a drop in replacement for RDKit OG FP with the following
// options branched=false, bitsPerHash=1 feature

// JWM: It should be easy to add 
// - custom atom invariants 
// - atomToBit mapping
// - multiple bits per hash (note: I think this is not useful)

#ifndef NMS_RDPATH_DROPIN_H
#define NMS_RDPATH_DROPIN_H

#include <stdlib.h>

#include <algorithm>

#include <GraphMol/GraphMol.h>
#include <DataStructs/ExplicitBitVect.h>
#include <RDGeneral/hash/hash.hpp>

#define MAX_PATH_CAP 16
#define DEFAULT_MIN_LEN 1
#define DEFAULT_MAX_LEN 7
#define DEFAULT_FP_SIZE 2048

struct RdPathFp {
  std::vector<unsigned int> bondHashes;
  const RDKit::Bond *bonds[MAX_PATH_CAP];
  const RDKit::ROMol *mol; 
  ExplicitBitVect *bit_vect;
  unsigned char *avisit = (unsigned char*)0;
  unsigned char *degree = (unsigned char*)0; 
  unsigned int cap = 0; 
  unsigned int len = 0;
  unsigned int min_len;
  unsigned int max_len;
  unsigned int fp_size;
  unsigned int fp_mask;

  RdPathFp(unsigned int min_len_ = DEFAULT_MIN_LEN,
           unsigned int max_len_ = DEFAULT_MAX_LEN,
           unsigned int fp_size_ = DEFAULT_FP_SIZE)
    : min_len(min_len_), max_len(max_len_), 
      fp_size(fp_size_), fp_mask(fp_size-1) {
  }

  ~RdPathFp() {
    free(avisit);
  }

  bool init(ExplicitBitVect *bit_vect_, const RDKit::ROMol *mol_)
  {
    if (min_len < 0 || min_len > max_len || max_len > MAX_PATH_CAP)
      return false;
    mol = mol_;  
    bit_vect = bit_vect_;
    if (mol->getNumAtoms() > cap) {
      cap = mol->getNumAtoms();
      avisit = (unsigned char*)realloc(avisit, 2*cap);
      degree = avisit + cap;      
      memset(avisit,0,2*cap);
    }
    return true;
  }

  // this encoding of path matches the generateBondHashes()
  // except, the atom degree is tracked as we goo,
  // we can derive the bnbor (bond neighbors) based
  // on the atom degree in O(1), likewise we can can track
  // the number of each degree in the path (degree histogram)
  // and determine the number of unique atoms from this
  void encode() 
  {    
    bondHashes.clear();
    unsigned char deg_hist[8] = {0};
    for (unsigned int i = 0; i < len; i++) {
      const RDKit::Bond *bond = bonds[i];

      RDKit::Atom *beg = bond->getBeginAtom();
      RDKit::Atom *end = bond->getEndAtom();
      
      unsigned int btype = bond->getBondType();
      if (bond->getIsAromatic())
        btype = RDKit::Bond::AROMATIC;
      unsigned int atype1 = 2*beg->getAtomicNum() + beg->getIsAromatic();
      unsigned int adeg1  = degree[beg->getIdx()];
      unsigned int atype2 = 2*end->getAtomicNum() + end->getIsAromatic();
      unsigned int adeg2  = degree[end->getIdx()];
      unsigned int bnbor  = (adeg1-1) + (adeg2-1);
      deg_hist[adeg1]++;
      deg_hist[adeg2]++;

      if (atype1 < atype2) {
        std::swap(atype1, atype2);
        std::swap(adeg1, adeg2);
      } else if (atype1 == atype2 && adeg1 < adeg2) {
        std::swap(adeg1, adeg2);
      }
      
      std::uint32_t hash = bnbor;
      gboost::hash_combine(hash, btype);
      gboost::hash_combine(hash, atype1);
      gboost::hash_combine(hash, adeg1);
      gboost::hash_combine(hash, atype2);
      gboost::hash_combine(hash, adeg2);
      bondHashes.push_back(hash);
    }

    std::hash_result_t seed;
    if (bondHashes.size() > 1) {
      std::sort(bondHashes.begin(), bondHashes.end());
      seed = 0;
      for (unsigned int i=0; i<bondHashes.size(); i++)
        gboost::hash_combine(seed, bondHashes[i]);
      
      // numUniqueAtoms = deg_hist[1] + deg_hist[2]/2 + deg_hist[3]/3 + ...
      // but for our specialized case of paths a simpler sum can be used
      unsigned char numUniqueAtoms = len + (deg_hist[1]!=0) - (deg_hist[3]!=0);

      gboost::hash_combine(seed, numUniqueAtoms);
    } else {
      seed = bondHashes[0];
    }
    bit_vect->setBit(seed & fp_mask);
  }

  bool push(const RDKit::Bond *bond, 
            bool do_encode=true) 
  {
    bonds[len++] = bond;
    degree[bond->getBeginAtomIdx()]++;
    degree[bond->getEndAtomIdx()]++;
    if (len >= min_len && do_encode)
      encode();
    return len < max_len;  
  }

  void pop() 
  {
    const RDKit::Bond *bond = bonds[--len];
    degree[bond->getBeginAtomIdx()]--;
    degree[bond->getEndAtomIdx()]--;
  }

  void encode_ring(const RDKit::Bond *bond) 
  {
    // jwm: cannot do an ordering check here as may gave D3 atoms
    push(bond);
    pop();
  }

  void traverse(const RDKit::Atom *root, 
                const RDKit::Atom *atom, 
                const RDKit::Bond *prev) 
  {
    // root < atom, lexicographic ordering hashes each 'path' once
    if (push(prev, root < atom)) {
      avisit[atom->getIdx()] = 1;
      for (const RDKit::Bond *bond : mol->atomBonds(atom)) {
        if (bond == prev)
          continue;
        RDKit::Atom *nbor = bond->getOtherAtom(atom);
        if (!avisit[nbor->getIdx()])
          traverse(root, nbor, bond);    
        else if (len+1 == max_len)
          encode_ring(bond);  
      }
      avisit[atom->getIdx()] = 0;
    }
    pop();
  }

  void encode_mol() 
  {
    for (const RDKit::Atom *atom : mol->atoms()) {
      avisit[atom->getIdx()] = 2;
      for (const RDKit::Bond *bond : mol->atomBonds(atom))
        traverse(atom, bond->getOtherAtom(atom), bond);       
      avisit[atom->getIdx()] = 0;
    }
  }
};

#endif