/* ==================================== *
 * Copyright (C) 2022 NextMove Software *
 * ==================================== */
// rdpath_hash.h
// see rdpath_dropin.h, but here we change the hash
// encoding of paths and capture all rings 3, 4, 5, ..
// this implementation produces a different fingerprint

// JWM: It should be easy to add 
// - custom atom invariants 
// - atomToBit mapping
// - multiple bits per hash (note: I think this is not useful)
#ifndef NMS_RDPATH_HASH_H
#define NMS_RDPATH_HASH_H

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
  const RDKit::Bond *bonds[MAX_PATH_CAP];
  unsigned int path[MAX_PATH_CAP];
  const RDKit::ROMol *mol; 
  ExplicitBitVect *bit_vect;
  unsigned char *avisit = (unsigned char*)0;
  unsigned char *bondHashes = (unsigned char*)0;
  unsigned int acap = 0; 
  unsigned int bcap = 0; 
  unsigned int len = 0;
  unsigned int min_len;
  unsigned int max_len;
  unsigned int fp_size;
  unsigned int fp_mask;

  RdPathFp(unsigned int min_len_ = DEFAULT_MIN_LEN,
           unsigned int max_len_ = DEFAULT_MAX_LEN,
           unsigned int fp_size_ = DEFAULT_FP_SIZE)
    : min_len(min_len_), max_len(max_len_), fp_size(fp_size_), fp_mask(fp_size-1) {
  }

  ~RdPathFp() {
    free(avisit);
    free(bondHashes);
  }

  bool init(ExplicitBitVect *bit_vect_, const RDKit::ROMol *mol_)
  {
    if (min_len < 0 || min_len > max_len || max_len > MAX_PATH_CAP)
      return false;
    mol = mol_;  
    bit_vect = bit_vect_;
    if (mol->getNumAtoms() > acap) {
      avisit = (unsigned char*)realloc(avisit, mol->getNumAtoms());
      memset(avisit,0,mol->getNumAtoms());
      acap = mol->getNumAtoms();
    }
    if (mol->getNumBonds() > bcap) {
      bondHashes = (unsigned char*)realloc(bondHashes, mol->getNumBonds());
      for (const RDKit::Bond *bond : mol->bonds())
        bondHashes[bond->getIdx()] = hash_bond(bond);      
    }
    return true;
  }

  static int comp(const unsigned int *arr, unsigned char i, unsigned char j)
  {
    return arr[j] - arr[i];
  }

  // when reversed is array lexicographically better (-1), same (0)
  // or worse (+1)
  static int comp_rev(const unsigned int *arr, const unsigned char len)
  {
    unsigned char i = 0, j = len-1;
    int cmp;
    while (i < j)
      if ((cmp = comp(arr, i++, j--)))
        return cmp;
    return 0; // symmetric
  }

  // hash array values from begin to end (foward)
  static unsigned int hash_fwd(const unsigned int *arr, const unsigned int len)
  {
    unsigned int hash = len;
    for (unsigned int i = 0; i < len; i++)
      gboost::hash_combine(hash, arr[i]);
    return hash;      
  }

  // hash array values from end to being (reverse)
  static unsigned int hash_rev(const unsigned int *arr, const unsigned int len)
  {
    unsigned int hash = len;
    for (unsigned int i = 1; i <= len; i++)
      gboost::hash_combine(hash, arr[len-i]);
    return hash;      
  }

  void encode_path() {
    unsigned int hash;
    if (comp_rev(path, len) > 0)
      hash = hash_fwd(path, len);
    else
      hash = hash_rev(path, len);
    bit_vect->setBit(hash & fp_mask);
  }

  uint32_t hash_bond(const RDKit::Bond *bond)
  {
    uint32_t hash = 13;
    const RDKit::Atom *beg = bond->getBeginAtom();
    const RDKit::Atom *end = bond->getEndAtom();
    uint32_t bhash = 2*beg->getAtomicNum() + beg->getIsAromatic();
    uint32_t ehash = 2*end->getAtomicNum() + end->getIsAromatic();
    uint32_t bdeg = beg->getDegree();
    uint32_t edeg = end->getDegree();
    if (bhash < ehash || (bhash == ehash && bdeg < edeg)) {
      gboost::hash_combine(hash, bhash);
      gboost::hash_combine(hash, bdeg);
      gboost::hash_combine(hash, ehash);
      gboost::hash_combine(hash, edeg);
    } else {
      gboost::hash_combine(hash, ehash);
      gboost::hash_combine(hash, edeg);
      gboost::hash_combine(hash, bhash);
      gboost::hash_combine(hash, bdeg);
    }
    if (bond->getIsAromatic())
      gboost::hash_combine(hash, RDKit::Bond::AROMATIC);
    else  
      gboost::hash_combine(hash, bond->getBondType());
    return hash;
  }

  bool push(const RDKit::Bond *bond, bool do_encode=true) {
    bonds[len] = bond;
    path[len] = bondHashes[bond->getIdx()];
    len++;
    if (len >= min_len && do_encode)
      encode_path();
    return len < max_len;  
  }

  void pop() {
    --len;
  }

  // determine if an array when rotated and optionally reversed is 
  // lexicographically better (-1), equal (0) or worse (+1).
  // * args a and b can be in the range: 1 < a/b < len+1 and -(len+1) < b < 1
  // * if -ve it indicates the comparison should be reversed
  static int comp_off(const unsigned int *ptr, const unsigned char len, 
                      const char a, const char b)
  { 
    int cmp;
    if (a < 0 && b < 0) {
      for (unsigned char i=0; i<len; i++) {
        if ((cmp = comp(ptr, (len-a-i-1)%len, (len-b-i-1)%len)))
          return cmp;
      }  
    } else if (a < 0) {
      for (unsigned char i=0; i<len; i++) {
        if ((cmp = comp(ptr, (len-a-i-1)%len, (i+b-1)%len)))
          return cmp;
      }  
    } else if (b < 0) { 
      for (unsigned char i=0; i<len; i++) {
        if ((cmp = comp(ptr, (i+a-1)%len, (len-b-i-1)%len)))
          return cmp;
      }  
    } else {
      for (unsigned char i=0; i<len; i++) {
        if ((cmp = comp(ptr, (i+a-1)%len, (i+b-1)%len)))
          return cmp;
      }  
    }
    return cmp;
  }

  // determine the canonical ring offset (see ring_cansmi)
  static char ring_offset(const unsigned int *ptr, const unsigned char len)
  {
    char best;
  #if 1
    // find the index of the max atom value to start from
    unsigned char max = 0;
    for (u_int8_t i=2; i<len; i+=2)
      if (ptr[i] > ptr[max])
        max = i;
    
    // determine which ordering is best
    if (comp_off(ptr, len, -(max+1), max+1) < 0)
      best = -(max+1);
    else
      best = max+1;
    
    // check for alternate starting points 
    for (u_int8_t i=max+2; i<len; i+=2) {
      if (ptr[i] == ptr[max]) {
        if (comp_off(ptr, len, i+1, best) < 0)
          best = i+1;
        if (comp_off(ptr, len, -(i+1), best) < 0)
          best = -(i+1);
      }
    }
  #else
    // do all the comparisons - old-school 
    best = comp_off(ptr, len, -1, +1) ? -1 : +1;
    for (u_int8_t i=2; i<len; i+=2) {
      if (comp_off(ptr, len, i+1, best) < 0)
        best = i+1;
      if (comp_off(ptr, len, -(i+1), best) < 0)
        best = -(i+1);
    }
  #endif  
    return best;
  }

  void encode_ring() {
    
#if 1 // disable to check ring offset works correctly
    // lexicographic check
    for (unsigned int i = 1; i < len; i++)
      if (bonds[i] > bonds[0])
        return;
    if (bonds[1] > bonds[len-1])    
      return;
#endif      

    const int8_t r = ring_offset(path, len);
    unsigned int hash = len;
    gboost::hash_combine(hash, 1); // differentiate ring/path
    if (r < 0) {
      gboost::hash_combine(hash, path[(len-r-1)%len]);
      for (uint8_t i=1; i<len; i++)
        gboost::hash_combine(hash, path[(len-r-i-1)%len]);      
    } else {
      gboost::hash_combine(hash, path[r-1]);
      for (uint8_t i=1; i<len; i++)
        gboost::hash_combine(hash, path[(r+i-1)%len]);
    }

    bit_vect->setBit(hash & fp_mask);
  }

  void encode_ring(const RDKit::Bond *bond) {
    if (len+1 < min_len)
      return;
    bonds[len] = bond;
    path[len] = bondHashes[bond->getIdx()];
    len++;
    encode_ring();
    len--;
  }

  void traverse(const RDKit::Atom *root, const RDKit::Atom *atom, const RDKit::Bond *prev) {
    if (push(prev, root < atom)) {
      avisit[atom->getIdx()] = 1;
      for (const RDKit::Bond *bond : mol->atomBonds(atom)) {
        if (bond == prev)
          continue;
        RDKit::Atom *nbor = bond->getOtherAtom(atom);
        if (!avisit[nbor->getIdx()])
          traverse(root, nbor, bond);    
        else if (avisit[nbor->getIdx()] == 2)
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