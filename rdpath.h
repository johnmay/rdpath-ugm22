/* =====================================
 * Copyright (C) 2022 NextMove Software 
 * ===================================== */
#ifndef NMS_RDPATH_H
#define NMS_RDPATH_H

#define NEW_HASH

#ifdef NEW_HASH
#include "rdpath_hash.h"
#else
#include "rdpath_dropin.h"
#endif

static bool EncodePathFp(ExplicitBitVect *bv, 
                         RDKit::ROMol &mol, 
                         unsigned int fp_size=2048) {
  if ((fp_size & (fp_size-1)) != 0)                          
    return false; // must be pow-of-two!
  RdPathFp s(1,7,fp_size);
  if (!s.init(bv, &mol))
    return false;
  s.encode_mol();
  return true;
} 

static ExplicitBitVect *EncodePathFp(RDKit::ROMol &mol, 
                                     unsigned int fp_size=2048) {
  ExplicitBitVect *bit_vect = new ExplicitBitVect(fp_size);
  EncodePathFp(bit_vect, mol, fp_size);
  return bit_vect;
}

#endif