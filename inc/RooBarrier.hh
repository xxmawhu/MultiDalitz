// Copyright (c) 2019-3-3 maxx
#ifndef Roo_barrierFactor_HH
#define Roo_barrierFactor_HH

#include <math.h>
#include <iostream>
#include <fstream>
#include "TNamed.h"

#if defined(USEROOT) || defined(__CINT__)
#else
#endif

namespace RooBarrier{
    Double_t Factor(const Int_t & angL, const Double_t *p4_1,
            const Double_t *p4_2, const Double_t &r = 5.0);
    Double_t BreakMomenta(const Double_t &mRSq, const Double_t& m1Sq,
            const Double_t& m2Sq);
    Double_t BreakMomentaSq(const Double_t &mRSq, const Double_t& m1Sq,
            const Double_t& m2Sq);
};  // namespace RooBarrier

#endif
