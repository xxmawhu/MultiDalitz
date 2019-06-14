// Copyright (c) 2019-3-3 maxx
#ifndef INC_ROOSPINFACTOR_H_
#define INC_ROOSPINFACTOR_H_

#include <math.h>
#include <iostream>
#include <fstream>
#include "TNamed.h"

#if defined(USEROOT) || defined(__CINT__)
#else
#endif

namespace RooSpinFactor{
    Double_t SpinFactor(const Int_t & angL, const Double_t *p4_1,
            const Double_t *p4_2, const Double_t *p4_3);
};  // namespace RooSpinFactor

#endif
