// Copyright (c) 2019-3-3 maxx
#ifndef DecayType_HH
#define DecayType_HH

#include <math.h>
#include <iostream>
#include <fstream>
#include "TNamed.h"

#if defined(USEROOT) || defined(__CINT__)
#else
#endif

namespace DecayType{
    enum  DecayType{
        AB = 3,
        AC = 2,
        BC = 1,
        BA = 3,
        CB = 1,
        CA = 2
    };
};

#endif
