// Copyright (c) 2019-3-3 maxx
#ifndef LineShape_HH
#define LineShape_HH

#include <iostream>
#if defined(USEROOT) || defined(__CINT__)
#else
#endif
namespace LineShape{
    enum Shape{
        Flat,
        RBW,
        a980_p,
        a980_p_3C,
        a980_0,
        K1430_p,
        K1430_0,
        GS
    };
};  // namespace LineShape

#endif
