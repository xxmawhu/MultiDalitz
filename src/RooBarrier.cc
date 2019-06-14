// Copyright (c) 2019-3-3 maxx
#include "RooP4Vector.hh"
#include "RooBarrier.hh"
#include "TComplex.h"

Double_t RooBarrier::BreakMomenta(const Double_t &Sr, const Double_t &Sa,
        const Double_t &Sb) {
    Double_t Qabcs = (Sr+Sa-Sb)*(Sr+Sa-Sb)/(4*Sr)-Sa;
    if (Qabcs < 0) Qabcs = 1e-16;
    return sqrt(Qabcs);
}
Double_t RooBarrier::BreakMomentaSq(const Double_t &Sr, const Double_t &Sa,
        const Double_t &Sb) {
    Double_t Qabcs = (Sr+Sa-Sb)*(Sr+Sa-Sb)/(4*Sr)-Sa;
    if (Qabcs < 0) Qabcs = 1e-16;
    return Qabcs;
}

Double_t RooBarrier::Factor(const Int_t &l, const Double_t p4_1[4],
        const Double_t p4_2[4], const Double_t &r) {
    Double_t p4R[4] = { p4_1[0]+ p4_2[0], p4_1[1]+p4_2[1],
        p4_1[2]+p4_2[2], p4_1[3]+p4_2[3] };
    Double_t q = RooBarrier::BreakMomenta(RooP4Vector::dot(p4R, p4R),
            RooP4Vector::dot(p4_1, p4_1),
            RooP4Vector::dot(p4_2, p4_2));
    q = sqrt(q);
    Double_t z = q*r;
    z = z*z;
    Double_t F = 1;
    if (l > 2) F = 0;
    if (l == 1) {
        F = sqrt((2*z)/(1+z));
    }
    if (l == 2) {
        F = sqrt((13*z*z)/(9+3*z+z*z));
    }
    return F;
}

