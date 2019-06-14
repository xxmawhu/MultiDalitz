// Copyright (c) 2019-3-3 maxx
#include "RooP4Vector.hh"
#include "RooSpinFactor.hh"
#include "TComplex.h"

Double_t RooSpinFactor::SpinFactor(const Int_t & Ang, const Double_t *p4_1,
        const Double_t *p4_2, const Double_t *p4_solo) {
    if (Ang == 0) {
        return 1;
    }
    Double_t temp_PDF = 0;
    Double_t pR[4], pD[4];
    for (Int_t i=0; i != 4; i++) {
        pR[i] = p4_1[i] + p4_2[i];
        pD[i] = pR[i] + p4_solo[i];
    }
    Double_t t1[4], T1[4];
    RooP4Vector::calt1(p4_1, p4_2, t1);
    RooP4Vector::calt1(pR, p4_solo, T1);
    for (Int_t i = 0; i < 4; i++) {
        temp_PDF += t1[i]*T1[i]*RooP4Vector::G[i][i];
    }
    if (Ang == 1) {
        return temp_PDF;
    }
    Double_t t2[4][4], T2[4][4];
    Double_t sa[2], sb[2], sc[2], B[2];
        RooP4Vector::calt2(p4_1, p4_2, t2);
        RooP4Vector::calt2(pR, p4_solo, T2);
        for (Int_t i = 0; i < 4; i++) {
            for (Int_t j = 0; j < 4; j++) {
                temp_PDF += t2[i][j] * T2[j][i] *
                    RooP4Vector::G[i][i] * RooP4Vector::G[j][j];
            }
        }
    if (Ang == 2) {
        return temp_PDF;
    }
    return temp_PDF;
}
