// Copyright (c) 2019-3-3 maxx
#include "RooP4Vector.hh"
#include "TComplex.h"
Double_t RooP4Vector::dot(const Double_t p4a[4], const Double_t p4b[4]) {
    return p4a[0]*p4b[0] - p4a[1]*p4b[1]
        - p4a[2]*p4b[2] - p4a[3]*p4b[3];
}

void RooP4Vector::calt1(const Double_t daug1[], const Double_t daug2[],
        Double_t t1[]) {
    Double_t pa[4] = { daug1[0] + daug2[0], daug1[1] + daug2[1],
        daug1[2] + daug2[2], daug1[3] + daug2[3] };
    Double_t qa[4] = { daug1[0] - daug2[0], daug1[1] - daug2[1],
        daug1[2] - daug2[2], daug1[3] - daug2[3] };

    Double_t pp, pq;
    pp = RooP4Vector::dot(pa, pa);
    pq = RooP4Vector::dot(pa, qa);
    t1[0] = qa[0] - pq/pp*pa[0];
    t1[1] = qa[1] - pq/pp*pa[1];
    t1[2] = qa[2] - pq/pp*pa[2];
    t1[3] = qa[3] - pq/pp*pa[3];
}

void RooP4Vector::calt2(const Double_t daug1[], const Double_t daug2[],
        Double_t t2[][4]) {
    Double_t pp, r;
    Double_t pa[4], t1[4];

    calt1(daug1, daug2, t1);
    r = RooP4Vector::dot(t1, t1);

    for (Int_t i=0; i != 4; i++) {
        pa[i] = daug1[i] + daug2[i];
    }

    pp = RooP4Vector::dot(pa, pa);
    for (Int_t i=0; i != 4; i++) {
        for (Int_t j=0; j != 4; j++) {
            t2[i][j] = t1[i]*t1[j] -
                1.0/3*r*(RooP4Vector::G[i][j]-pa[i]*pa[j]/pp);
        }
    }
}
