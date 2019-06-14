// Copyright (c) 2019-3-3 maxx
#ifndef Roo_Propagator_HH
#define Roo_Propagator_HH

#include <math.h>
#include <iostream>

#if defined(USEROOT) || defined(__CINT__)
#else
#endif
#include "TComplex.h"
#include <vector>
#include "LineShape.hh"
using std::vector;
namespace Propagator {
    TComplex getVal(const LineShape::Shape &lineshape, const Double_t &s,
            const Double_t &Sa, const Double_t &Sb,
            const vector<Double_t> &parameters);


    Int_t GetParaSize(const LineShape::Shape &shape);


namespace BW {
    // parameters order:
    // spin, mass, width
    TComplex getVal(const Double_t &s,
            const vector<Double_t> &parameters);
};  // namespace BW


namespace RBW {
    // parameters order: spin, radius, mass, width
    TComplex getVal(const Double_t &s, const Double_t & Sa,
            const Double_t &Sb, const vector<Double_t> &parameters);
    Double_t  width(const Double_t &mass, const Double_t &sa,
            const Double_t &sb, const Double_t &sc, const Double_t &r,
            const Int_t &l);
};  // namespace RBW


namespace Flatte {
    TComplex getVal(const Double_t &m0Sq, const Double_t &s,
            const Double_t &g1,
            const Double_t &mSqa1,
            const Double_t &mSqb1,
            const Double_t g2,
            const Double_t &mSqa2,
            const Double_t &mSqb2,
            const Double_t &g3 = 0.0,
            const Double_t &mSqa3 = 0.1,
            const Double_t &mSqb3 = 0.1);
    TComplex rho(const Double_t &Sr, const Double_t &Sa,
            const Double_t &Sb);
};  // namespace Flatte


namespace K1430_0 {
        // parameters order:
        // spin, r, mass, g(KPi)^2, g(Pi Eta)^2
        TComplex getVal(const Double_t &s,
                const vector<Double_t> &parameters);
};  // namespace K1430_0


namespace K1430_p {
        // parameters order:
        // spin, r, mass, g(KPi)^2, g(Pi Eta)^2
        TComplex getVal(const Double_t &s,
                const vector<Double_t> &parameters);
};  // namespace K1430_p


namespace a980_p {
    // parameters order:
    // spin, r, mass, g(pi eta)^2, g(K K0)^2
    TComplex getVal(const Double_t &s,
                const vector<Double_t> &parameters);
};  // namespace a980_p

namespace a980_p_3C {
        // parameters order:
        // spin, r, mass, g(pi eta)^2, g(K K0)^2, g(pi+ eta')
        TComplex getVal(const Double_t &s,
                const vector<Double_t> &parameters);
};  // namespace a980_p_3C


namespace a980_0 {
    // parameters order:
    // spin, r, mass, g(pi eta)^2, g(K K0)^2
    TComplex getVal(const Double_t &s,
            const vector<Double_t> &parameters);
};  // namespace a980_0


namespace GS {
        // parameters order
        // spin, r, mass, width
        Double_t h(const Double_t &, const Double_t&);
        Double_t dh(const Double_t &m0, const Double_t& q0);
        Double_t f(const Double_t& s, const Double_t& m0, const Double_t
                &width, const Double_t &q0, const Double_t& _q);
        Double_t d(const Double_t& mass, const Double_t &q0);
        TComplex getVal(const Double_t& s, const Double_t &Sa, const Double_t
                &Sb,
                const vector<Double_t> &parameters);
};  // namespace GS


};  // namespace Propagator

#endif
