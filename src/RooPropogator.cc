// Copyright (c) 2019-3-3 maxx
#include <math.h>
#include "RooPropogator.hh"
#include "LineShape.hh"
#include "RooP4Vector.hh"
#include "RooBarrier.hh"
#include "LauConstants.hh"
#include "TComplex.h"
TComplex Propagator::getVal(const LineShape::Shape &lineShpe,
        const Double_t &s, const Double_t &Sa, const Double_t &Sb,
       const vector<Double_t> &parameters) {
    switch (lineShpe) {
        case LineShape::RBW:
            return Propagator::RBW::getVal(s, Sa, Sb, parameters);
        case LineShape::a980_p:
            return Propagator::a980_p::getVal(s, parameters);
        case LineShape::a980_p_3C:
            return Propagator::a980_p_3C::getVal(s, parameters);
        case LineShape::a980_0:
            return Propagator::a980_0::getVal(s, parameters);
        case LineShape::K1430_p:
            return Propagator::K1430_p::getVal(s, parameters);
        case LineShape::K1430_0:
            return Propagator::K1430_0::getVal(s, parameters);
        case LineShape::GS:
            return Propagator::GS::getVal(s, Sa,  Sb,  parameters);
        default:  // flat, noneresoance
            return TComplex(1.0, 0);
    }
    return TComplex(1.0, 0);
}
Int_t Propagator::GetParaSize(const LineShape::Shape &shape) {
    switch (shape) {
        case LineShape::a980_0:
            return 5;
        case LineShape::a980_p:
            return 5;
        case LineShape::K1430_0:
            return 6;
        case LineShape::K1430_p:
              return 6;
        case LineShape::Flat:
            return 1;
        case LineShape::RBW:
            return 4;
        case LineShape::a980_p_3C:
            return 6;
        default:
            return 4;
    }
}

Double_t Propagator::RBW::width(const Double_t & mass, const Double_t &
        Sr, const Double_t & Sa, const Double_t & Sb, const Double_t & r,
        const Int_t &l) {
    Double_t widm(0.), q(0.), q0(0.);
    // // cout<<__LINE__<<__func__<<endl;
    Double_t Sr0 = mass*mass;
    Double_t m = sqrt(Sr);
    q = RooBarrier::BreakMomentaSq(Sr, Sa, Sb);
    q0 = RooBarrier::BreakMomentaSq(Sr0, Sa, Sb);
    // // cout<<__LINE__<<__func__<<endl;

    Double_t z, z0;
    z = q*r*r;
    z0 = q0*r*r;
    Double_t t = q/q0;
    Double_t F(0.);
    if (l == 0) F = 1;
    if (l == 1) F = sqrt ( (1+z0)/ (1+z));
    if (l == 2) F = sqrt ( (9+3*z0+z0*z0)/ (9+3*z+z*z));
    // // cout<<__LINE__<<__func__<<endl;
    widm = pow(t, l+0.5)*mass/m*F*F;
    return widm;
}

TComplex Propagator::RBW::getVal(const Double_t &s, const Double_t &Sa,
        const Double_t &Sb,
        const vector<Double_t> & paras) {
    // // cout<<"paras size: "<<paras.size()<<endl;
    TComplex ci(0, 1);
    // // cout<<__LINE__<<__func__<<endl;
    Double_t width = Propagator::RBW::width(paras[2],
            s, Sa, Sb, paras[1], paras[0]);
    // // cout<<__LINE__<<__func__<<endl;
    TComplex prop = 1.0 / (paras[2] * paras[2] - s -
            ci * paras[2] * width * paras[3]);
    // // cout<<__LINE__<<__func__<<endl;
    return prop;
}
// 2018-03-28 20:18
TComplex Propagator::Flatte::rho(const Double_t &Sr,
        const Double_t &Sa, const Double_t &Sb) {
    Double_t q = (Sr+Sa-Sb)*(Sr+Sa-Sb)/(4*Sr)-Sa;
    TComplex rho;
    TComplex ci(0, 1);
    if (q > 0) rho = TComplex::One ()*sqrt (q/Sr);
    if (q < 0) rho = ci*sqrt (-q/Sr);
    rho = 2.0*rho;
    return rho;
}
// 2018-03-28 20:20
TComplex Propagator::Flatte::getVal(const Double_t &m0Sq, const Double_t &s,
                const Double_t &g1,
                const Double_t &mSqa1,
                const Double_t &mSqb1,
                const Double_t g2,
                const Double_t &mSqa2,
                const Double_t &mSqb2,
                const Double_t &g3,
                const Double_t &mSqa3,
                const Double_t &mSqb3) {
    TComplex ci(0, 1);
    TComplex rho1 = Propagator::Flatte::rho(s, mSqa1, mSqb1);
    TComplex rho2 = Propagator::Flatte::rho(s, mSqa2, mSqb2);
    TComplex rho3 = Propagator::Flatte::rho(s, mSqa3, mSqb3);
    // cout<<"rho1: "<<rho1<<endl;
    // cout<<"rho2: "<<rho2<<endl;
    // cout<<"rho3: "<<rho3<<endl;
    // cout<<"g1*g2: "<<g1*g1<<endl;
    // cout<<"g2:g2: "<<g2*g2<<endl;
    // cout<<"g3:g3: "<<g3*g3<<endl;
    TComplex prop = 1.0/(m0Sq - s - ci*(g1*g1*rho1 + g2*g2*rho2 +
                g3*g3*rho3) );
    return prop;
}
TComplex Propagator::a980_p::getVal(const Double_t &s,
        const vector<Double_t> & paras) {
    return Propagator::Flatte::getVal(paras[2]*paras[2], s,
            paras[3], LauConstants::mPiSq, LauConstants::mEtaSq,
            paras[4], LauConstants::mKSq, LauConstants::mK0Sq);
}
TComplex Propagator::a980_0::getVal(const Double_t &s,
        const vector<Double_t> & paras) {
    return Propagator::Flatte::getVal(
            paras[2]*paras[2], s,
            paras[3], LauConstants::mPi0Sq, LauConstants::mEtaSq,
            paras[4], LauConstants::mKSq, LauConstants::mKSq);
}
TComplex Propagator::K1430_p::getVal(const Double_t &s,
        const vector<Double_t> & paras) {
    // cout<<__func__ <<endl;
    // cout<<"value"<<
    return Propagator::Flatte::getVal(
            paras[2]*paras[2], s,
            paras[3], LauConstants::mKSq, LauConstants::mPi0Sq,
            paras[4], LauConstants::mKSq, LauConstants::mEtaPrimeSq,
            paras[5], LauConstants::mKSq, LauConstants::mEtaSq);
}
// 2018-03-28 21:04
TComplex Propagator::BW::getVal(const Double_t &s,
        const vector<Double_t> &paras) {
    TComplex ci(0, 1);
    TComplex pro = 1.0/(sqrt(s) - paras[1] - 0.5*ci * paras[2]);
    return pro;
}
// 2018-03-28 21:14
TComplex Propagator::K1430_0::getVal(const Double_t &s,
        const vector<Double_t> &paras) {
    return Propagator::Flatte::getVal(
            paras[0]*paras[0], s,
            paras[1], LauConstants::mKSq, LauConstants::mPiSq,
            paras[2], LauConstants::mPi0Sq, LauConstants::mEtaSq);
}
// 2018-03-30 16:12
// used to paramerize rho -> K K
// q is the break momenta depend on the s
Double_t Propagator::GS::h(const Double_t &m, const Double_t &q) {
    Double_t h(0.);
    h = 2 * q / (TMath::Pi() * m )
        * log((m+2*q) / (LauConstants::mK + LauConstants::mK0) );
    return h;
}
Double_t Propagator::GS::dh(const Double_t &m0, const Double_t& q0) {
    Double_t h = Propagator::GS::h(m0, q0);
    return h*(1.0/(8*q0*q0)-1.0/(2*m0*m0))+1.0/(2*LauConstants::pi*m0*m0);
}
Double_t Propagator::GS::f(const Double_t& s, const Double_t& m0,
        const Double_t &width, const Double_t &q0, const Double_t& _q) {
    Double_t mK = 0.5*(LauConstants::mK + LauConstants::mK0);
    Double_t tmp = 1;
    tmp *= width * m0*m0 /(q0*q0*q0);
    tmp *=  _q*_q * (Propagator::GS::h(sqrt(s), _q)
            - Propagator::GS::h(m0, q0) )
        + q0*q0 * (m0*m0 - mK*mK) * Propagator::GS::dh(m0, q0);
    return tmp;
}
Double_t Propagator::GS::d(const Double_t& mass, const Double_t &q0) {
    Double_t mK = 0.5*(LauConstants::mK + LauConstants::mK0);
    Double_t item1 = 3.0 / LauConstants::pi * mK * mK / (q0*q0) *
        log(0.5* (mass + 2*q0)/mK);
    Double_t item2 = 0.5 * mass/(LauConstants::pi*q0);
    Double_t item3 = mK*mK *mass /(LauConstants::pi * q0*q0*q0);
    return item1 + item2 - item3;
}
// spin, r, mass, width
TComplex Propagator::GS::getVal(const Double_t& s,
        const Double_t &Sa, const Double_t &Sb,
        const vector<Double_t> & paras) {
    // // cout<<__func__<<__LINE__<<endl;
    Int_t spin = paras[0];
    Double_t rRes = paras[1];
    Double_t mass = paras[2];
    Double_t width0 = paras[3];
    // // cout<<"spin: "<<spin<<endl;
    // // cout<<"rRes: "<<rRes<<endl;
    // // cout<<"mass: "<<mass<<endl;
    // // cout<<"width0: "<<width0<<endl;

    TComplex ci(0, 1);
    Double_t q = RooBarrier::BreakMomenta(s, Sa, Sb);
    Double_t Sr0 = mass*mass;
    Double_t q0 = RooBarrier::BreakMomenta(Sr0, Sa, Sb);
    Double_t width = Propagator::RBW::width(mass,  s,  Sa,  Sb, rRes,  spin);
    // cal pro
    Double_t upitem = 1 + Propagator::GS::d(mass, q0);
    Double_t donitem1 = mass*mass - s;
    Double_t donitem2 = Propagator::GS::f(mass, s, width0, q0, q);
    TComplex donitem3 = ci * mass *width0 * width;
    TComplex pro = upitem / (donitem1 + donitem2 - donitem3);
    // // cout<<"d: "<<upitem<<endl;
    // // cout<<"donitem1:"<<donitem1<<endl;
    // // cout<<"donitem2:"<<donitem2<<endl;
    // // cout<<"donitem3:"<<donitem3<<endl;
    return pro;
}
/*


TComplex RooPropogator::propagatorf600(Double_t mass, Double_t width, Double_t sx)const{
    TComplex ci(0,1);
    //----------form used by E791-------------------------------
    Double_t factor = sqrt(1-(4*mpi*mpi)/sx);
    TComplex prop = 1.0/(mass*mass-sx-ci*width*factor);
    return prop;
}
TComplex RooPropogator::rho4Pi(Double_t Sr)const{
    TComplex rho(0,0);
    TComplex ci(0,1);
    Double_t temp = 1-16*mpi*mpi/Sr;
    if(temp > 0) rho = TComplex::One()*sqrt(temp)/(1+exp(9.8-3.5*Sr));
    if(temp < 0) rho = ci*sqrt(-temp)/(1+exp(9.8-3.5*Sr));
    return rho;
}
TComplex RooPropogator::propagatorsigma500(Double_t Sr, Double_t Sa, Double_t Sb)const{
    //----------form Int_troduced by Zou and Bugg-----------------
    TComplex ci(0,1);
    Double_t f = 0.5843+1.6663*Sr;
    Double_t M = 0.9264;
    Double_t mpi2 = mpi*mpi;
    Double_t mass2 = M*M;
    Double_t g1 = f*(Sr-mpi2/2)/(mass2-mpi2/2)*exp((mass2-Sr)/1.082);
    TComplex rho1s = Flatte_rhoab(Sr,Sa,Sb);
    TComplex rho1M = Flatte_rhoab(mass2,Sa,Sb);
    TComplex rho2s = rho4Pi(Sr);
    TComplex rho2M = rho4Pi(mass2);
    TComplex prop = 1.0/(M*M-Sr-ci*M*(g1*rho1s/rho1M+0.0024*rho2s/rho2M));
    return prop;
}

TComplex RooPropogator::FlatteK1430_0(Double_t m0, Double_t *g, Double_t sx)const{
    // m0 is the nominal mass of resonance, i.e m0
    // but sx is the mass square: sx = m(R)*m(R)
    // for K(1430): there are two channels:
    // K(1430)0 -> K+ pi- and K(1430)0 -> K0 eta'
    // and Sa[0] = m(K)^2, Sb[0] = m(pi)^2
    // similarly Sa[1] = m(K0)^2, Sb[1] = m(eta')^2
    Double_t mKp = 0.493766;
    Double_t mK0 = 0.497611;
    Double_t mEtap =  0.95778;
    Double_t mPip = 0.13957061;
    Double_t Sa[2] = {mKp*mKp, mK0*mK0};
    Double_t Sb[2] = {mPip*mPip, mEtap*mEtap};
    TComplex prop = Flatte(m0, g, sx, Sa, Sb);
    return prop;
}
TComplex RooPropogator::FlatteK1430_p(Double_t m0, Double_t *g, Double_t sx)const{
    // for K(1430)+: there are two channels:
    // K(1430)+ -> K+ pi0
    // K(1430)+ -> K+ eta'
    Double_t mKp = 0.493766;
    Double_t mK0 = 0.497611;
    Double_t mEtap =  0.95778;
    Double_t mPi0 = 0.134977;
    Double_t Sa[2] = {mKp*mKp, mKp*mKp};
    Double_t Sb[2] = {mPi0*mPi0, mEtap*mEtap};
    TComplex prop = Flatte(m0, g, sx, Sa, Sb);
    return prop;
}
TComplex RooPropogator::Flattea980_p(Double_t m0, Double_t *g, Double_t sx)const{
    // for K(1430)+: there are two channels:
    // a(980)+ -> K+ K0
    // a(980)+ -> pi+ eta
    Double_t mKp = 0.493766;
    Double_t mK0 = 0.497611;
    Double_t mEta = 0.547862;
    Double_t mPip = 0.13957061;
    Double_t Sa[2] = {mKp*mKp, mPip * mPip};
    Double_t Sb[2] = {mK0*mK0, mEta*mEta};
    TComplex prop = Flatte(m0, g, sx, Sa, Sb);
    return prop;
}
TComplex RooPropogator::Flattea980_0(Double_t m0, Double_t *g, Double_t sx)const{
    // for K(1430)+: there are two channels:
    // a(980)0 -> K+ K-
    // a(980)0 -> pi0 eta
    Double_t mKp = 0.493766;
    Double_t mEta = 0.547862;
    Double_t mPi0 = 0.134977;
    Double_t Sa[2] = {mKp*mKp, mPi0 * mPi0};
    Double_t Sb[2] = {mKp*mKp, mEta*mEta};
    TComplex prop = Flatte(m0, g, sx, Sa, Sb);
    return prop;
}
*/
// 13:18 2018-06-12
TComplex Propagator::a980_p_3C::getVal(const Double_t &s,
        const vector<Double_t> & paras) {
    return Propagator::Flatte::getVal(paras[2]*paras[2], s,
            paras[3], LauConstants::mPi0Sq, LauConstants::mEtaSq,
            paras[4], LauConstants::mKSq, LauConstants::mKSq,
            paras[5], LauConstants::mPi0Sq, LauConstants::mEtaPrimeSq);
}
