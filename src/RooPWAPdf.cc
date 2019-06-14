// Copyright (c) 2019-3-3 maxx
#include <fstream>
#include "RooPWAPdf.hh"
#include "RooBarrier.hh"
#include "RooSpinFactor.hh"
#include "RooP4Vector.hh"
#include "RooPropogator.hh"
#include "LineShape.hh"

#include "RooArgList.h"
#include "TTree.h"
#include "TRandom.h"
#include "TGenPhaseSpace.h"
#include "TSystem.h"

using std::cout;
using std::endl;
using std::ios;

ClassImp(RooPWAPdf)
    RooPWAPdf::RooPWAPdf(const char *name, const char *title,
            RooAbsReal& _p11,  RooAbsReal& _p12,   RooAbsReal& _p13,
            RooAbsReal& _p14, RooAbsReal& _p21,  RooAbsReal& _p22,
            RooAbsReal& _p23,   RooAbsReal& _p24, RooAbsReal& _p31,
            RooAbsReal& _p32,   RooAbsReal& _p33,   RooAbsReal& _p34,
            const TString& PHSPdat):
        RooAbsPdf(name, title),
        p11("p11", "p11", this, _p11),
        p12("p12", "p12", this, _p12),
        p13("p13", "p13", this, _p13),
        p14("p14", "p14", this, _p14),
        p21("p21", "p21", this, _p21),
        p22("p22", "p22", this, _p22),
        p23("p23", "p23", this, _p23),
        p24("p24", "p24", this, _p24),
        p31("p31", "p31", this, _p31),
        p32("p32", "p32", this, _p32),
        p33("p33", "p33", this, _p33),
        p34("p34", "p34", this, _p34),
        _ParameterCol("ParameterCol", "", this),
        _OthersCol("OthersCol", "", this),
        _rhoList("rhoList", "", this),
        _phiList("phiList", "", this) {
    _parsItr = _ParameterCol.createIterator();
    _others = _OthersCol.createIterator();
    _rhoItr = _rhoList.createIterator();
    _phiItr = _phiList.createIterator();
    // _rhoList = RooArgList("rhoList");
    // _phiList = RooArgList("phiList");
    _PHSPDat = PHSPdat;
    _SoloNotIdenPar = 0;
    initialize();
    _paraSize = new Int_t[50];
    for (Int_t i  = 0; i< 50; i++) {
        _paraSize[i] = 0;
    }
    _parameters = new Double_t*[50];
    for (Int_t i = 0; i < 50; ++i) {
      _parameters[i] = new Double_t[10];
    }
    _cofV = new TComplex[50];
    for (Int_t index = 0; index < 50; index++) {
        _cofV[index] = TComplex(0);
    }
    m_fixbr = false;
    _freeShape = false;
    m_fixFFa0 = 0.266;
}
RooPWAPdf::RooPWAPdf(const RooPWAPdf& other, const char* name):
    RooAbsPdf(other, name),
    p11("p11", this, other.p11),
    p12("p12", this, other.p12),
    p13("p13", this, other.p13),
    p14("p14", this, other.p14),
    p21("p21", this, other.p21),
    p22("p22", this, other.p22),
    p23("p23", this, other.p23),
    p24("p24", this, other.p24),
    p31("p31", this, other.p31),
    p32("p32", this, other.p32),
    p33("p33", this, other.p33),
    p34("p34", this, other.p34),
    _ParameterCol("ParameterCol", this, other._ParameterCol),
    _OthersCol("OthersCol", this, other._OthersCol),
    _rhoList("rhoList", this, other._rhoList),
    _phiList("phiList", this, other._phiList) {
    _parsItr = _ParameterCol.createIterator();
    _others  = _OthersCol.createIterator();
    _phiItr  = _phiList.createIterator();
    _rhoItr  = _rhoList.createIterator();
    _rD = other._rD;
    _fracDat = other._fracDat;
    _PHSPDat = other._PHSPDat;
    _SoloNotIdenPar = other._SoloNotIdenPar;
    // _ParaList = other._ParaList;
    _NameList = other._NameList;
    _LineShapeList = other._LineShapeList;
    _DecayNumList = other._DecayNumList;
    // _NeedSpin = other._NeedSpin;
    // _NeedrRes = other._NeedrRes;
    /*
       _mcp1 = other._mcp1;
       _mcp2 = other._mcp2;
       _mcp3 = other._mcp3;
       _weight = other._weight;
       */


    m_fixbr = other.m_fixbr;
    _freeShape = other._freeShape;
    m_fixFFa0 = other.m_fixFFa0;
    Nmc = other.Nmc;
    _mcp1 = new Double_t*[Nmc];
    _mcp2 = new Double_t*[Nmc];
    _mcp3 = new Double_t*[Nmc];
    _weight = new Double_t[Nmc];
    _Couple = new TComplex*[50];
    _ParAmpInt = new TComplex*[50];
    for (int i = 0; i < 50; i++) {
        _Couple[i] = new TComplex[50];
        _ParAmpInt[i] = new TComplex[50];
    }
    for (Int_t i = 0; i < Nmc; i++) {
        _mcp1[i] = new Double_t[4];
        _mcp2[i] = new Double_t[4];
        _mcp3[i] = new Double_t[4];
    }
    for (Int_t i = 0; i < Nmc; i++) {
        _mcp1[i][0] = other._mcp1[i][0];
        _mcp1[i][1] = other._mcp1[i][1];
        _mcp1[i][2] = other._mcp1[i][2];
        _mcp1[i][3] = other._mcp1[i][3];
        _mcp2[i][0] = other._mcp2[i][0];
        _mcp2[i][1] = other._mcp2[i][1];
        _mcp2[i][2] = other._mcp2[i][2];
        _mcp2[i][3] = other._mcp2[i][3];
        _mcp3[i][0] = other._mcp3[i][0];
        _mcp3[i][1] = other._mcp3[i][1];
        _mcp3[i][2] = other._mcp3[i][2];
        _mcp3[i][3] = other._mcp3[i][3];
        _weight[i] = other._weight[i];
    }
    for (int i  = 0; i < 50; i++) {
        for (int j  = 0; j < 50; j++) {
            _Couple[i][j] = other._Couple[i][j];
            _ParAmpInt[i][j] = other._ParAmpInt[i][j];
        }
    }
    // // cout<<__func__<<" cpoy done"<<endl;
    _paraSize = new Int_t[50];
    for (Int_t i  = 0; i< 50; i++) {
        _paraSize[i] = other._paraSize[i];
    }
    _parameters = new Double_t*[50];
    for (Int_t i = 0; i < 50; ++i) {
      _parameters[i] = new Double_t[10];
    }
    for (Int_t i = 0; i < 50; ++i) {
        for (Int_t j = 0; j < 10; ++j) {
            _parameters[i][j] = other._parameters[i][j];
        }
    }
    _cofV = new TComplex[50];
    for (Int_t index = 0; index < 50; index++) {
        _cofV[index] = TComplex(0);
    }
}
void RooPWAPdf::initialize() {
    Nmc = 1E6;
    _mcp1 = new Double_t*[Nmc];
    _mcp2 = new Double_t*[Nmc];
    _mcp3 = new Double_t*[Nmc];
    _weight = new Double_t[Nmc];
    _Couple = new TComplex*[50];
    _ParAmpInt = new TComplex*[50];
    for (int i = 0; i < 50; i++) {
        _Couple[i] = new TComplex[50];
        _ParAmpInt[i] = new TComplex[50];
    }
    for (Int_t i = 0; i < Nmc; i++) {
        _mcp1[i] = new Double_t[4];
        _mcp2[i] = new Double_t[4];
        _mcp3[i] = new Double_t[4];
    }
    Double_t fx1, fy1, fz1, ft1, fx2, fy2, fz2, ft2, fx3, fy3, fz3, ft3;
    Double_t weight;
    FILE *fp;
    if ((fp=fopen(_PHSPDat, "r")) == NULL) {
        printf("can't open input file");
        exit(0);
    }
    Int_t i = 0;
    while (fscanf (fp, "%lf%lf%lf%lf\n%lf%lf%lf%lf\n%lf%lf%lf%lf%lf\n",
                &ft1 , &fx1 , &fy1 , &fz1 ,
                &ft2 , &fx2 , &fy2 , &fz2 ,
                &ft3 , &fx3 , &fy3 , &fz3 , &weight) != EOF) {
        if ( i >= Nmc) break;
        _mcp1[i][0] = ft1;
        _mcp1[i][1] = fx1;
        _mcp1[i][2] = fy1;
        _mcp1[i][3] = fz1;
        _mcp2[i][0] = ft2;
        _mcp2[i][1] = fx2;
        _mcp2[i][2] = fy2;
        _mcp2[i][3] = fz2;
        _mcp3[i][0] = ft3;
        _mcp3[i][1] = fx3;
        _mcp3[i][2] = fy3;
        _mcp3[i][3] = fz3;
        _weight[i] = weight;
        i++;
    }
    fclose(fp);
    Nmc = i;

    _parameters = new Double_t*[50];
    for (Int_t i = 0; i < 50; ++i) {
      _parameters[i] = new Double_t[10];
    }
    for (Int_t i = 0; i < 50; ++i) {
        for (Int_t j = 0; j < 10; ++j) {
            _parameters[i][j] = 0;
        }
    }
}


Double_t RooPWAPdf::evaluate() const {
    Double_t p4_1[4] = {p11, p12, p13, p14};
    Double_t p4_2[4] = {p21, p22, p23, p24};
    Double_t p4_3[4] = {p31, p32, p33, p34};
    FillParameter();
    return calEva(p4_1, p4_2, p4_3)*ConstrFactor(m_fixFFa0);
}

Double_t RooPWAPdf::evaluate(Double_t _p11, Double_t _p12, Double_t
        _p13, Double_t _p14, Double_t _p21, Double_t _p22, Double_t
        _p23, Double_t _p24, Double_t _p31, Double_t _p32, Double_t
        _p33, Double_t _p34) const {
    Double_t p4_1[4] = {_p11, _p12, _p13, _p14};
    Double_t p4_2[4] = {_p21, _p22, _p23, _p24};
    Double_t p4_3[4] = {_p31, _p32, _p33, _p34};
    FillParameter();
    return calEva(p4_1, p4_2, p4_3);
}

Int_t RooPWAPdf::getAnalyticalIntegral(RooArgSet& allVars,
       RooArgSet& analVars, const char* rangeName) const {
    RooArgSet theSet1, theSet2;
    theSet1.add(RooArgSet(p11.arg(), p12.arg(), p13.arg(), p14.arg(),
                p21.arg(), p22.arg(), p23.arg(), p24.arg()));
    theSet2.add(RooArgSet(p31.arg(), p32.arg(), p33.arg(), p34.arg()));
    RooArgSet theSet(theSet1, theSet2, " ");
    if (matchArgs(allVars, analVars, theSet)) {
        return 1;
    }
    return 0;
}

Double_t RooPWAPdf::analyticalIntegral(Int_t code, const char* rangeName)
    const {
    assert(code == 1);
    Double_t sum = 0;
    // for(Int_t i=0;i<Nmc;i++){
    //     Double_t eva = calEva(_mcp1[i], _mcp2[i], _mcp3[i]);
    //     sum = sum + eva*_weight[i];
    //  }
    //
    Int_t nAmps = _NameList.size();
    if (_freeShape) {
        for (int nn = 0; nn < nAmps; ++nn) {
            for (Int_t mm = 0; mm < nAmps ; ++mm) {
                _ParAmpInt[nn][mm] = 0;
            }
        }
        TComplex *partiAmpList = new TComplex[nAmps];
        for (Int_t i  = 0; i < Nmc; ++i) {
            Double_t *p4_1 = _mcp1[i];
            Double_t *p4_2 = _mcp2[i];
            Double_t *p4_3 = _mcp3[i];
            Double_t weight = _weight[i];
            for (Int_t index  = 0; index < nAmps; ++index) {
                partiAmpList[index] = getFFVal(index, p4_1, p4_2, p4_3);
            }

            for (int nn = 0; nn < nAmps; ++nn) {
                for (Int_t mm = 0; mm < nAmps ; ++mm) {
                    _ParAmpInt[nn][mm]  += 1.0/Nmc * weight
                        * partiAmpList[nn]
                        * TComplex::Conjugate(partiAmpList[mm]);
                }
            }
        }
    }
    TComplex amptot(0);
    for (int i = 0; i < nAmps; ++i) {
        for (int j = 0; j < nAmps; ++j) {
            amptot += _cofV[i]
                *TComplex::Conjugate(_cofV[j])*_ParAmpInt[i][j];
        }
    }
    return amptot.Rho();
}
void RooPWAPdf::setFracDat(const TString &dat) {
    _fracDat = dat;
}
void RooPWAPdf::project(const char* fname) {
    FillParameter();
    TFile f(fname, "recreate");

    TTree t("project", "the project of dalitz");

    Double_t m12, m13, m23, eva;
    Int_t nAmps = _NameList.size();
    const Int_t length = nAmps*nAmps;
    Double_t *FF = new Double_t[length];

    Double_t m_p4_1[4], m_p4_2[4], m_p4_3[4];
    t.Branch("eva", &eva, "eva/D");
    t.Branch("m12", &m12, "m12/D");
    t.Branch("m13", &m13, "m13/D");
    t.Branch("m23", &m23, "m23/D");
    t.Branch("p1", m_p4_1, "p1[4]/D");
    t.Branch("p2", m_p4_2, "p2[4]/D");
    t.Branch("p3", m_p4_3, "p3[4]/D");
    Int_t waves = length;
    t.Branch("length", &waves, "length/I");
    t.Branch("FF", FF, "FF[length]/D");

    for (Int_t i = 0;  i< Nmc; i++) {
        TLorentzVector p4_1(_mcp1[i][1], _mcp1[i][2], _mcp1[i][3], _mcp1[i][0]);
        TLorentzVector p4_2(_mcp2[i][1], _mcp2[i][2], _mcp2[i][3], _mcp2[i][0]);
        TLorentzVector p4_3(_mcp3[i][1], _mcp3[i][2], _mcp3[i][3], _mcp3[i][0]);

        // std::////////////////cout<<"mass:\t"<<p4Ks.M()<<"\t"<<p4Pion.M()<<"\t"<<p4Eta.M()<<std::endl;
        m12 = (p4_1 + p4_2).M();
        m13 = (p4_1 + p4_3).M();
        m23 = (p4_2 + p4_3).M();
        for (int j = 0; j < 4; ++j) {
            m_p4_1[j] = _mcp1[i][j];
            m_p4_2[j] = _mcp2[i][j];
            m_p4_3[j] = _mcp3[i][j];
        }

        // calculate the weight according to PDF
        // // // // // // // cout<<__func__<<endl;
        eva =  _weight[i] * calEva(_mcp1[i],  _mcp2[i],  _mcp3[i]);
        waves = nAmps * nAmps;
        //  //cout<<__func__<<endl;
        //  //cout<<"eva: "<<eva<<endl;
        //  TComplex totAmp(0);
        //  for(Int_t ii =0; ii<nAmps; ++ii){
        //      totAmp += _cofV[ii]
        //          * getFFVal(ii,, _mcp1[i], _mcp2[i], _mcp3[i]);
        //  }
        // // cout<<"totAmp:" <<totAmp<<" rh02 = "<< totAmp.Rho2()<<endl;

        for (Int_t ii = 0; ii < nAmps; ++ii) {
            for (Int_t jj = 0; jj < nAmps; ++jj) {
                Int_t index = ii*nAmps + jj;
                TComplex ampii = _cofV[ii]
                    * getFFVal(ii, _mcp1[i], _mcp2[i], _mcp3[i]);
                TComplex ampjj = _cofV[jj]
                    * getFFVal(jj, _mcp1[i], _mcp2[i], _mcp3[i]);
                FF[index] = _weight[i] * 0.5
                    *(ampii * TComplex::Conjugate(ampjj)
                    + ampjj * TComplex::Conjugate(ampii) );
            }
        }
        t.Fill();
    }
    t.Write();
    f.Close();
    delete[] FF;
}

void RooPWAPdf::fitFractions(const RooArgList& newParList,
        Bool_t prInt_t, ostream& os) {
    // recover from newPar
    Int_t nAmps = _NameList.size();
    TIterator *itrnewPar = newParList.createIterator();
    RooRealVar *newVar(0);
    RooRealVar *aRho(0), *aPhi(0);
    Double_t *oldrhoList = new Double_t[nAmps];
    TIterator *itrRho = _rhoList.createIterator();
    TIterator *itrPhi = _phiList.createIterator();
    itrRho->Reset();
    itrPhi->Reset();
    for (Int_t j = 0; j < nAmps; j++) {
        aRho = reinterpret_cast<RooRealVar*>(itrRho->Next());
        oldrhoList[j] = aRho->getVal();
    }
    for (Int_t i = 0; i< newParList.getSize(); i++) {
        newVar = reinterpret_cast<RooRealVar*>(itrnewPar->Next());
        // // cout<<"The "<<newVar->GetName()<<" is Reset"<<endl;
        itrRho->Reset();
        itrPhi->Reset();
        for (Int_t j = 0; j < nAmps; j++) {
            aRho = reinterpret_cast<RooRealVar*>(itrRho->Next());
            //   //cout<<"The "<<newVar->GetName()<<" need Reset"<<endl;
            //   //cout<<"tmp "<<aRho->GetName()<<" "<<aRho->getVal()<<endl;
            if (TString(newVar->GetName()) == TString(aRho->GetName())) {
                // //cout<<"fitFractions:: "<<endl;
                // //cout<<aRho->GetName()<<" old value"
                //     <<aRho->getVal()<<endl;
                aRho->setVal(newVar->getVal() );
                // // cout<<aRho->GetName()<<" new value"
                //    <<aRho->getVal()<<endl;
            }
            aPhi = reinterpret_cast<RooRealVar*>(itrPhi->Next());
            if (TString(newVar->GetName()) == TString(aPhi->GetName())) {
                aPhi->setVal(newVar->getVal() );
            }
        }
    }
    FillParameter();
    TComplex totamp(0);
    for (int i = 0; i< nAmps; ++i) {
        for (int j = 0; j < nAmps; ++j) {
            totamp += _Couple[i][j] * _cofV[i]
                *TComplex::Conjugate(_cofV[j]);
        }
    }


    Double_t *parAmps = new Double_t[nAmps];
    for (int i = 0; i < nAmps; ++i) {
        parAmps[i]  = _Couple[i][i].Rho() * _cofV[i].Rho2();
        // // cout<<"parAmps["<<i<<"]"<<parAmps[i]<<endl;
    }

    Double_t total = 0;
    for (Int_t i = 0; i < nAmps; i++) {
        const TString &name = _NameList[i];
        Double_t FF = parAmps[i]/totamp.Rho();
        os<< "Fraction:: " << name << " BF: "
            <<FF << endl;
        total += FF;
    }
    os << "Fraction:: total BF: " << total << endl;
    // reset to teh inital value
    itrRho->Reset();
    for (Int_t j = 0; j < nAmps; j++) {
        aRho = reinterpret_cast<RooRealVar*>(itrRho->Next());
        aRho->setVal(oldrhoList[j]);
        // // cout<<aRho->GetName()<<" old value:"<<oldrhoList[j]<<endl;
    }
    delete itrRho;
    delete itrnewPar;
    delete[] oldrhoList;
}
Int_t RooPWAPdf::setPar(const RooArgList& newPar) {
    return 0;
}
bool RooPWAPdf::addResonance(const TString &name,
        const LineShape::Shape &lineShape,
        const DecayType::DecayType &modetype,
        RooRealVar &rho, RooRealVar &phi,
        RooArgList &params) {
    //  RooRealVar *arho = new RooRealVar ("Rho"+name, "", 1, 0, 20);
    //  RooRealVar *aphi = new RooRealVar ("Phi"+name, "", 0, -15, 15);
    if (params.getSize() != Propagator::GetParaSize(lineShape)) {
        std::cerr << "Error: Need " << Propagator::GetParaSize(lineShape)
             <<  " parameters for lineShape of " << name
             << ", but " << params.getSize() << " provided!!!"
            <<endl;
        exit(0);
    }
    // fix the first amplitude rho and phi
    if (_rhoList.getSize() == 0) {
        rho.setConstant();
        phi.setConstant();
    }
    // _ParameterCol.add(*arho);
    // _ParameterCol.add(*aphi);

    if (std::find (_NameList.begin (), _NameList.end (), name)
            != _NameList.end()) {
        std::cerr << "Error:: The resonance : " << name
            <<" is already exist" <<std::endl;
        exit(0);
        return false;
    } else {
        // std::// // // // // // cout<<"Inf:: addResonance,"
        //   <<" name:"<< name
        //  <<" LineShape:"<<lineShape
        //  <<" modetype: "<<modetype<<endl;
    }
    _rhoList.add(rho);
    _phiList.add(phi);
    _ParameterCol.add(params);
    // _ParaList.push_back( params );
    _LineShapeList.push_back(lineShape);
    _NameList.push_back(name);
    _DecayNumList.push_back(modetype);
    RooRealVar *spin = reinterpret_cast<RooRealVar*>(params.at(0));
    // delete arho;
    // delete aphi;
    cout << "Inf:: addResonance," << endl;
    cout << "name: "             << name           << endl;
    cout << "LineShape: "        << lineShape      << endl;
    cout << "modetype: "         << modetype       << endl;
    cout << "Spin: "              << spin->getVal() << endl;
    return true;
}
void RooPWAPdf::setrD(const Double_t &rD) {
    _rD = rD;
}
void RooPWAPdf::setIdenticalParticle(const Int_t& mode) {
    if (mode == 12 || mode == 21) {
        _SoloNotIdenPar = 3;
    } else if (mode == 23 || mode == 32) {
        _SoloNotIdenPar = 1;
    } else if (mode == 13 || mode == 31) {
        _SoloNotIdenPar = 2;
    } else {
        _SoloNotIdenPar = 0;
    }
}

/*
 * the parWname is the name of partial wave. The name is set at the
 * function addResonance(name...). The name is used to get the
 * parameter of this resonance, which is a "RooArgList".
 * The parameters are stored in a List by order, exactly:
 * spin, radius, mass, ...
 * RBW or BW | --- width
 * Flatte    | --- g1, m1_1, m2_1, g2, m1_2, m2_2,sA
 * a0+       | --- gKK, gPiEta
 */
TComplex RooPWAPdf::PartialAmplitude(const Int_t & index,
        const vector<Double_t> &paras,
        const Double_t p4_1[4], const Double_t p4_2[4],
        const Double_t p4_solo[4]) const {
    Int_t Ang = paras[0];
    Double_t rRes = 3.0;
    if (paras.size() > 1) {
        rRes = paras[1];
    }
    // cout<<__func__<<endl;

    TComplex _pro(1, 0);
    Double_t Sa = RooP4Vector::dot(p4_1, p4_1);
    // cout<<"Sa:"<<Sa<<endl;
    Double_t Sb = RooP4Vector::dot(p4_2, p4_2);
    // cout<<"Sb: "<<Sb<<endl;
    Double_t s = 2*RooP4Vector::dot(p4_1, p4_2) + Sa + Sb;
    LineShape::Shape lineShape = _LineShapeList[index];
    _pro = Propagator::getVal(lineShape, s, Sa, Sb, paras);
    // cout<<"S: "<<s<<endl;

    Double_t p4R[4] = { p4_1[0]+ p4_2[0], p4_1[1]+p4_2[1],
        p4_1[2]+p4_2[2], p4_1[3]+p4_2[3] };
    Double_t barrFactorD = RooBarrier::Factor(Ang, p4R, p4_solo, _rD);
    Double_t barrFactorRes = RooBarrier::Factor(Ang, p4_1, p4_2, rRes);
    Double_t spinFactor = RooSpinFactor::SpinFactor(Ang, p4_1, p4_2,
            p4_solo);
    // cout<<"spinFactor: "<<spinFactor<<endl;
    // cout<<"rRes: "<<rRes<<endl;
    // cout<<"pro: "<<_pro<<endl;
    // cout<<"barrFactorRes: "<<barrFactorRes<<endl;
    // cout<<"barrFactorD: "<<barrFactorD<<endl;
    return spinFactor * _pro * barrFactorRes * barrFactorD;
}
void RooPWAPdf::init() {
    FillParameter();
    Int_t nAmps = _NameList.size();
    for (int i  = 0 ; i < nAmps; ++i) {
        Double_t *pars = _parameters[i];
        // cout<<"waves: "<<i<<endl;
        for (int j  = 0; j< _paraSize[i] ; j++) {
            // cout<<" par: "<<pars[j]<<endl;
        }
    }
    // cout<<"_SoloNotIdenPar: "<<_SoloNotIdenPar<<endl;
    PartAmpInt();
    // Double_t
    Int_t i = 0;
    Double_t _sum = 0;
    for (int nn = 0; nn < nAmps; ++nn) {
        for (Int_t mm = 0; mm < nAmps ; ++mm) {
            _ParAmpInt[nn][mm] = 0;
        }
    }
    TComplex *partiAmpList = new TComplex[nAmps];
    for (Int_t i  = 0; i < Nmc; ++i) {
        Double_t *p4_1 = _mcp1[i];
        Double_t *p4_2 = _mcp2[i];
        Double_t *p4_3 = _mcp3[i];
        Double_t weight = _weight[i];
        for (Int_t index  = 0; index < nAmps; ++index) {
            partiAmpList[index] = getFFVal(index, p4_1, p4_2, p4_3);
        }

        for (int nn = 0; nn < nAmps; ++nn) {
            for (Int_t mm = 0; mm < nAmps ; ++mm) {
                _ParAmpInt[nn][mm]  += 1.0/Nmc * weight
                    * partiAmpList[nn]
                    * TComplex::Conjugate(partiAmpList[mm]);
            }
        }
    }

    for (int nn = 0; nn < nAmps; ++nn) {
        for (Int_t mm = 0; mm < nAmps ; ++mm) {
            cout << "ParAmpInt :" << nn << " and " << mm
                 << " is " << _ParAmpInt[nn][mm] << endl;
        }
    }
    delete[] partiAmpList;
}
void RooPWAPdf::test() {
    // for(
    //         vector<RooArgList>::iterator itr = _ParaList.begin(),
    //         end = _ParaList.end();
    //
    //         itr != end;
    //         itr++){
    //     //////////////cout<<"resonance name"<<itr->first<<endl;
    //     TIterator *parLitr = (*itr).createIterator();
    //     RooRealVar *par(0);
    //     while( 0 !=( par= (RooRealVar*) parLitr->Next())){
    //         ////////////cout<<" parameter: "<<par->GetName()
    //          //  <<" val: "<<par->getVal()
    //          //  <<endl;
    //     }
    //     delete parLitr;
    // }
    // // // // // // // // cout<<"TEST calEva:\t"<<eva<<endl;
    MCIntG();
}
Double_t RooPWAPdf::calEva(const Double_t p4_1[4],
        const Double_t p4_2[4], const Double_t p4_3[4]) const {
    Int_t nAmps = _NameList.size();
    FillParameter();
    TComplex _Amplitude(0);
    // cout<<__func__<<endl;
    for (int index  = 0; index < nAmps; ++index) {
        // cout<<"index: "<<index<<endl;
        _Amplitude += _cofV[index] * getFFVal(index, p4_1, p4_2, p4_3);
        // cout<<"rho: "<<_cofV[index]<<endl;
        // cout<<"getFFVal: "<<getFFVal(index, p4_1, p4_2, p4_3)<<endl;
    }
    Double_t value = _Amplitude.Rho2();
    return (value <= 0) ? 1e-20 : value;
}
Double_t RooPWAPdf::MCIntG() {
    ifstream f;
    f.open(_fracDat);
    Int_t i = 0;
    Double_t _sum = 0;
    while (1) {
        if (f.eof () != 0) break;
        Double_t p4_1[4], p4_2[4], p4_3[4], weight;
        f >> p4_1[0] >> p4_1[1] >> p4_1[2] >> p4_1[3];
        f >> p4_2[0] >> p4_2[1] >> p4_2[2] >> p4_2[3];
        f >> p4_3[0] >> p4_3[1] >> p4_3[2] >> p4_3[3];
        f>> weight;
        _sum += weight * calEva(p4_1, p4_2, p4_3);
        i += 1;
    }
    return _sum/i;
}
// 2018-03-27
void RooPWAPdf::setPHSPDat(const TString &dat) {
    _PHSPDat = dat;
}

// 2018-03-28 21:16
Double_t RooPWAPdf::getFFVal(const Int_t &ii, const Int_t &jj,
        const  Double_t p4_1[4], const Double_t p4_2[4],
        const Double_t p4_3[4]) {
    TComplex _Ampii = getFFVal(ii, p4_1, p4_2, p4_3);
    TComplex _Ampjj = getFFVal(jj, p4_1, p4_2, p4_3);
    TComplex module = _Ampii*TComplex::Conjugate(_Ampjj)
        + _Ampjj * TComplex::Conjugate(_Ampii);
    return 0.5*module.Re();
}
// 2018-03-30 20:05
inline  RooPWAPdf::~RooPWAPdf() {
    for (Int_t i  = 0 ; i < Nmc;  i++) {
        delete[] _mcp1[i];
        delete[] _mcp2[i];
        delete[] _mcp3[i];
    }
    delete[] _weight;
    for (int i  = 0 ; i < 50 ; ++i) {
    delete[] _Couple[i];
    delete[] _ParAmpInt[i];
    }
    TIterator *itrRho = _rhoList.createIterator();
    TIterator *itrPhi = _phiList.createIterator();
    RooRealVar *aRho(0);
    RooRealVar *aPhi(0);
    // // // // // // // cout<<"_rhoList size: "<<_rhoList.getSize()<<endl;
    Int_t nAmps = _rhoList.getSize();
    for (Int_t index  = 0 ; index < nAmps; ++index) {
        TComplex cof, partAmp(0);
        aRho = reinterpret_cast<RooRealVar*>(itrRho->Next());
        aPhi = reinterpret_cast<RooRealVar*>(itrPhi->Next());
    }
}
// 2018-04-25 18:07
Double_t RooPWAPdf::ConstrFactor(const Double_t &FFa0) const {
    if (!m_fixbr) return 1.0;
    _rhoItr->Reset();
    _phiItr->Reset();
    int nAmps = _NameList.size();
    TComplex *RhoList = new TComplex[nAmps];
    for (Int_t index = 0; index < nAmps; ++index) {
        RooRealVar *aRho = reinterpret_cast<RooRealVar*>(_rhoItr->Next());
        RooRealVar *aPhi = reinterpret_cast<RooRealVar*>(_phiItr->Next());
        Double_t rho = aRho->getVal();
        Double_t phi = aPhi->getVal();
        RhoList[index] = TComplex(rho*TMath::Cos(phi), rho*TMath::Sin(phi));
    }
    TComplex donmin(0, 0);
    for (Int_t ii = 0; ii < nAmps; ++ii) {
        for (Int_t jj = 0; jj< nAmps; ++jj) {
            donmin += _Couple[ii][jj] * RhoList[ii] *
                TComplex::Conjugate(RhoList[jj]);
        }
    }
    // // cout<<"donminer:"<<donmin<<endl;
    Double_t AmpA0 = _Couple[0][0] * RhoList[0] *
        TComplex::Conjugate(RhoList[0]);
    Double_t g = fabs(AmpA0/(Double_t)donmin - FFa0);
    // // cout<<" g "<<g<<endl;
    _others->Reset();
    RooRealVar *aLm = reinterpret_cast<RooRealVar*>(_others->Next());
    Double_t aLmVal = aLm->getVal();
    delete[] RhoList;
    return TMath::Exp(-aLmVal * g  - 100 * g * g);
}
TComplex RooPWAPdf::PartAmp(const Int_t &index,
        const vector<Double_t> &paraLst,
        const Double_t p4_1[4],
        const Double_t p4_2[4],
        const Double_t p4_3[4]) const {
    TString name = _NameList[index];
    DecayType::DecayType modeNum = _DecayNumList[index];
    TComplex partAmp;
    if (modeNum == 1) {
        partAmp = PartialAmplitude(index, paraLst,
                p4_2, p4_3, p4_1);
        if (_SoloNotIdenPar == 0) {
            return partAmp;
        }
        if (_SoloNotIdenPar == 1) {
            partAmp += PartialAmplitude(index, paraLst,
                    p4_3, p4_2, p4_1);
            return partAmp;
        }
        if (_SoloNotIdenPar == 2) {
            partAmp += PartialAmplitude(index, paraLst,
                    p4_2, p4_1, p4_3);
            return partAmp;
        }
        if (_SoloNotIdenPar == 3) {
            partAmp += PartialAmplitude(index, paraLst, p4_1, p4_3, p4_2);
            return partAmp;
        }
    } else if (modeNum == 2) {
        partAmp = PartialAmplitude(index, paraLst, p4_3, p4_1, p4_2);
        switch (_SoloNotIdenPar) {
            case 1:
                partAmp += PartialAmplitude(index, paraLst, p4_2, p4_1, p4_3);
                break;
            case 2:
                partAmp += PartialAmplitude(index, paraLst, p4_3, p4_1, p4_2);
                break;
            case 3:
                partAmp += PartialAmplitude(index, paraLst, p4_3, p4_2, p4_1);
                break;
            default:
                break;
        }
    } else if (modeNum == 3) {
        partAmp = PartialAmplitude(index, paraLst, p4_1, p4_2, p4_3);
        switch (_SoloNotIdenPar) {
            case 1:
                partAmp += PartialAmplitude(index, paraLst, p4_1, p4_3, p4_2);
                break;
            case 2:
                partAmp += PartialAmplitude(index, paraLst, p4_3, p4_2, p4_1);
                break;
            case 3:
                partAmp += PartialAmplitude(index, paraLst, p4_2, p4_1, p4_3);
                break;
            default:
                break;
        }
    }
    return partAmp;
}
void RooPWAPdf::PartAmpInt() {
    // Double_t
    int nAmps = _NameList.size();
    RooRealVar *aPara(0);
    Int_t i = 0;
    Double_t _sum = 0;
    FillParameter();
    for (int nn = 0; nn < nAmps; ++nn) {
        for (Int_t mm = 0; mm < nAmps ; ++mm) {
            _Couple[nn][mm] = 0;
        }
    }
    TComplex *partiAmpList = new TComplex[nAmps];
    ifstream f;
    f.open(_fracDat);
    while (1) {
        if (f.eof () != 0) break;
        Double_t p4_1[4], p4_2[4], p4_3[4], weight;
        f >> p4_1[0] >> p4_1[1] >> p4_1[2] >> p4_1[3];
        f >> p4_2[0] >> p4_2[1] >> p4_2[2] >> p4_2[3];
        f >> p4_3[0] >> p4_3[1] >> p4_3[2] >> p4_3[3];
        f>> weight;
        i += 1;
        if (i > 4E5) break;
        _parsItr->Reset();
        for (Int_t index  = 0; index < nAmps; ++index) {
            partiAmpList[index] = getFFVal(index, p4_1, p4_2, p4_3);
        }
        for (int nn = 0; nn < nAmps; ++nn) {
            for (Int_t mm = 0; mm < nAmps ; ++mm) {
                _Couple[nn][mm]  += 1.0/Nmc * weight * partiAmpList[nn]
                    *TComplex::Conjugate(partiAmpList[mm]);
            }
        }
    }
    for (int nn = 0; nn < nAmps; ++nn) {
        for (Int_t mm = 0; mm < nAmps ; ++mm) {
            cout << "Couple of " << nn << " and " << mm << " is "
                 << _Couple[nn][mm] << endl;
        }
    }
    delete[] partiAmpList;
}
void RooPWAPdf::setLambda(RooRealVar &aLambda, const Double_t &FF) {
    m_fixbr = true;
    m_fixFFa0 = FF;
    _OthersCol.add(aLambda);
}
void RooPWAPdf::FreeLineShape() {
    _freeShape = true;
}
void RooPWAPdf::FillParameter() const {
    Int_t nAmps = _NameList.size();
    _parsItr->Reset();

    for (Int_t index  = 0 ; index < nAmps; ++index) {
        // get parameters
        RooRealVar *aPara(0);
        LineShape::Shape lineshape = _LineShapeList[index];
        Int_t size   = Propagator::GetParaSize(lineshape);
        _paraSize[index] = size;
        for (Int_t i  = 0; i < size; i++) {
            aPara = reinterpret_cast<RooRealVar*>(_parsItr->Next());
            _parameters[index][i] = aPara->getVal();
        }
    }

    _rhoItr->Reset();
    _phiItr->Reset();
    for (Int_t index  = 0 ; index < nAmps; ++index) {
        RooRealVar *aRho = reinterpret_cast<RooRealVar*>(_rhoItr->Next());
        RooRealVar *aPhi = reinterpret_cast<RooRealVar*>(_phiItr->Next());
        Double_t rho = aRho->getVal();
        Double_t phi = aPhi->getVal();
        _cofV[index] = TComplex(rho*TMath::Cos(phi),
                rho*TMath::Sin(phi));
    }
    return;
}

TComplex RooPWAPdf::getFFVal(const Int_t & index,
        const  Double_t p4_1[4],
        const Double_t p4_2[4],
        const Double_t p4_3[4]) const {
    TComplex partAmp;
    DecayType::DecayType modeNum = _DecayNumList[index];
    vector<Double_t> paras;
    for (int i = 0; i < _paraSize[index]; i++) {
        paras.push_back(_parameters[index][i]);
    }
    if (modeNum == 1) {
        partAmp =  PartialAmplitude(index,  paras, p4_2,  p4_3,  p4_1);
        switch (_SoloNotIdenPar) {
            case 1:
                partAmp += PartialAmplitude(index, paras, p4_3, p4_2, p4_1);
                break;
            case 2:
                partAmp += PartialAmplitude(index, paras, p4_2, p4_1, p4_3);
                break;
            case 3:
                partAmp += PartialAmplitude(index, paras, p4_1, p4_3, p4_2);
                break;
            default:
                break;
        }
    } else if (modeNum == 2) {
        partAmp = PartialAmplitude(index,  paras, p4_3,  p4_1,  p4_2);
        switch (_SoloNotIdenPar) {
            case 1:
                partAmp += PartialAmplitude(index, paras, p4_2, p4_1, p4_3);
                break;
            case 2:
                partAmp += PartialAmplitude(index, paras, p4_3, p4_1, p4_2);
                break;
            case 3:
                partAmp += PartialAmplitude(index, paras, p4_3, p4_2, p4_1);
                break;
            default:
                break;
        }
    } else if (modeNum == 3) {
        partAmp = PartialAmplitude(index,  paras, p4_1,  p4_2,  p4_3);
        switch (_SoloNotIdenPar) {
            case 1:
                partAmp += PartialAmplitude(index, paras, p4_1, p4_3, p4_2);
                break;
            case 2:
                partAmp += PartialAmplitude(index, paras, p4_3, p4_2, p4_1);
                break;
            case 3:
                partAmp += PartialAmplitude(index, paras, p4_2, p4_1, p4_3);
                break;
            default:
                break;
        }
    }
    // cout<<__func__<<endl;
    // cout<<"partal amplitude: "<<partAmp<<endl;
    return partAmp;
}

void RooPWAPdf::DIYMC(const Int_t& events, const TString& fout,
        const Int_t &sed) {
    ifstream f;
    f.open(_fracDat);
    TLorentzVector p4_1, p4_2, p4_3;
    while (1) {
        if (f.eof () != 0) break;
        Double_t E, Px, Py, Pz;
        f >> E >> Px >> Py >> Pz;
        p4_1.SetE(E); p4_1.SetPx(Px); p4_1.SetPy(Py); p4_1.SetPz(Pz);
        f >> E >> Px >> Py >> Pz;
        p4_2.SetE(E); p4_2.SetPx(Px); p4_2.SetPy(Py); p4_2.SetPz(Pz);
        f >> E >> Px >> Py >> Pz;
        p4_3.SetE(E); p4_3.SetPx(Px); p4_3.SetPy(Py); p4_3.SetPz(Pz);
        break;
    }
    f.close();
    Double_t M = (p4_1 + p4_2 + p4_3).M();
    FillParameter();
    // cout<<__func__<<endl;
    // cout<<"mass of mother: "<<M<<endl;
    TLorentzVector p4D(0, 0, 0, M);
    Double_t m1 =  p4_1.M();
    Double_t m2 =  p4_2.M();
    Double_t m3 =  p4_3.M();
    // cout<<"mass of child 1: "<<m1<<endl;
    // cout<<"mass of child 2: "<<m2<<endl;
    // cout<<"mass of child 3: "<<m3<<endl;
    Double_t masses[3] = {m1, m2, m3};
    TGenPhaseSpace event;
    event.SetDecay(p4D, 3, masses);
    Double_t weight = 0;
    int n = 0;
    TRandom3 random;
    random.SetSeed(sed);
    Double_t maxAmp = MaxAmp();
    // cout<<"Max Amplitude: "<<maxAmp<<endl;
    ofstream outdat(fout);
    outdat.precision(8);
    outdat.setf(ios::fixed);
    while (1) {
        weight =  event.Generate();
        if (random.Uniform (0, 1) > weight) continue;
        TLorentzVector *p1 = event.GetDecay(0);
        TLorentzVector *p2 = event.GetDecay(1);
        TLorentzVector *p3 = event.GetDecay(2);
        Double_t p4_1[4]={p1->E(), p1->Px(), p1->Py(), p1->Pz() };
        Double_t p4_2[4]={p2->E(), p2->Px(), p2->Py(), p2->Pz() };
        Double_t p4_3[4]={p3->E(), p3->Px(), p3->Py(), p3->Pz() };
        Double_t amp = calEva(p4_1, p4_2, p4_3);
        if (random.Uniform (0, maxAmp*1.1) >amp) continue;
        outdat << p1->E() <<"\t"
             << p1->Px() << "\t"
             << p1->Py() << "\t"
             << p1->Pz() << endl;
        outdat << p2->E() <<"\t"
             << p2->Px() << "\t"
             << p2->Py() << "\t"
             << p2->Pz() << endl;
        outdat << p3->E() <<"\t"
             << p3->Px() << "\t"
             << p3->Py() << "\t"
             << p3->Pz() << "\t1.00" << endl;
        n++;
        if (n >= events) break;
    }
    outdat.close();
}
Double_t RooPWAPdf::MaxAmp() {
    ifstream f;
    f.open(_fracDat);
    TLorentzVector p4_1, p4_2, p4_3;
    while (1) {
        if (f.eof () != 0) break;
        f >> p4_1[3] >> p4_1[0] >> p4_1[1] >> p4_1[2];
        f >> p4_2[3] >> p4_2[0] >> p4_2[1] >> p4_2[2];
        f >> p4_3[3] >> p4_3[0] >> p4_3[1] >> p4_3[2];
        break;
    }
    f.close();
    FillParameter();
    Double_t M  = (p4_1.E() + p4_2.E() + p4_3.E());
    TLorentzVector p4D(0, 0, 0, M);
    Double_t m1 =  p4_1.M();
    Double_t m2 =  p4_2.M();
    Double_t m3 =  p4_3.M();
    Double_t masses[3] = {m1, m2, m3};
    TGenPhaseSpace event;
    event.SetDecay(p4D, 3, masses);
    Double_t weight;
    Int_t events = 10000;
    int n = 0;
    TRandom *random;
    Double_t maxAmp(1);
    for (int i  = 0; i< _NameList.size(); i++) {
     // cout<<_NameList[i]<<" : "<<_cofV[i]<<endl;
    }
    while (1) {
        event.Generate();
        TLorentzVector *p1 = event.GetDecay(0);
        TLorentzVector *p2 = event.GetDecay(1);
        TLorentzVector *p3 = event.GetDecay(2);
        Double_t p4_1[4]={p1->E(), p1->Px(), p1->Py(), p1->Pz() };
        Double_t p4_2[4]={p2->E(), p2->Px(), p2->Py(), p2->Pz() };
        Double_t p4_3[4]={p3->E(), p3->Px(), p3->Py(), p3->Pz() };
        Double_t amp = calEva(p4_1, p4_2, p4_3);
        if (amp > maxAmp) {
            maxAmp = amp;
        }
        if (++n > events) break;
    }
    return maxAmp;
}
