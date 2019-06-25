// Copyright (c) 2019-Jun-13 maxx
#include <fstream>
#include "MultiDalitPdf.hh"
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

ClassImp(MultiDalitPdf)
    MultiDalitPdf::MultiDalitPdf(const char *name, const char *title,
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
        _motherParameters("motherParameters", "", this),
        _OthersCol("OthersCol", "", this),
        _rhoList("rhoList", "", this),
        _phiList("phiList", "", this) {
    _parsItr = _ParameterCol.createIterator();
    _mothItr = _motherParameters.createIterator();
    _others = _OthersCol.createIterator();
    _rhoItr = _rhoList.createIterator();
    _phiItr = _phiList.createIterator();
    // _rhoList = RooArgList("rhoList");
    // _phiList = RooArgList("phiList");
    _PHSPDat = PHSPdat;
    initialize();
}


MultiDalitPdf::MultiDalitPdf(const MultiDalitPdf& other, const char* name):
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
    _motherParameters("motherParameters", this, other._motherParameters),
    _OthersCol("OthersCol", this, other._OthersCol),
    _rhoList("rhoList", this, other._rhoList),
    _phiList("phiList", this, other._phiList) {
    _parsItr = _ParameterCol.createIterator();
    _mothItr = _motherParameters.createIterator();
    _others  = _OthersCol.createIterator();
    _phiItr  = _phiList.createIterator();
    _rhoItr  = _rhoList.createIterator();
    _fracDat = other._fracDat;
    _PHSPDat = other._PHSPDat;
    _NameList = other._NameList;
    _LineShapeList = other._LineShapeList;
    _DecayNumList = other._DecayNumList;
    // _resoParsEnd = other._resoParsEnd;
    _sourceList = other._sourceList;
    _angL = other._angL;
    _motherName = other._motherName;
    _motherShapeList = other._motherShapeList;
    // _NeedSpin = other._NeedSpin;
    // _NeedrRes = other._NeedrRes;
    /*
       _mcp1 = other._mcp1;
       _mcp2 = other._mcp2;
       _mcp3 = other._mcp3;
       _weight = other._weight;
       */


    // m_fixbr = other.m_fixbr;
    // _freeShape = other._freeShape;
    //m_fixFFa0 = other.m_fixFFa0;
    Nmc = other.Nmc;
    _mcp1 = new Double_t*[Nmc];
    _mcp2 = new Double_t*[Nmc];
    _mcp3 = new Double_t*[Nmc];
    _weight = new Double_t[Nmc];
    // _Couple = new TComplex*[50];
    _ParAmpInt = new TComplex*[50];
    for (int i = 0; i < 50; i++) {
        // _Couple[i] = new TComplex[50];
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
            //_Couple[i][j] = other._Couple[i][j];
            _ParAmpInt[i][j] = other._ParAmpInt[i][j];
        }
    }
    _resoParsBegin = new Int_t[50];
    _resoParsBegin[0] = 0;
    // _resoParsEnd = new Int_t[50];
    _mothParsBegin = new Int_t[50];
    _mothParsBegin[0] = 0;
    // _mothParsEnd = new Int_t[50];
    for (Int_t i =0; i < 50; ++i) {
        _resoParsBegin[i] = other._resoParsBegin[i];
        //  _resoParsEnd[i] = other._resoParsEnd[i];
        _mothParsBegin[i] = other._mothParsBegin[i];
        // _mothParsEnd[i] = other._mothParsEnd[i];
    }
}
void MultiDalitPdf::initialize() {
    Nmc = 1E6;
    _mcp1 = new Double_t*[Nmc];
    _mcp2 = new Double_t*[Nmc];
    _mcp3 = new Double_t*[Nmc];
    _weight = new Double_t[Nmc];
    // _Couple = new TComplex*[50];
    _resoParsBegin = new Int_t[50];
    // _resoParsEnd = new Int_t[50];
    _mothParsBegin = new Int_t[50];
    // _mothParsEnd = new Int_t[50];

    _ParAmpInt = new TComplex*[50];
    for (int i = 0; i < 50; i++) {
        // _Couple[i] = new TComplex[50];
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
        cout << "Error::can't open input file " << _PHSPDat << endl;
        exit(0);
    } else {
        cout << "Info:: open " << _PHSPDat << endl;
        cout << "the largest number of events is " << Nmc << endl;
        cout << "at lest 100 times of data is recommend" << endl;
        cout << "1st line four-momentum of K+" << endl;
        cout << "2nd line four-momentum of K-" << endl;
        cout << "3rd line four-momentum of pi0, weight" << endl;
        cout << "the order of four-momentum should be " << endl
            << "E, px, py, pz" << endl ;
    }

    Int_t i = 0;
    while (fscanf (fp, "%lf%lf%lf%lf\n%lf%lf%lf%lf\n%lf%lf%lf%lf%lf\n",
                &ft1 , &fx1 , &fy1 , &fz1 ,
                &ft2 , &fx2 , &fy2 , &fz2 ,
                &ft3 , &fx3 , &fy3 , &fz3 , &weight) != EOF) {
        if (i >= Nmc) break;
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

    cout << "Read from " << _PHSPDat << "successful!, the total number of PHSP events is "
        << Nmc << endl << endl << endl;
}


Double_t MultiDalitPdf::evaluate() const {
    Double_t p4_1[4] = {p11, p12, p13, p14};
    Double_t p4_2[4] = {p21, p22, p23, p24};
    Double_t p4_3[4] = {p31, p32, p33, p34};
    return calTotalWidth(p4_1, p4_2, p4_3);
}

Double_t MultiDalitPdf::evaluate(Double_t _p11, Double_t _p12, Double_t _p13, 
        Double_t _p14, Double_t _p21, Double_t _p22, Double_t _p23, 
        Double_t _p24, Double_t _p31, Double_t _p32, Double_t _p33,
        Double_t _p34) const {
    Double_t p4_1[4] = {_p11, _p12, _p13, _p14};
    Double_t p4_2[4] = {_p21, _p22, _p23, _p24};
    Double_t p4_3[4] = {_p31, _p32, _p33, _p34};
    return calTotalWidth(p4_1, p4_2, p4_3);
}

Int_t MultiDalitPdf::getAnalyticalIntegral(RooArgSet& allVars,
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

Double_t MultiDalitPdf::analyticalIntegral(Int_t code, const char* rangeName)
    const {
    assert(code == 1);

    Double_t sum = 0;
    for (int i = 0; i < Nmc; i++) {
        double eva = calTotalWidth(_mcp1[i], _mcp2[i], _mcp3[i]);
        sum +=  eva*_weight[i] / Nmc;
    }
    return sum;
}


void MultiDalitPdf::setFracDat(const TString &dat) {
    _fracDat = dat;
}


void MultiDalitPdf::project(const char* proFileName) {
    for (int i = 0; i < _rhoList.getSize(); i++) {
        cout << i << "-th waves" << endl;
        cout <<"\trho = "
            << reinterpret_cast<RooRealVar*>(_rhoList.at(i))->getVal() << endl;
        cout <<"\tphi = "
            << reinterpret_cast<RooRealVar*>(_phiList.at(i))->getVal() << endl;
    }
    TFile f(proFileName, "recreate");
    TTree t("project", "the project of the mulitdalitz plot fit result");

    Int_t nAmps = _NameList.size();
    const Int_t length = nAmps*nAmps;
    Double_t *FF = new Double_t[length];

    //eva total amplitude^2
    Double_t eva;
    t.Branch("eva", &eva, "eva/D");

    TH1D hkk("hkk", "mass of K+ K-", 40, 0.96, 1.34);
    TH1D hkstp("hkstp", "mass of K+ pi0", 40, 0.62, 1.05);
    TH1D hkstm("hkstm", "mass of K- pi0", 40, 0.62, 1.05);
    TH1D heta("heta", "mass of K+ K- pi0", 40, 1.3, 1.6);
    
    //  m12: mkk, m13: mkstp, m23: mkstm 
    Double_t m12, m13, m23;
    t.Branch("mkk", &m12, "mkk/D");
    t.Branch("mkstp", &m13, "mkstp/D");
    t.Branch("mkstm", &m23, "mkstm/D");

    // four momentum
    Double_t m_p4_1[4], m_p4_2[4], m_p4_3[4];
    t.Branch("p4kp", m_p4_1, "p4kp[4]/D");
    t.Branch("p4km", m_p4_2, "p4km[4]/D");
    t.Branch("p4pi0", m_p4_3, "p4pi0[4]/D");

    // each partial wave contribution
    Int_t waves = length;
    t.Branch("length", &waves, "length/I");
    t.Branch("FF", FF, "FF[length]/D");
    
    TComplex *ffval = new TComplex[nAmps];
    double max=0.0;
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

        eva =  _weight[i] * calTotalWidth(_mcp1[i],  _mcp2[i],  _mcp3[i]);
        if (eva > max) {
            max = eva;
            cout << "project max = " << max << endl;
        }
        waves = nAmps * nAmps;
        for (int ii = 0; ii < nAmps; ii++) {
            ffval[ii] = getPartAmp(ii, _mcp1[i], _mcp2[i], _mcp3[i]);
        }

        for (Int_t ii = 0; ii < nAmps; ++ii) {
            for (Int_t jj = 0; jj < nAmps; ++jj) {
                Int_t index = ii*nAmps + jj;
                FF[index] = 0.5*( ffval[ii] * TComplex::Conjugate(ffval[jj]) 
                       + ffval[jj] * TComplex::Conjugate(ffval[ii]) );
            }
        }
        hkk.Fill(m12, eva);
        hkstp.Fill(m13, eva);
        hkstm.Fill(m23, eva);
        heta.Fill((p4_1 + p4_2 + p4_3).M(), eva);
        t.Fill();
    }
    t.Write();
    hkk.Write();
    hkstp.Write();
    hkstm.Write();
    heta.Write();
    delete[] FF;
    // delete cof;
    delete []ffval;
    f.Close();
}


void MultiDalitPdf::fitFractions(ostream&os=cout) {
    //FillParameter();
    double p4_1[4], p4_2[4], p4_3[4], weight;
    int nWaves = _rhoList.getSize();
    double **FF = new double*[nWaves];
    for (int i = 0; i < nWaves; i++) {
        FF[i] = new double[nWaves];
    }
    for (int i = 0; i < nWaves; i++) {
        for (int j = 0; j < nWaves; j++) {
            FF[i][j] = 0;
        }
    }
    TComplex *tmp = new TComplex[nWaves];
    int nEvents = 0.0;
    
    for (int i = 0; i < _rhoList.getSize(); i++) {
        cout << i << "-th waves" << endl;
        cout <<"\trho = "
            << reinterpret_cast<RooRealVar*>(_rhoList.at(i))->getVal() << endl;
        cout <<"\tphi = "
            << reinterpret_cast<RooRealVar*>(_phiList.at(i))->getVal() << endl;
    }

    FILE *fp;
    fp=fopen(_fracDat, "r");
    double max = 0.0;
    while (fscanf (fp, "%lf%lf%lf%lf\n%lf%lf%lf%lf\n%lf%lf%lf%lf%lf\n",
                &p4_1[0], &p4_1[1], &p4_1[2], &p4_1[3],
                &p4_2[0], &p4_2[1], &p4_2[2], &p4_2[3],
                &p4_3[0], &p4_3[1], &p4_3[2], &p4_3[3],
                &weight) != EOF) {
        TComplex amp(0);
        for (int i = 0; i < nWaves; i++) {
            tmp[i] = this->getPartAmp(i, p4_1, p4_2, p4_3);
            amp += tmp[i];
        }
        if (amp.Rho2() > max) {
            max = amp.Rho2();
            cout << "fitFractions:: amp = " << amp.Rho2() << endl;
        }
        for (int i = 0; i < nWaves; i++) {
            for (int j = 0; j < nWaves; j++) {
                FF[i][j] += 0.5*( tmp[i] *TComplex::Conjugate(tmp[j])
                       + tmp[j] * TComplex::Conjugate(tmp[i])).Re();
            }
        }
        nEvents += 1;
    }
    cout << "Cal Fit Fraction, total events " << nEvents << endl;
    double sum = 0.0;

    for (int i = 0; i < nWaves; i++) {
        for (int j = 0; j < nWaves; j++) {
            FF[i][j] = FF[i][j]/nEvents;
            sum += FF[i][j];
        }
    }
    cout << "the partial fraction is" << endl << endl;
    for (int i = 0; i < nWaves; i++) {
        cout << "\t" <<std::setw(10) <<  _NameList[i] << " fraction is " 
            << std::setw(4)<< std::setprecision(2) << FF[i][i]/sum*100 << "%"
            << endl;
    }
    fclose(fp);

    delete []FF;
    delete []tmp;
}



Int_t MultiDalitPdf::setPar(const RooArgList& newPar) {
    return 0;
}
bool MultiDalitPdf::configMother(const TString &name, const TString &title,
        const LineShape::Shape &shape, RooArgList &parameters) {
    // name: the mother particle name, as well as title
    if (std::find (_motherName.begin (), _motherName.end (), name) != _motherName.end()) {
        std::cerr << "Error:: The mother particle : " << name
            <<" is already exist" <<std::endl;
        exit(0);
        return false;
    } else {
        _motherName.push_back(name);
    } 
    cout << "Inf::set for " << name << endl;
    cout << std::setw(20) << "name" << std::setw(20) << name << endl;
    cout << std::setw(20) << "LineShape";
    switch (shape) {
        case LineShape::Flat :
            cout << std::setw(20) << "Flat" << endl;
            break;
        case LineShape::RBW :
            cout << std::setw(20)<< "RBW" << endl;
            break;
        case LineShape::a980_0 :
            cout << std::setw(20)<< "a980_0" << endl;
            break;
        case LineShape::a980_p :
            cout << std::setw(20)<< "a980_p" << endl;
            break;
        case LineShape::a980_p_3C :
            cout << std::setw(20)<< "a980_p_3C" << endl;
            break;
        case LineShape::K1430_p :
            cout << std::setw(20)<< "K1430_p" << endl;
            break;
        case LineShape::K1430_0 :
            cout << std::setw(20)<< "K1430_0" << endl;
            break;
        case LineShape::GS :
            cout << std::setw(20)<< "GS" << endl;
            break;
        default:
            cout << std::setw(20)<< "NULL" << endl;
            break;
    }
    _motherShapeList.push_back(shape);
    _motherParameters.add(parameters);
    int index = _motherShapeList.size() -1;
    _mothParsBegin[index+1] = _mothParsBegin[index] + parameters.getSize();
    cout << endl << endl;
    return true;
    // if (index == 0) {
    //     _mothParsBegin[0] = 0;
    //     _mothParsEnd[0] = parameters.getSize();
    // } else {
    //     _mothParsBegin[index] = _mothParsBegin[index-1] + parameters.getSize();
    //     _mothParsEnd[index] = _mothParsBegin[index] + parameters.getSize();
    // }
}


bool MultiDalitPdf::addResonance(const TString &name, const TString &mothername,
        const LineShape::Shape &lineShape,
        const DecayType::DecayType &modetype,
        RooAbsReal &rho, RooAbsReal &phi,
        RooArgList &params, const Int_t &angL) {
    // mothername used to match the mother, the "mothername" must be add in the
    // configMother
    // fix the first amplitude rho and phi
    //if (_rhoList.getSize() == 0) {
    //    reinterpret_cast<RooRealVar*>(rho).setConstant();
    //    reinterpret_cast<RooRealVar*>(phi).setConstant();
    //}

    // check the resonance list
    if (std::find (_NameList.begin (), _NameList.end (), name)!= _NameList.end()) {
        std::cerr << "Error:: The resonance : " << name
            <<" is already exist" <<std::endl;
        exit(0);
        return false;
    } 
    cout << "Inf::Now add a resonance in the decay of " << mothername << endl;
    _rhoList.add(rho);
    _phiList.add(phi);
    _ParameterCol.add(params);
    _angL.push_back(angL);

    int index = _rhoList.getSize();
    _resoParsBegin[index] = _resoParsBegin[index-1] + params.getSize();
    // if (index = 0){
    //     _resoParsBegin[index] = 0;
    //     _resoParsEnd[index] = params.getSize();
    // } else {
    //     _resoParsBegin[index] = _resoParsEnd[index-1];
    //     _resoParsEnd[index] = _resoParsBegin[index] + params.getSize();
    // }

    // fille the source
    vector<TString>::iterator mitr = std::find(_motherName.begin(), 
            _motherName.end(), mothername);
    if (mitr == _motherName.end()) { 
        std::cerr << "Error::no such mother particle \"" << mothername 
            << "\", please check it!" << endl;
        exit(0);
    }
    int source = mitr - _motherName.begin(); 
    cout <<std::setw(3) << _LineShapeList.size()+1 << ") the resonce comes from " 
        << source << ", i.e " << mothername <<  endl;
    _sourceList.push_back(source);
    // _ParaList.push_back(params);
    _LineShapeList.push_back(lineShape);
    _NameList.push_back(name);
    _DecayNumList.push_back(modetype);
    RooRealVar *spin = reinterpret_cast<RooRealVar*>(params.at(0));
    // delete arho;
    // delete aphi;
    cout << std::setw(20) << "name" << std::setw(10) << name << endl;
    cout << std::setw(20) << "LineShape";
    switch (lineShape) {
        case LineShape::Flat :
            cout << std::setw(10) << "Flat" << endl;
            break;
        case LineShape::RBW :
            cout << std::setw(10)<< "RBW" << endl;
            break;
        case LineShape::a980_0 :
            cout << std::setw(10)<< "a980_0" << endl;
            break;
        case LineShape::a980_p :
            cout << std::setw(10)<< "a980_p" << endl;
            break;
        case LineShape::a980_p_3C :
            cout << std::setw(10)<< "a980_p_3C" << endl;
            break;
        case LineShape::K1430_p :
            cout << std::setw(10)<< "K1430_p" << endl;
            break;
        case LineShape::K1430_0 :
            cout << std::setw(10)<< "K1430_0" << endl;
            break;
        case LineShape::GS :
            cout << std::setw(10)<< "GS" << endl;
            break;
        default:
            cout << std::setw(10)<< "NULL" << endl;
            break;
    }
    cout << std::setw(20) << "DecayType" ;
    switch (modetype) {
        case DecayType:: AB:
           cout << std::setw(10)  << "AB"  << endl;
           break;
        case DecayType:: BC:
           cout << std::setw(10)  << "BC"  << endl;
           break;
        case DecayType:: AC:
           cout << std::setw(10) << "AC"  << endl;
           break;
        default :
           cout  << "NULL"  << endl;
    }
    cout << std::setw(20) << "Spin" << std::setw(10) << spin->getVal() << endl;
    cout << std::setw(20) << "Angular Moment" << std::setw(10) << angL << endl;
    cout << endl << endl;
    return true;
}


TComplex MultiDalitPdf::PartialAmplitude(const Int_t & index,
        const vector<Double_t> &mothpars,
        const vector<Double_t> &resparas,
        const Double_t p4_1[4], const Double_t p4_2[4],
        const Double_t p4_solo[4]) const {
    // obtain the config for the mother
    int mothIdx = _sourceList[index];
    // cout << "PartialAmplitude:: mothIdx = " << _sourceList[index] << endl;
    LineShape::Shape motherShape = _motherShapeList[mothIdx];
    // cout << "PartialAmplitude:: motherShape = " << motherShape << endl;
    

    int resParsBegin = _resoParsBegin[index];
    int resParsEnd = _resoParsBegin[index+1];
    // order spin, effective radius, mass, width,
    Double_t reta = mothpars[1];
    // vector<Double_t> mothPars;
    // for (int i = mothParsBegin; i < mothParsEnd; ++i) {
    //     mothPars.push_back(_motherParameters.at(i)->getVal());
    // }
    Double_t p4total[4] = {
        p4_1[0] + p4_2[0] + p4_solo[0],
        p4_1[1] + p4_2[1] + p4_solo[1],
        p4_1[2] + p4_2[2] + p4_solo[2],
        p4_1[3] + p4_2[3] + p4_solo[3]
    };

    Int_t AngEta = _angL[index];
    Int_t Ang =  resparas[0];
    // cout << "PartialAmplitude:: angEta " << AngEta << endl;
    // cout << "PartialAmplitude:: angL " << Ang << endl;
    Double_t rRes = 3.0;
    if (resparas.size() > 1) {
        rRes = resparas[1];
    }
    double rEta  = mothpars[1];
    // cout<<__func__<<endl;

    TComplex _pro(1, 0);
    Double_t Sa = RooP4Vector::dot(p4_1, p4_1);
    // cout<<"Sa:"<<Sa<<endl;
    Double_t Sb = RooP4Vector::dot(p4_2, p4_2);
    // cout<<"Sb: "<<Sb<<endl;
    Double_t s = 2*RooP4Vector::dot(p4_1, p4_2) + Sa + Sb;
    LineShape::Shape lineShape = _LineShapeList[index];
    _pro = Propagator::getVal(lineShape, s, Sa, Sb, resparas);

    // cout << "PartialAmplitude resonance pro = " << _pro << endl;
    
    TComplex propMoth = Propagator::getVal(motherShape, 
            RooP4Vector::dot(p4total, p4total), 
            s, // mass^2 of resonance
            RooP4Vector::dot(p4_solo, p4_solo), // mass^2 of solo particle
            mothpars);
    
    // cout << "PartialAmplitude mother pro = " << propMoth << endl;

    // cout<<"S: "<<s<<endl;

    Double_t p4R[4] = {p4_1[0]+ p4_2[0], p4_1[1]+p4_2[1],
        p4_1[2]+p4_2[2], p4_1[3]+p4_2[3] };

    // amplitude for eta(1405)->a_0 pi0, a_0->K+ K-
    // M = RBW_eta  Barrier_eta * Barrier_a0 * SpinFactor


    Double_t barrFactorEta = RooBarrier::Factor(AngEta, p4R, p4_solo, rEta);
    Double_t barrFactorRes = RooBarrier::Factor(Ang, p4_1, p4_2, rRes);
    Double_t spinFactor = RooSpinFactor::SpinFactor(AngEta, p4_1, p4_2,
            p4_solo);
    // cout<<"spinFactor: "<<spinFactor<<endl;
    // cout<<"rRes: "<<rRes<<endl;
    // cout<<"pro: "<<_pro<<endl;
    // cout<<"barrFactorRes: "<<barrFactorRes<<endl;
    // cout<<"barrFactorEta: "<<barrFactorEta<<endl;
    return propMoth * spinFactor * _pro * barrFactorRes * barrFactorEta;
}
void MultiDalitPdf::init() {
}
void MultiDalitPdf::test() {
    cout << "total partical waves " << _rhoList.getSize() << endl;
    cout << "total parameters " << _ParameterCol.getSize() << endl;
    for (int i = 0; i < _rhoList.getSize(); i++) {
        cout <<"check wave " << i+1 << endl;
        int begin = _resoParsBegin[i];
        int end = _resoParsBegin[i+1];
        cout << "\tTotal #parameters " << end - begin << endl;
        RooRealVar *par(0);
        for (int i = begin; i < end; i++) {
            par = reinterpret_cast<RooRealVar*>(_ParameterCol.at(i));
            // aPara = reinterpret_cast<RooRealVar*>(_parsItr->Next());
            cout << std::setw(15) << par->GetName() 
                << std::setw(10) << par->getVal() << endl ;
        }
    }
    cout << endl << endl;
    cout << "Inf::fill parameters" << endl << endl << endl;
    cout << "total MC sample: "<< Nmc << endl;
    cout << "mass of particle 1 is " << TMath::Sqrt(RooP4Vector::dot(_mcp1[0], _mcp1[0])) << endl;
    cout << "mass of particle 2 is " << TMath::Sqrt(RooP4Vector::dot(_mcp2[0], _mcp2[0])) << endl;
    cout << "mass of particle 3 is " << TMath::Sqrt(RooP4Vector::dot(_mcp3[0], _mcp3[0])) << endl;
    for (int i = 0; i < 10; ++i) {
        cout << "total Width " << calTotalWidth(_mcp1[i], _mcp2[i], _mcp3[i]) << endl;
        for (int j = 0; j < _rhoList.getSize(); j++) {
            TComplex parAmp = getPartAmp(j, _mcp1[i], _mcp2[i], _mcp3[i]);
            cout << "\tparAmp " << parAmp << endl;
        }
    }
    // check the rho of i-th partial waves    
    RooRealVar *rho(0);
    RooRealVar *phi(0);
    for (int i = 0; i < _rhoList.getSize(); i++) {
        cout << i << "-th waves" << endl;
        cout <<"\trho = "
            << reinterpret_cast<RooRealVar*>(_rhoList.at(i))->getVal() << endl;
        cout <<"\tphi = "
            << reinterpret_cast<RooRealVar*>(_phiList.at(i))->getVal() << endl;
    }
    return;
}
Double_t MultiDalitPdf::calTotalWidth(const Double_t p4_1[4],
        const Double_t p4_2[4], const Double_t p4_3[4]) const {
    Int_t nAmps = _NameList.size();
    TComplex _Amplitude(0);
    for (int index  = 0; index < nAmps; ++index) {
        _Amplitude += getPartAmp(index, p4_1, p4_2, p4_3);
    }
    Double_t value = _Amplitude.Rho2();
    return (value <= 0) ? 1e-20 : value;
}

Double_t MultiDalitPdf::MCIntG() {
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
        _sum += weight * calTotalWidth(p4_1, p4_2, p4_3);
        i += 1;
    }
    return _sum/i;
}


void MultiDalitPdf::setPHSPDat(const TString &dat) {
    _PHSPDat = dat;
}

Double_t MultiDalitPdf::getCrossWidth(const Int_t &ii, const Int_t &jj,
        const  Double_t p4_1[4], const Double_t p4_2[4],
        const Double_t p4_3[4]) {
    // return  (A_i * conj[A_j] + conj[A_i] * A_j)/2
    TComplex _Ampii = getPartAmp(ii, p4_1, p4_2, p4_3);
    TComplex _Ampjj = getPartAmp(jj, p4_1, p4_2, p4_3);
    TComplex module = _Ampii*TComplex::Conjugate(_Ampjj)
        + _Ampjj * TComplex::Conjugate(_Ampii);
    return 0.5*module.Re();
}

inline  MultiDalitPdf::~MultiDalitPdf() {
    for (Int_t i  = 0 ; i < Nmc;  i++) {
        delete[] _mcp1[i];
        delete[] _mcp2[i];
        delete[] _mcp3[i];
    }
    delete[] _weight;

}


/*
Double_t MultiDalitPdf::ConstrFactor(const Double_t &FFa0) const {
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
*/


/*
void MultiDalitPdf::PartAmpInt() {
    // Double_t
    int nAmps = _NameList.size();
    RooRealVar *aPara(0);
    Int_t i = 0;
    Double_t _sum = 0;
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
*/


// void MultiDalitPdf::setLambda(RooRealVar &aLambda, const Double_t &FF) {
//    m_fixbr = true;
//    m_fixFFa0 = FF;
//    _OthersCol.add(aLambda);
// }

// void MultiDalitPdf::FreeLineShape() {
//    _freeShape = true;
// }

TComplex MultiDalitPdf::getPartAmp(const Int_t & index,
        const  Double_t p4_1[4],
        const Double_t p4_2[4],
        const Double_t p4_3[4]) const {
    TComplex partAmp(0);

    int source = _sourceList[index];

    // cout << "getPartAmp:: source = " << source << endl; 
    vector<Double_t> motherparas;
    int mIndex = _sourceList[index];
    for (int i = _mothParsBegin[mIndex]; i < _mothParsBegin[mIndex+1]; i++) {
        double tmp = reinterpret_cast<RooRealVar*>(_motherParameters.at(i))->getVal();
        motherparas.push_back(tmp);
    }
    vector<Double_t> resparas;
    for (int i = _resoParsBegin[index]; i < _resoParsBegin[index+1]; i++) {
        double tmp = reinterpret_cast<RooRealVar*>(_ParameterCol.at(i))->getVal();
        resparas.push_back(tmp);
    }

    DecayType::DecayType modeNum = _DecayNumList[index];
    // cout << "getPartAmp modeNum = " << modeNum << endl;
    // cout << "getPartAmp index = " << index << endl;
    if (modeNum == 1) {
        partAmp = PartialAmplitude(index, motherparas, resparas, p4_2, p4_3, p4_1);
    } else if (modeNum == 2) {
        partAmp = PartialAmplitude(index, motherparas, resparas, p4_3, p4_1, p4_2);
    } else {
        partAmp = PartialAmplitude(index, motherparas, resparas, p4_1, p4_2, p4_3);
    }
    
    // obtain the rho and phi
    double rho = reinterpret_cast<RooRealVar*>(_rhoList.at(index))->getVal();
    double phi = reinterpret_cast<RooRealVar*>(_phiList.at(index))->getVal();
    // cout << "getPartAmp:: rho = " << rho << endl;
    // cout << "getPartAmp:: phi = " << phi << endl;
    // cout<<"getPartAmp partal amplitude: "<<partAmp<<endl;
    return partAmp * TComplex(rho*TMath::Cos(phi), rho*TMath::Sin(phi));
}

void MultiDalitPdf::DIYMC(const Int_t& events, const TString& fout,
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
    //FillParameter();
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
        Double_t p4_1[4]={p1->E(), p1->Px(), p1->Py(), p1->Pz()};
        Double_t p4_2[4]={p2->E(), p2->Px(), p2->Py(), p2->Pz()};
        Double_t p4_3[4]={p3->E(), p3->Px(), p3->Py(), p3->Pz()};
        Double_t amp = calTotalWidth(p4_1, p4_2, p4_3);
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
Double_t MultiDalitPdf::MaxAmp() {
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
    //FillParameter();
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
    while (1) {
        event.Generate();
        TLorentzVector *p1 = event.GetDecay(0);
        TLorentzVector *p2 = event.GetDecay(1);
        TLorentzVector *p3 = event.GetDecay(2);
        Double_t p4_1[4]={p1->E(), p1->Px(), p1->Py(), p1->Pz() };
        Double_t p4_2[4]={p2->E(), p2->Px(), p2->Py(), p2->Pz() };
        Double_t p4_3[4]={p3->E(), p3->Px(), p3->Py(), p3->Pz() };
        Double_t amp = calTotalWidth(p4_1, p4_2, p4_3);
        if (amp > maxAmp) {
            maxAmp = amp;
        }
        if (++n > events) break;
    }
    return maxAmp;
}
