// Copyright (c) 2019-Jun-13 maxx
#ifndef MultiDalitPdf_HH
#define MultiDalitPdf_HH

#include <math.h>
#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include "TComplex.h"
#include "TString.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TStopwatch.h"
#include "TH1F.h"
#if defined(USEROOT) || defined(__CINT__)
#include "RooStringVar.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooFitResult.h"
#include "RooArgSet.h"
#include "RooRealConstant.h"
// #include "TSpline.h"
#include "TTree.h"
// #include "TMap.h"
#else
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooAbsData.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooRandom.h"
#include "RooArgSet.h"
#include "RooRealConstant.h"
#endif
#include "RooPropogator.hh"
#include "LineShape.hh"
#include "DecayType.hh"
class MultiDalitPdf : public RooAbsPdf {
 public:
        MultiDalitPdf(const char *name, const char *title,
                RooAbsReal& _p11,
                RooAbsReal& _p12,
                RooAbsReal& _p13,
                RooAbsReal& _p14,
                RooAbsReal& _p21,
                RooAbsReal& _p22,
                RooAbsReal& _p23,
                RooAbsReal& _p24,
                RooAbsReal& _p31,
                RooAbsReal& _p32,
                RooAbsReal& _p33,
                RooAbsReal& _p34,
                const TString &PHSPDat);
        MultiDalitPdf(const MultiDalitPdf& other, const char* name = 0);
        virtual TObject* clone(const char* newname) const {
            return new MultiDalitPdf(*this, newname);
        }
        inline virtual ~MultiDalitPdf();

        void project(const char* fname);
        void setFracDat(const TString &dat);
        void setPHSPDat(const TString &dat);
        void setIdenticalParticle(const Int_t & ii = 0);
        void setrD(const Double_t &);
        void init();
        void setLambda(RooRealVar &aLambda, const Double_t &FF = 0.266);
        void FreeLineShape();

        // parameters is the list of paramter, the order: spin,
        // effective radius, mass, width( or g1, g2 ...)
        // for RBW, the last parameter must be width
        // for flatte, the last two parameter is g1, g2

        bool configMother(const TString &name, const TString &title,
                const LineShape::Shape &shape, RooArgList &parameters);
        bool addResonance(const TString &name, const TString &mothername,
                const LineShape::Shape &shape, 
                const DecayType::DecayType &modetype,
                RooRealVar &rho, RooRealVar &phi,
                RooArgList &params,
                const Int_t & angL);
        void fitFractions(const RooArgList& newPar, Bool_t
                prInt_t = kFALSE, ostream& os = std::cout);
        Int_t setPar(const RooArgList& newPar);
        Double_t MCIntG();
        void DIYMC(const Int_t& events, const TString& fout,
                const Int_t&sed);

        void test();
        Double_t calEva(const Double_t *p1, const Double_t *p2,
                const Double_t *p3) const;
        Double_t getFFVal(const Int_t &ii, const Int_t &jj,
                const Double_t _p4_1[4], const Double_t _p4_2[4],
                const Double_t _p4_3[4]);
        TComplex getFFVal(const Int_t & index,
                const  Double_t p4_1[4],
                const Double_t p4_2[4],
                const Double_t p4_3[4]) const;

 protected:
    RooRealProxy p11;
    RooRealProxy p12;
    RooRealProxy p13;
    RooRealProxy p14;
    RooRealProxy p21;
    RooRealProxy p22;
    RooRealProxy p23;
    RooRealProxy p24;
    RooRealProxy p31;
    RooRealProxy p32;
    RooRealProxy p33;
    RooRealProxy p34;
    Double_t evaluate() const;
    Double_t evaluate(Double_t _p11, Double_t _p21, Double_t _p31,
            Double_t _p41, Double_t _p12, Double_t _p22, Double_t
            _p32, Double_t _p42, Double_t _p13, Double_t _p23,
            Double_t _p33, Double_t _p43) const;

    Int_t getAnalyticalIntegral(RooArgSet& allVars,
            RooArgSet& analVars, const char* rangeName = 0) const;
    Double_t analyticalIntegral(Int_t code, const char* rangeName) const;

 private:
    Int_t *_resoParsBegin;
    // Int_t *_resoParsEnd;
    Int_t *_mothParsBegin;
    // Int_t *_mothParsEnd;
    RooListProxy _ParameterCol;
    RooListProxy _motherParameters;
    RooListProxy _OthersCol;
    RooListProxy _rhoList, _phiList;
    TIterator *_parsItr;
    TIterator *_mothItr;
    TIterator *_others;
    TIterator *_rhoItr;
    TIterator *_phiItr;
    void initialize();
    void FillParameter() const;
    Double_t ConstrFactor(const Double_t &FFa0 = 0.266) const;
    TComplex PartialAmplitude(const Int_t &index,
            const vector<Double_t> &motherparas, 
            const vector<Double_t> &paras,
            const Double_t p4_1[4],
            const Double_t p4_2[4],
            const Double_t p4_solo[4]) const;
    TComplex PartAmp(const Int_t &index,
            const vector<Double_t> &parameters,
            const Double_t p4_1[4],
            const Double_t p4_2[4],
            const Double_t p4_3[4]) const;
    void PartAmpInt();
    Double_t MaxAmp();
    vector<RooArgList> _ParaList;
    vector<TString> _NameList;
    vector<Int_t> _sourceList;
    vector<Int_t> _angL;
    vector<LineShape::Shape> _LineShapeList;
    vector<DecayType::DecayType> _DecayNumList;
    vector<TString> _motherName;
    vector<LineShape::Shape> _motherShapeList;

    // the data are used to calculate the Branchting fraction of
    // each partial waves
    TString  _fracDat , _PHSPDat;
    Int_t _SoloNotIdenPar;
    Double_t _rD;
    // vector< vector<Double_t> > _resparas;
    // max: 500
    Double_t *_resparaList;
    // max: 50
    Double_t *_motheparaList; 
    Int_t *_paraSize;
    // vector<Double_t> _rhoV;
    // vector<Double_t> _phiV;
    // vector<TComplex> _cofV;
    TComplex *_cofV;
    TComplex **_ParAmpInt;
    // vector<LineShape::Shape> _NeedSpin;
    // vector<LineShape::Shape> _NeedrRes;

    Double_t **_mcp1;
    Double_t **_mcp2;
    Double_t **_mcp3;
    Double_t *_weight;
    Int_t Nmc;
    // bool m_fixbr;
    // bool _freeShape;
    // Double_t m_fixFFa0;
    // TComplex **_Couple;

    ClassDef(MultiDalitPdf, 1)
};
#endif

