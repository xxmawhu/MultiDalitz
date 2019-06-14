// Copyright (c) 2019-3-3 maxx
#ifndef PWAPlot_H
#define PWAPlot_H
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include <iostream>
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TColor.h"
#include <vector>

using std::cout;
using std::endl;
class PWAPlot {
 public:
    PWAPlot();
    PWAPlot(TString projectfile, TString Datadat);
    void Plotm12();
    void Plotm13();
    void Plotm23();
    void Plotdalitz();
    void PlotData(TString);
    void PlotFitResult(TString);
    void SetBinWidth(Double_t width);
    void SetTitle(TString, TString);
    void SetParName(const TString & name1,
            const TString &name2,
            const TString &name3);
    void SetLength(const TString & legStr);
    void SetFFName(const TString &FFname);
    void SetResName(const int& index, const TString &name);
    void Draw(TString name = "c");
    void test();

 private:
    void SetP11(TH1 *h);
    void SetP12(TH1 *h);
    void SetP22(TH1 *h);
    void FormatData();
    void ReadData();
    void SetBound();
    void ReadMC();
    Double_t GetMass(Int_t i);
    void PlotLengend(const std::vector<TH1D>& hists);
    void PlotEfficiency(TString name = "c");
    void PlotWavesM12(TString title = "c");
    void PlotWavesM13(TString title = "c");
    void PlotWavesM23(TString title = "c");
    TString str(Int_t kk);
    Int_t GetWaves();
    TString _projectfile, _Datadat;
    TH1D _m12_data, _m23_data, _m13_data;
    TH1D _m12_bkg, _m23_bkg, _m13_bkg;
    TH1D _m12_MC, _m23_MC, _m13_MC;
    TH2D _dalitz_MC, _dalitz_data, _dalitz_bkg;
    Double_t _binwidth;
    Double_t _m12_bounds[2], _m23_bounds[2], _m13_bounds[2];
    TString _title[5];
    TString _name[3];
    Int_t _bins;
    TString m_lengStr, m_FFName, m_ResNames[50];
    bool _isFormated, _isSetBound;
    Double_t _scaleFactor;
    void SetColor();
};
#endif
