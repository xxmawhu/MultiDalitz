// Copyright (c) 2019-3-3 Ma Xinxin
#include <fstream>
#include <map>
#include <string>
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TPad.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TColor.h"
#include <iostream>
#include <sstream>
#include "PWAPlot.hh"
#include "bes3plotstyle.hh"
#include "LauConstants.hh"
using std::ifstream;
using std::cout;
using std::endl;
using std::string;

PWAPlot::PWAPlot(TString project, TString Datadat) {
    _projectfile = project;
    _Datadat = Datadat;
    _binwidth = 20E-3;  // 30 MeV
    _bins = 0;
    _isFormated = false;
    _isSetBound = false;
    _m12_data.Reset();
    _m23_data.Reset();
    _m13_data.Reset();
    _m12_bkg.Reset();
    _m23_bkg.Reset();
    _m13_bkg.Reset();
    _m12_MC.Reset();
    _m23_MC.Reset();
    _m13_MC.Reset();
    _dalitz_data.Reset();
    _m12_bounds[0] = 0;
    _m12_bounds[1] = 0;
    _m13_bounds[0] = 0;
    _m13_bounds[1] = 0;
    _m23_bounds[0] = 0;
    _m23_bounds[1] = 0;
    _name[0] = TString("#bar{K}^{0}_{1}");
    _name[1] = TString("#bar{K}^{0}_{2}");
    _name[2] = TString("K^{+}");
    m_lengStr = "length";
    m_FFName = "FF";
    for (int i = 0; i < 50; ++i) {
        m_ResNames[i] = TString("");
    }
}
Double_t PWAPlot::GetMass(Int_t i) {
    ifstream f;
    f.open(_Datadat);
    TLorentzVector p4;
    for (Int_t j = 0; j < i; j++) {
        f >> p4[3] >> p4[0] >> p4[1] >> p4[2];
    }
    f.close();
    return p4.M();
}
void PWAPlot::SetBinWidth(Double_t width) {
    _binwidth = width;
}

void PWAPlot::SetBound() {
    if (_isSetBound) return;
    Double_t m1 = GetMass(1);
    Double_t m2 = GetMass(2);
    Double_t m3 = GetMass(3);

    Double_t totwidth = (LauConstants::mD - m1 - m2 - m3)*1.3;

    _bins = totwidth/_binwidth +1;
    Double_t dM = (_binwidth * _bins - LauConstants::mD + m1 + m2 + m3)/2;

    _m12_bounds[0] = m1 + m2 - dM;
    _m12_bounds[1] = LauConstants::mD - m3 + dM;


    _m13_bounds[0] = m1 + m3 - dM;
    _m13_bounds[1] = LauConstants::mD - m2 + dM;

    _m23_bounds[0] = m3 + m2 - dM;
    _m23_bounds[1] = LauConstants::mD - m1 + dM;
    _isSetBound = true;
}
void PWAPlot::test() {
    cout << "mass 1 "  << GetMass(1) << endl;
    cout << "mass 2 "  << GetMass(2) << endl;
    cout << "mass 3 "  << GetMass(3) << endl;
    ReadData();
    _m12_data.Print("v");
    cout << "mean:" << _m12_data.GetMean() << endl;
    ReadMC();
    bes3plotstyle::SetStyle();
    TCanvas c("c", "c",  800,  600);
    c.Divide(2, 2);
    c.cd(1);
    Plotm12();
    c.Update();
    c.cd(2);
    Plotm13();
    c.Update();
    c.cd(3);
    Plotm23();
    c.Update();
    c.cd(4);
    Plotdalitz();
    c.Update();
    c.SaveAs("c.eps");
}

void PWAPlot::ReadData() {
    ifstream f;
    f.open(_Datadat);
    TLorentzVector p4_1, p4_2, p4_3;
    if (_m12_bounds[1] < 0.001) SetBound ();
    _m12_data = TH1D("_m12_data", "", _bins,  _m12_bounds[0], _m12_bounds[1]);
    _m13_data = TH1D("_m13_data", "", _bins,  _m13_bounds[0], _m13_bounds[1]);
    _m23_data = TH1D("_m23_data", "", _bins,  _m23_bounds[0], _m23_bounds[1]);
    _m12_bkg = TH1D("_m12_bkg", "", _bins,  _m12_bounds[0], _m12_bounds[1]);
    _m13_bkg = TH1D("_m13_bkg", "", _bins,  _m13_bounds[0], _m13_bounds[1]);
    _m23_bkg = TH1D("_m23_bkg", "", _bins,  _m23_bounds[0], _m23_bounds[1]);
    _dalitz_data = TH2D("_dalitz_data", "",
            _bins, pow(_m13_bounds[0], 2), pow(_m13_bounds[1], 2),
            _bins, pow(_m23_bounds[0], 2), pow(_m23_bounds[1], 2) );
    _dalitz_bkg = TH2D("_dalitz_bkg", "",
            _bins, pow(_m13_bounds[0], 2), pow(_m13_bounds[1], 2),
            _bins, pow(_m23_bounds[0], 2), pow(_m23_bounds[1], 2) );
    Double_t weight;
    while (!f.eof()) {
        f >> p4_1[3] >> p4_1[0] >> p4_1[1] >> p4_1[2];
        f >> p4_2[3] >> p4_2[0] >> p4_2[1] >> p4_2[2];
        f >> p4_3[3] >> p4_3[0] >> p4_3[1] >> p4_3[2] >> weight;
        if (weight >0) {
            _m12_data.Fill((p4_1+p4_2).M(), weight);
            _m23_data.Fill((p4_2+p4_3).M(), weight);
            _m13_data.Fill((p4_1+p4_3).M(), weight);
            _dalitz_data.Fill(pow((p4_1+p4_3).M(), 2)
                    , pow((p4_2+p4_3).M(), 2), weight);
        } else if (weight < 0) {
            _m12_bkg.Fill((p4_1+p4_2).M(), -weight);
            _m13_bkg.Fill((p4_1+p4_3).M(), -weight);
            _m23_bkg.Fill((p4_3+p4_2).M(), -weight);
            _dalitz_bkg.Fill(pow((p4_1+p4_3).M(), 2),
                    pow((p4_2+p4_3).M(), 2),  -weight);
        }
    }
    _isFormated = false;
}
void PWAPlot::SetLength(const TString & legStr) {
  m_lengStr = legStr;
}
void PWAPlot::SetFFName(const TString &FFname) {
  m_FFName = FFname;
}
Int_t PWAPlot::GetWaves() {
    TFile *f = new TFile(_projectfile, "read");
    TTree *tree = reinterpret_cast<TTree*>(f->Get("project"));
    Int_t length;
    tree->SetBranchAddress(m_lengStr, &length);
    tree->GetEntry(0);
    Int_t waves = sqrt(length);
    delete tree;

    f->Close();
    delete f;
    return waves;
}

void PWAPlot::ReadMC() {
    Int_t waves = GetWaves();

    TFile *f = new TFile(_projectfile, "read");
    TTree *tree = reinterpret_cast<TTree*>(f->Get("project"));
    Double_t m12, m13, m23, eva;
    tree->SetBranchAddress("m12", &m12);
    tree->SetBranchAddress("m23", &m23);
    tree->SetBranchAddress("m13", &m13);
    tree->SetBranchAddress("eva", &eva);

    if (!_isSetBound) SetBound ();
    _m12_MC = TH1D("_m12_MC", "", _bins,  _m12_bounds[0],
            _m12_bounds[1]);
    _m13_MC = TH1D("_m13_MC", "", _bins, _m13_bounds[0],
            _m13_bounds[1]);
    _m23_MC = TH1D("_m23_MC", "", _bins,  _m23_bounds[0],
            _m23_bounds[1]);
    _dalitz_MC = TH2D("_dalitz_MC", "",
            _bins, pow(_m13_bounds[0], 2), pow(_m13_bounds[1], 2),
            _bins, pow(_m23_bounds[0], 2), pow(_m23_bounds[1], 2) );
    for (Int_t i = 0 ; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        _m12_MC.Fill(m12, eva);
        _m13_MC.Fill(m13, eva);
        _m23_MC.Fill(m23, eva);
        _dalitz_MC.Fill(pow(m13, 2), pow(m23, 2), eva);
    }
    return;
    delete tree;
    f->Close();
    delete f;
    _isFormated = false;
    return;
}
void PWAPlot::FormatData() {
    if (_isFormated) return;
    _m12_data.SetMarkerSize(0.8);
    _m12_data.SetMarkerColor(1);
    _m12_MC.SetLineColor(kBlue);
    _m12_MC.SetLineWidth(2);
    _m12_MC.SetMarkerSize(0.001);
    _m12_bkg.SetLineColor(kRed);
    _m12_bkg.SetLineWidth(1);

    _m23_data.SetMarkerSize(0.8);
    _m23_data.SetMarkerColor(1);
    _m23_MC.SetLineColor(kBlue);
    _m23_MC.SetLineWidth(2);
    _m23_MC.SetMarkerSize(0.001);
    _m23_bkg.SetLineColor(kRed);
    _m23_bkg.SetLineWidth(1);

    _m13_data.SetMarkerSize(0.8);
    _m13_data.SetMarkerColor(1);
    _m13_MC.SetLineColor(kBlue);
    _m13_MC.SetLineWidth(2);
    _m13_MC.SetMarkerSize(0.001);
    _m13_bkg.SetLineColor(kRed);
    _m13_bkg.SetLineWidth(1);

    _dalitz_bkg.SetMarkerSize(0.001);
    _dalitz_MC.SetMarkerSize(0.001);

    _dalitz_MC.SetLineWidth(kBlue);
    _dalitz_MC.SetLineWidth(0.5);
    _dalitz_MC.SetLineStyle(3);

    _dalitz_data.SetMarkerSize(0.001);
    _dalitz_data.SetLineWidth(kBlack);
    _dalitz_data.SetLineWidth(0.5);
    _dalitz_data.SetLineStyle(2);

    bes3plotstyle::FormatYAxis(_m13_data.GetYaxis());
    bes3plotstyle::FormatYAxis(_m23_data.GetYaxis());

    bes3plotstyle::FormatAxis(_m12_data.GetXaxis());
    bes3plotstyle::FormatAxis(_m13_data.GetXaxis());
    bes3plotstyle::FormatAxis(_m23_data.GetXaxis());
    _isFormated = true;
}
void PWAPlot::Plotm12() {
    FormatData();
    TPad *pad_1 = new TPad("pad_1", "", 0.0, 0.00,  1.0 ,  1.0);
    bes3plotstyle::Format(pad_1);
    pad_1->Draw();
    pad_1->cd();
    Double_t dataNum = _m12_data.Integral();
    Double_t bkgNum = _m12_bkg.Integral();
    Double_t mcNum = _m12_MC.Integral();
    _scaleFactor = (dataNum - bkgNum)/mcNum;
    _m12_MC.Scale(_scaleFactor);
    _m12_MC.Add(&_m12_bkg);

    Double_t max = _m12_data.GetMaximum();
    if (_m12_MC.GetMaximum ()> max) max = _m12_MC.GetMaximum ();
    _m12_data.SetMaximum(max/0.75);
    _m12_data.SetMinimum(0.001);
    TString xtitle = "M("+_name[0] +_name[1]+") GeV/c^{2}";
    _m12_data.GetXaxis()->SetTitle(xtitle);
    char title_12[200];
    sprintf(title_12, "Evetns/%.0f MeV", _binwidth*1E3);
    _m12_data.GetYaxis()->SetTitle(title_12);
    bes3plotstyle::FormatYAxis(_m12_data.GetYaxis());
    bes3plotstyle::FormatAxis(_m12_data.GetXaxis());
    _m12_data.Draw("E");
    _m12_bkg.Draw("same");
    _m12_MC.Draw("same");

    // TLegend *leg = new TLegend(0.22,0.7,0.45,0.8 );
    // leg->AddEntry(p11,"a(980)^{+} K_{S}");
    // leg->AddEntry(p12,"#bar{K}^{*}(1430)#pi^{+}");
    // leg->AddEntry(p22,"Interference");
    // Format(leg);
    // leg->Draw();

    Double_t chisqnum;
    Int_t nbinnum;
    Int_t isgood;
    _m12_data.Chi2TestX(&_m12_MC, chisqnum, nbinnum, isgood, "UU");
    TLatex lt;
    lt.SetNDC();
    lt.SetTextAngle(0);
    lt.SetTextSize(0.04);
    lt.DrawLatex(0.50,  0.90,  Form("#chi^{2}/nbin = %.2f/%d = %.2f",
                chisqnum, nbinnum,  chisqnum/nbinnum));
    pad_1->Update();
}
void PWAPlot::Plotm13() {
    FormatData();

    TPad *pad_2 = new TPad("pad_2", "", 0.0, 0.00,  1.0 ,  1.0);
    bes3plotstyle::Format(pad_2);
    pad_2->Draw();
    pad_2->cd();
    Double_t dataNum = _m13_data.Integral();
    Double_t bkgNum = _m13_bkg.Integral();
    Double_t mcNum = _m13_MC.Integral();
    _scaleFactor = (dataNum - bkgNum)/mcNum;
    _m13_MC.Scale(_scaleFactor);
    _m13_MC.Add(&_m13_bkg);

    Double_t max = _m13_data.GetMaximum();
    if (_m13_MC.GetMaximum ()> max) max = _m13_MC.GetMaximum ();
    _m13_data.SetMaximum(max/0.75);
    _m13_data.SetMinimum(0.001);

    char ytitle[200];
    sprintf(ytitle, "Evetns/%.0f MeV", _binwidth*1E3);
    _m13_data.GetYaxis()->SetTitle(ytitle);

    TString xtitle = "M("+_name[0] +_name[2]+") GeV/c^{2}";
    _m13_data.GetXaxis()->SetTitle(xtitle);

    bes3plotstyle::FormatYAxis(_m13_data.GetYaxis());
    bes3plotstyle::FormatAxis(_m13_data.GetXaxis());
    _m13_data.Draw("E");
    _m13_bkg.Draw("same");
    _m13_MC.Draw("same");

    // TLegend *leg = new TLegend(0.22,0.7,0.45,0.8 );
    // leg->AddEntry(p11,"a(980)^{+} K_{S}");
    // leg->AddEntry(p12,"#bar{K}^{*}(1430)#pi^{+}");
    // leg->AddEntry(p22,"Interference");
    // Format(leg);
    // leg->Draw();

    Double_t chisqnum;
    Int_t nbinnum;
    Int_t isgood;
    _m13_data.Chi2TestX(&_m13_MC, chisqnum, nbinnum, isgood, "UU");
    TLatex lt;
    lt.SetNDC();
    lt.SetTextAngle(0);
    lt.SetTextSize(0.04);
    lt.DrawLatex(0.50,  0.90,  Form("#chi^{2}/nbin = %.2f/%d = %.2f",
                chisqnum, nbinnum,  chisqnum/nbinnum));
    pad_2->Update();
}
void PWAPlot::Plotm23() {
    FormatData();

    TPad *pad_3 = new TPad("pad_3", "", 0.0, 0.00,  1.0 ,  1.0);
    bes3plotstyle::Format(pad_3);
    pad_3->Draw();
    pad_3->cd();
    Double_t dataNum = _m23_data.Integral();
    Double_t bkgNum = _m23_bkg.Integral();
    Double_t mcNum = _m23_MC.Integral();
    _scaleFactor = (dataNum - bkgNum)/mcNum;
    _m23_MC.Scale(_scaleFactor);
    _m23_MC.Add(&_m23_bkg);

    Double_t max = _m23_data.GetMaximum();
    if (_m23_MC.GetMaximum ()> max) max = _m23_MC.GetMaximum ();
    _m23_data.SetMaximum(max/0.75);
    _m23_data.SetMinimum(0.001);

    char ytitle[200];
    sprintf(ytitle, "Evetns/%.0f MeV", _binwidth*1E3);
    _m23_data.GetYaxis()->SetTitle(ytitle);

    TString xtitle = "M("+_name[1] +_name[2]+") GeV/c^{2}";
    _m23_data.GetXaxis()->SetTitle(xtitle);
    bes3plotstyle::FormatYAxis(_m23_data.GetYaxis());
    bes3plotstyle::FormatAxis(_m23_data.GetXaxis());

    _m23_data.Draw("E");
    _m23_bkg.Draw("same");
    _m23_MC.Draw("same");

    // TLegend *leg = new TLegend(0.22,0.7,0.45,0.8 );
    // leg->AddEntry(p11,"a(980)^{+} K_{S}");
    // leg->AddEntry(p12,"#bar{K}^{*}(1430)#pi^{+}");
    // leg->AddEntry(p22,"Interference");
    // Format(leg);
    // leg->Draw();

    Double_t chisqnum;
    Int_t nbinnum;
    Int_t isgood;
    _m23_data.Chi2TestX(&_m23_MC, chisqnum, nbinnum, isgood, "UU");
    TLatex lt;
    lt.SetNDC();
    lt.SetTextAngle(0);
    lt.SetTextSize(0.04);
    lt.DrawLatex(0.50,  0.90,  Form("#chi^{2}/nbin = %.2f/%d = %.2f",
                chisqnum, nbinnum,  chisqnum/nbinnum));
    pad_3->Update();
}
void PWAPlot::Plotdalitz() {
    FormatData();
    TPad *p_dali = new TPad("p_dali", "", 0.0, 0.00,  1.0 ,  1.0);
    bes3plotstyle::Format(p_dali);
    p_dali->Draw();
    p_dali->cd();
    Double_t dataNum = _dalitz_data.Integral();
    Double_t bkgNum = _dalitz_bkg.Integral();
    Double_t mcNum = _dalitz_MC.Integral();
    _dalitz_MC.Scale((dataNum - bkgNum)/mcNum);
    _dalitz_MC.Add(&_dalitz_bkg);

    TString xtitle = "M^{2}("+_name[0] +_name[2]+") GeV/c^{2}";
    TString ytitle = "M^{2}("+_name[1] +_name[2]+") GeV/c^{2}";
    _dalitz_MC.GetXaxis()->SetTitle(xtitle);
    _dalitz_data.GetXaxis()->SetTitle(xtitle);
    _dalitz_MC.GetYaxis()->SetTitle(ytitle);
    _dalitz_data.GetYaxis()->SetTitle(ytitle);
    bes3plotstyle::FormatYAxis(_dalitz_MC.GetYaxis());
    bes3plotstyle::FormatAxis(_dalitz_MC.GetXaxis());
    bes3plotstyle::FormatYAxis(_dalitz_data.GetYaxis());
    bes3plotstyle::FormatAxis(_dalitz_data.GetXaxis());
    _dalitz_data.SetMarkerSize(0.001);

    Double_t chisqnum;
    Int_t nbinnum;
    Int_t isgood;
    _dalitz_data.Chi2TestX(&_dalitz_MC, chisqnum, nbinnum, isgood, "UU");
    _dalitz_data.Add(&_dalitz_MC, -1);
    // diff.Draw("contz");
    _dalitz_data.Draw("colz");
    // _dalitz_bkg.Draw("same");

    // TLegend *leg = new TLegend(0.22,0.7,0.45,0.8 );
    // leg->AddEntry(p11,"a(980)^{+} K_{S}");
    // leg->AddEntry(p12,"#bar{K}^{*}(1430)#pi^{+}");
    // leg->AddEntry(p22,"Interference");
    // Format(leg);
    // leg->Draw();

    TLatex lt;
    lt.SetNDC();
    lt.SetTextAngle(0);
    lt.SetTextSize(0.04);
    lt.DrawLatex(0.50,  0.90,  Form("#chi^{2}/nbin = %.2f/%d = %.2f",
                chisqnum, nbinnum,  chisqnum/nbinnum));
    p_dali->Update();
}
void PWAPlot::PlotData(TString title) {
    FormatData();
    TCanvas c("c", "",  800,  600);
    TPad *p_dali = new TPad("p_data", "", 0.0, 0.00,  1.0 ,  1.0);
    bes3plotstyle::Format(p_dali);
    p_dali->Draw();
    p_dali->cd();

    TString xtitle = "M^{2}("+_name[0] +_name[2]+") GeV/c^{2}";
    TString ytitle = "M^{2}("+_name[1] +_name[2]+") GeV/c^{2}";
    _dalitz_data.GetXaxis()->SetTitle(xtitle);
    _dalitz_data.GetYaxis()->SetTitle(ytitle);
    bes3plotstyle::FormatYAxis(_dalitz_data.GetYaxis());
    bes3plotstyle::FormatAxis(_dalitz_data.GetXaxis());
    _dalitz_data.SetMarkerSize(0.001);
    _dalitz_data.SetMinimum(0);

    _dalitz_data.Draw("colz");
    p_dali->Update();
    p_dali->Draw();
    c.SaveAs(title+".dataDalitz.eps");
    delete p_dali;
}
void PWAPlot::PlotFitResult(TString title) {
    FormatData();
    TCanvas c("c", "", 800, 600);
    TPad *p_dali = new TPad("p_fitresult", "", 0.0, 0.00,  1.0 ,  1.0);
    bes3plotstyle::Format(p_dali);
    p_dali->Draw();
    p_dali->cd();
    Double_t dataNum = _dalitz_data.Integral();
    Double_t bkgNum = _dalitz_bkg.Integral();
    Double_t mcNum = _dalitz_MC.Integral();
    // _dalitz_MC.Scale();

    TString xtitle = "M^{2}("+_name[0] +_name[2]+") GeV/c^{2}";
    TString ytitle = "M^{2}("+_name[1] +_name[2]+") GeV/c^{2}";
    _dalitz_MC.GetXaxis()->SetTitle(xtitle);
    // _dalitz_data.GetXaxis()->SetTitle(xtitle);
    _dalitz_MC.GetYaxis()->SetTitle(ytitle);
    // _dalitz_data.GetYaxis()->SetTitle(ytitle);
    bes3plotstyle::FormatYAxis(_dalitz_MC.GetYaxis());
    bes3plotstyle::FormatAxis(_dalitz_MC.GetXaxis());
    // bes3plotstyle::FormatYAxis(_dalitz_data.GetYaxis());
    // bes3plotstyle::FormatAxis(_dalitz_data.GetXaxis());
    _dalitz_MC.SetMarkerSize(0.001);

    // diff.Draw("contz");
    _dalitz_MC.Draw("colz");
    // _dalitz_bkg.Draw("same");

    // TLegend *leg = new TLegend(0.22,0.7,0.45,0.8 );
    // leg->AddEntry(p11,"a(980)^{+} K_{S}");
    // leg->AddEntry(p12,"#bar{K}^{*}(1430)#pi^{+}");
    // leg->AddEntry(p22,"Interference");
    // Format(leg);
    // leg->Draw();
    p_dali->Update();
    p_dali->Draw();
    c.SaveAs(title+".fitResult.eps");
    delete p_dali;
}
void PWAPlot::PlotEfficiency(TString name) {
    TFile *f = new TFile(_projectfile, "read");
    TTree *tree = reinterpret_cast<TTree*>(f->Get("project"));
    TH2D *h2 = new TH2D("efficiency", "",
            _bins, pow(_m13_bounds[0], 2), pow(_m13_bounds[1], 2),
            _bins, pow(_m23_bounds[0], 2), pow(_m23_bounds[1], 2) );
    tree->Project("efficiency", "m23**2:m13**2");
    TString xtitle = "M^{2}("+_name[0] +_name[2]+") GeV^{2}";
    TString ytitle = "M^{2}("+_name[1] +_name[2]+") GeV^{2}";
    h2->GetXaxis()->SetTitle(xtitle);
    h2->GetYaxis()->SetTitle(ytitle);
    bes3plotstyle::FormatYAxis(h2->GetYaxis());
    bes3plotstyle::FormatAxis(h2->GetXaxis());
    bes3plotstyle::SetStyle();
    TCanvas *cEff = new TCanvas("cEff", "",  800,  600);
    h2->Draw("box");
    cEff->SaveAs(name+".efficiency.eps");
    delete h2;
    delete cEff;
    delete tree;
    f->Close();
    delete f;
}
void PWAPlot::Draw(TString name) {
    cout << "Inf: read data..." << endl;
    ReadData();
    cout << "Inf: read fit result project..." << endl;
    ReadMC();
    cout << "Inf: Set bes3plotstyle..." << endl;
    bes3plotstyle::SetStyle();
    PlotData(name);
    PlotFitResult(name);
    TCanvas c("c", "c",  800,  600);
    c.Divide(2, 2);
    c.cd(1);
    Plotm12();
    c.Update();
    c.cd(2);
    Plotm13();
    c.Update();
    c.cd(3);
    Plotm23();
    c.Update();
    c.cd(4);
    Plotdalitz();
    c.Update();
    cout << "Inf: Save the figure:" << name << ".compare4.eps" << endl;
    c.SaveAs(name+".compare4.eps");
    PlotEfficiency(name);
    SetColor();
    PlotWavesM12(name);
    PlotWavesM13(name);
    PlotWavesM23(name);
}
void PWAPlot::PlotLengend(const std::vector<TH1D>& hists) {
    TLegend *leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    for (int i = 0; i< 50; ++i) {
        if (m_ResNames[i] != "") {
            leg->AddEntry(&hists[i], m_ResNames[i]);
        }
    }
    bes3plotstyle::Format(leg);
    leg->Draw();
}
void PWAPlot::PlotWavesM12(TString title) {
    Int_t waves = GetWaves();
    TString *hist_name = new TString[waves];
    std::vector<TH1D> MC_histos;
    if (!_isSetBound) SetBound ();
    for (Int_t i  = 0 ; i < waves; i++) {
        hist_name[i] = "_mc_"+str(i);
        TH1D h(hist_name[i], "",  _bins,  _m12_bounds[0],  _m12_bounds[1]);
        MC_histos.push_back(h);
    }

    TFile *f = new TFile(_projectfile, "read");
    TTree *tree = reinterpret_cast<TTree*>(f->Get("project"));
    Double_t m12;
    Double_t FF[900];
    tree->SetBranchAddress("m12", &m12);
    tree->SetBranchAddress(m_FFName, FF);
    Double_t eva;
    tree->SetBranchAddress("eva", &eva);

    for (Int_t i = 0 ; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        for (Int_t kk = 0 ; kk < waves; kk++) {
            MC_histos[kk].Fill(m12, FF[kk*waves+kk]);
           // cout<<"FF["<<kk<<"]"<<FF[kk*waves+kk]<<endl;
        }
        // cout<<"Eva: "<<eva<<endl;
    }
    FormatData();
    delete tree;
    delete f;
    TCanvas c("c_com", "", 800, 600);
    TPad *pad_1 = new TPad("pad_n", "", 0.0, 0.00,  1.0 ,  1.0);
    bes3plotstyle::Format(pad_1);
    pad_1->Draw();
    pad_1->cd();

    Double_t dataNum = _m12_data.Integral();
    Double_t bkgNum = _m12_bkg.Integral();
    Double_t mcNum = _m12_MC.Integral();
    Double_t scale_factor = (dataNum - bkgNum) / mcNum;
    _m12_MC.Scale(scale_factor);
    _m12_MC.Add(&_m12_bkg);

    Double_t max = _m12_data.GetMaximum();
    if (_m12_MC.GetMaximum ()> max) max = _m12_MC.GetMaximum ();
    // _m12_data.SetMaximum(max/0.85);
    _m12_data.SetMinimum(0.001);
    TString xtitle = "M("+_name[0] +_name[1]+") GeV/c^{2}";
    _m12_data.GetXaxis()->SetTitle(xtitle);
    char title_12[200];
    sprintf(title_12, "Evetns/%.0f MeV", _binwidth*1E3);
    _m12_data.GetYaxis()->SetTitle(title_12);
    bes3plotstyle::FormatYAxis(_m12_data.GetYaxis());
    bes3plotstyle::FormatAxis(_m12_data.GetXaxis());
    _m12_data.Draw("E");
    _m12_bkg.Draw("same");
    _m12_MC.Draw("same");
    Double_t chisqnum;
    Int_t nbinnum;
    Int_t isgood;
    _m12_data.Chi2TestX(&_m12_MC, chisqnum, nbinnum, isgood, "UU");
    TLatex lt;
    lt.SetNDC();
    lt.SetTextAngle(0);
    lt.SetTextSize(0.04);
    lt.DrawLatex(0.50,  0.90,  Form("#chi^{2}/nbin = %.2f/%d = %.2f",
                chisqnum, nbinnum,  chisqnum/nbinnum));
    // pad_1->Update();
    for (Int_t kk = 0;  kk < MC_histos.size(); kk++) {
        MC_histos[kk].Scale(_scaleFactor);
        MC_histos[kk].SetLineColor(2+kk);
        MC_histos[kk].SetLineStyle(2+kk);
        MC_histos[kk].SetLineWidth(2);
        MC_histos[kk].SetMarkerSize(0);
        MC_histos[kk].Draw("same");
    }
    PlotLengend(MC_histos);
    pad_1->Update();
    c.SaveAs(title+".projectm12.eps");
    delete pad_1;
    delete [] hist_name;
}
TString PWAPlot::str(Int_t kk) {
    std::stringstream s;
    s << kk;
    return s.str();
}

void PWAPlot::PlotWavesM13(TString title) {
    Int_t waves = GetWaves();
    TString *hist_name = new TString[waves];
    std::vector<TH1D> MC_histos;
    if (!_isSetBound) SetBound ();
    for (Int_t i  = 0 ; i < waves; i++) {
        hist_name[i] = "_mc_"+str(i);
        TH1D h(hist_name[i], "",  _bins,  _m13_bounds[0],  _m13_bounds[1]);
        MC_histos.push_back(h);
    }

    TFile *f = new TFile(_projectfile, "read");
    TTree *tree = reinterpret_cast<TTree*>(f->Get("project"));
    Double_t m13;
    Double_t FF[900];
    tree->SetBranchAddress("m13", &m13);
    tree->SetBranchAddress(m_FFName, FF);
    Double_t eva;
    tree->SetBranchAddress("eva", &eva);

    for (Int_t i = 0 ; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        for (Int_t kk = 0 ; kk < waves; kk++) {
            MC_histos[kk].Fill(m13, FF[kk*waves+kk]);
           // cout<<"FF["<<kk<<"]"<<FF[kk*waves+kk]<<endl;
        }
        // cout<<"Eva: "<<eva<<endl;
    }
    FormatData();
    delete tree;
    delete f;
    TCanvas c("c_com", "", 800, 600);
    TPad *pad_1 = new TPad("pad_n", "", 0.0, 0.00,  1.0 ,  1.0);
    bes3plotstyle::Format(pad_1);
    pad_1->Draw();
    pad_1->cd();

    Double_t dataNum = _m13_data.Integral();
    Double_t bkgNum = _m13_bkg.Integral();
    Double_t mcNum = _m13_MC.Integral();
    Double_t scale_factor = (dataNum - bkgNum) / mcNum;
    _m13_MC.Scale(scale_factor);
    _m13_MC.Add(&_m13_bkg);

    Double_t max = _m13_data.GetMaximum();
    if (_m13_MC.GetMaximum ()> max) max = _m13_MC.GetMaximum ();
    // _m13_data.SetMaximum(max/0.85);
    _m13_data.SetMinimum(0.001);
    TString xtitle = "M("+_name[0] +_name[2]+") GeV/c^{2}";
    _m13_data.GetXaxis()->SetTitle(xtitle);
    char title_12[200];
    sprintf(title_12, "Evetns/%.0f MeV", _binwidth*1E3);
    _m13_data.GetYaxis()->SetTitle(title_12);
    bes3plotstyle::FormatYAxis(_m13_data.GetYaxis());
    bes3plotstyle::FormatAxis(_m13_data.GetXaxis());
    _m13_data.Draw("E");
    _m13_bkg.Draw("same");
    _m13_MC.Draw("same");
    Double_t chisqnum;
    Int_t nbinnum;
    Int_t isgood;
    _m13_data.Chi2TestX(&_m13_MC, chisqnum, nbinnum, isgood, "UU");
    TLatex lt;
    lt.SetNDC();
    lt.SetTextAngle(0);
    lt.SetTextSize(0.04);
    lt.DrawLatex(0.50,  0.90,  Form("#chi^{2}/nbin = %.2f/%d = %.2f",
                chisqnum, nbinnum,  chisqnum/nbinnum));
    // pad_1->Update();
    for (Int_t kk = 0;  kk < MC_histos.size(); kk++) {
        MC_histos[kk].Scale(_scaleFactor);
        MC_histos[kk].SetLineColor(2+kk);
        MC_histos[kk].SetLineStyle(2+kk);
        MC_histos[kk].SetLineWidth(2);
        MC_histos[kk].SetMarkerSize(0);
        MC_histos[kk].Draw("same");
    }
    PlotLengend(MC_histos);
    pad_1->Update();
    c.SaveAs(title+".projectm13.eps");
    delete pad_1;
    delete [] hist_name;
}

void PWAPlot::PlotWavesM23(TString title) {
    Int_t waves = GetWaves();
    TString *hist_name = new TString[waves];
    std::vector<TH1D> MC_histos;
    if (!_isSetBound) SetBound ();
    for (Int_t i  = 0 ; i < waves; i++) {
        hist_name[i] = "_mc_"+str(i);
        TH1D h(hist_name[i], "",  _bins,  _m23_bounds[0],  _m23_bounds[1]);
        MC_histos.push_back(h);
    }

    TFile *f = new TFile(_projectfile, "read");
    TTree *tree = reinterpret_cast<TTree*>(f->Get("project"));
    Double_t m23;
    Double_t FF[900];
    tree->SetBranchAddress("m23", &m23);
    tree->SetBranchAddress(m_FFName, FF);
    Double_t eva;
    tree->SetBranchAddress("eva", &eva);

    for (Int_t i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        for (Int_t kk = 0 ; kk < waves; kk++) {
            MC_histos[kk].Fill(m23, FF[kk*waves+kk]);
           // cout<<"FF["<<kk<<"]"<<FF[kk*waves+kk]<<endl;
        }
        // cout<<"Eva: "<<eva<<endl;
    }
    FormatData();
    delete tree;
    delete f;
    TCanvas c("c_com", "", 800, 600);
    TPad *pad_1 = new TPad("pad_n", "", 0.0, 0.00,  1.0 ,  1.0);
    bes3plotstyle::Format(pad_1);
    pad_1->Draw();
    pad_1->cd();

    Double_t dataNum = _m23_data.Integral();
    Double_t bkgNum = _m23_bkg.Integral();
    Double_t mcNum = _m23_MC.Integral();
    Double_t scale_factor = (dataNum - bkgNum) / mcNum;
    _m23_MC.Scale(scale_factor);
    _m23_MC.Add(&_m23_bkg);

    Double_t max = _m23_data.GetMaximum();
    if (_m23_MC.GetMaximum ()> max) max = _m23_MC.GetMaximum ();
    // _m23_data.SetMaximum(max/0.85);
    _m23_data.SetMinimum(0.001);
    TString xtitle = "M("+_name[1] +_name[2]+") GeV/c^{2}";
    _m23_data.GetXaxis()->SetTitle(xtitle);
    char title_12[200];
    sprintf(title_12, "Evetns/%.0f MeV", _binwidth*1E3);
    _m23_data.GetYaxis()->SetTitle(title_12);
    bes3plotstyle::FormatYAxis(_m23_data.GetYaxis());
    bes3plotstyle::FormatAxis(_m23_data.GetXaxis());
    _m23_data.Draw("E");
    _m23_bkg.Draw("same");
    _m23_MC.Draw("same");
    Double_t chisqnum;
    Int_t nbinnum;
    Int_t isgood;
    _m23_data.Chi2TestX(&_m23_MC, chisqnum, nbinnum, isgood, "UU");
    TLatex lt;
    lt.SetNDC();
    lt.SetTextAngle(0);
    lt.SetTextSize(0.04);
    lt.DrawLatex(0.50,  0.90,  Form("#chi^{2}/nbin = %.2f/%d = %.2f",
                chisqnum, nbinnum,  chisqnum/nbinnum));
    // pad_1->Update();
    for (Int_t kk = 0;  kk < MC_histos.size(); kk++) {
        MC_histos[kk].Scale(_scaleFactor);
        MC_histos[kk].SetLineColor(2+kk);
        MC_histos[kk].SetLineStyle(2+kk);
        MC_histos[kk].SetLineWidth(2);
        MC_histos[kk].SetMarkerSize(0);
        MC_histos[kk].Draw("same");
    }
    PlotLengend(MC_histos);
    pad_1->Update();
    c.SaveAs(title+".projectm23.eps");
    delete pad_1;
    delete [] hist_name;
}
void PWAPlot::SetColor() {
    gROOT->GetColor(2)->SetRGB(0.70 , 0.16 , 0.14);
    gROOT->GetColor(3)->SetRGB(0.20 , 0.60 , 0.20);
    gROOT->GetColor(4)->SetRGB(0.10 , 0.13 , 0.77);
    gROOT->GetColor(5)->SetRGB(0.40 , 0.20 , 0.40);
    gROOT->GetColor(6)->SetRGB(0.60 , 0.17 , 0.23);
    gROOT->GetColor(7)->SetRGB(0.23 , 0.17 , 0.60);
    gROOT->GetColor(8)->SetRGB(0.23 , 0.60 , 0.17);
    gROOT->GetColor(9)->SetRGB(0.40 , 0.40 , 0.20);
    gROOT->GetColor(10)->SetRGB(0.20 , 0.40 , 0.20);
}
void PWAPlot::SetParName(const TString & name1,
        const TString &name2,
        const TString &name3) {
    _name[0] = name1;
    _name[1] = name2;
    _name[2] = name3;
}

void PWAPlot::SetResName(const int& index, const TString &name ) {
    if (index <=0) {
        cout << "Error::PWAPlot::SetResName(), the index must be larger than 1"
             << endl;
    }
    if (index > GetWaves()) {
        cout<< "Error:: PWAPlot::SetResName():"
             <<   "The total waves is: " << GetWaves() << endl;
    }
    m_ResNames[index-1] = name;
}
