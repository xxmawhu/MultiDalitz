// Copyright (c) 2019-6-24 maxx
/////////////////////////////////////////////////////////////////////////
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "TCanvas.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include <iostream>
#include <map>
#include <string>

// use this order for safety on library loading
using namespace RooFit;
using namespace RooStats;
using std::cout;
using std::endl;
using std::map;
using std::string;

void AddRes(TString, map<string, int> &);

void Fit(int kk = 1, TString resonanceList  =  "",
        TString Datadat = "/besfs/users/maxx/commu/hsuhz/Sample/Data/Data.dat",
        TString Fracdat = "/besfs/users/maxx/commu/hsuhz/Sample/PHSP/Frac.dat",
        TString PHSPdat = "/besfs/users/maxx/commu/hsuhz/Sample/PHSP/PHSP.dat"
        ) {
    // Fracdat used for compute the partial fraction

    gSystem->Load("libPhysics");
    gSystem->Load("/besfs/users/maxx/commu/hsuhz/MultiDalitz/lib/libKsKsK.so");

    map<string, int> AddResonance;
    AddRes(resonanceList, AddResonance);
    // model
    double high = 1.8865;
    double low = 0-high;
    const double pi = 3.1415926*4;
    RooRealVar v11("v11", "v11", low, high);
    RooRealVar v12("v12", "v12", low, high);
    RooRealVar v13("v13", "v13", low, high);
    RooRealVar v14("v14", "v14", low, high);
    RooRealVar v21("v21", "v21", low, high);
    RooRealVar v22("v22", "v22", low, high);
    RooRealVar v23("v23", "v23", low, high);
    RooRealVar v24("v24", "v24", low, high);
    RooRealVar v31("v31", "v31", low, high);
    RooRealVar v32("v32", "v32", low, high);
    RooRealVar v33("v33", "v33", low, high);
    RooRealVar v34("v34", "v34", low, high);

    MultiDalitPdf ksksKpdf("ksksKpdf", "ksksKpdf",
            v11, v12, v13, v14, v21, v22, v23, v24, v31, v32, v33, v34,
            PHSPdat);
    ksksKpdf.setFracDat(Fracdat);

    // -----------------Orbital angular momentum-------------------
    RooRealVar spin0("spin0", "spin0",  0);
    RooRealVar spin1("spin1", "spin1",  1);
    RooRealVar spin2("spin2", "spin2",  2);
    const int S_WAVE(0);
    const int P_WAVE(1);
    spin0.setConstant();
    spin1.setConstant();
    spin2.setConstant();
    RooRealVar rRes("rRes", "", 3);
    rRes.setConstant();

    // http:// pdglive.lbl.gov/Particle.action?init=0&node=M027&home=MXXX005
    RooRealVar meta1405("meta1405", "meta1405", 1413.9E-3);
    meta1405.setConstant();
    RooRealVar weta1405("weta1405", "meta1405", 48E-3);
    weta1405.setConstant();
    ksksKpdf.configMother("eta(1405)", "eta(1405)", LineShape::RBW,
            RooArgList(spin0, rRes, meta1405, weta1405));

    // http:// pdglive.lbl.gov/Particle.action?init=0&node=M175&home=MXXX005
    RooRealVar meta1475("meta1475", "meta1475", 1475E-3);
    meta1475.setConstant();
    RooRealVar weta1475("weta1475", "weta1475", 90E-3, 10e-3, 200e-3);
    //weta1475.setConstant();
    ksksKpdf.configMother("eta(1475)", "eta(1475)", LineShape::RBW,
            RooArgList(spin0, rRes, meta1475, weta1475));

    // ----------------add partial wave----------------------------
    if (AddResonance["a0_980"]) {
        RooRealVar ma0("ma0",    "ma0" , 0.990 , 0.88,   1.08);
        ma0.setConstant();
        RooRealVar gPiEtaa0 ("gPiEtaa0", "",  sqrt(0.341));
        gPiEtaa0.setConstant();
        RooRealVar gKK2a0 ("gKK2a0", "",  sqrt(0.892*0.341));
        gKK2a0.setConstant();
        RooRealVar Rhoa0_980("Rhoa0_980", "",  1,  0,  10);
        RooRealVar Phia0_980("Phia0_980", "",  0,  0,  15);
        ksksKpdf.addResonance("a0_980", "eta(1405)",
                LineShape::a980_0, DecayType::AB,
                Rhoa0_980, Phia0_980,
                RooArgList(spin0, rRes, ma0, gPiEtaa0 , gKK2a0), S_WAVE);
    }
    if (AddResonance["K*+"]) {
        RooRealVar mkstp("mkstp",    "mkstp" , 0.89176);
        mkstp.setConstant();
        RooRealVar widkstp ("widkstp", "",  50.3e-3);
        widkstp.setConstant();
        RooRealVar r_kstp("r_kstp", "",  3.0);
        r_kstp.setConstant();
        RooRealVar Rhokstp("Rhokstp", "",  10,  0,  20);
        Rhokstp.setConstant();
        RooFormulaVar Rhokstm("Rhokstm", "", "-Rhokstp",  Rhokstp);
        RooRealVar Phikstp("Phikstp", "",  0.3,  0,  15);
        Phikstp.setConstant();
        ksksKpdf.addResonance("K*+", "eta(1475)",
                LineShape::RBW, DecayType::AC,
                Rhokstp, Phikstp,
                RooArgList(spin1, r_kstp, mkstp, widkstp), P_WAVE);
        // add C-conjuge wave
        ksksKpdf.addResonance("K*-", "eta(1475)",
                LineShape::RBW, DecayType::BC,
                Rhokstm, Phikstp,
                RooArgList(spin1, r_kstp, mkstp, widkstp), P_WAVE);
    }
    if (AddResonance["kappa"]) {
        RooRealVar mkappap("mkappap",    "mkappap" , 0.824);
        mkappap.setConstant();
        RooRealVar widkappap ("widkappap", "",  478e-3);
        widkappap.setConstant();
        RooRealVar r_kappap("r_kappap", "",  3.0);
        r_kappap.setConstant();
        RooRealVar Rhokappap("Rhokappap", "",  1,  0,  20);
        RooFormulaVar Rhokappam("Rhokappam", "", "-Rhokappap",  Rhokappap);
        RooRealVar Phikappap("Phikappap", "",  0,  0,  15);
        ksksKpdf.addResonance("kappa+", "eta(1475)",
                LineShape::RBW, DecayType::AC,
                Rhokappap, Phikappap,
                RooArgList(spin0, r_kappap, mkappap, widkappap), S_WAVE);
        // add C-conjuge wave
        ksksKpdf.addResonance("kappa-", "eta(1475)",
                LineShape::RBW, DecayType::BC,
                Rhokappam, Phikappap,
                RooArgList(spin0, r_kappap, mkappap, widkappap), S_WAVE);
    }
    if (AddResonance["NoneRes"]) {
        RooRealVar mNoneRes("mNoneRes",  "mNoneRes" , 0.982 , 0.88,   1.08);
        mNoneRes.setConstant();
        RooRealVar GNoneRes("GNoneRes", "GNoneRes",  47.4E-3,   10E-3,   200E-3);
        GNoneRes.setConstant();
        RooRealVar RhoNoneRes("RhoNoneRes", "",  1,  0,  10);
        RooRealVar PhiNoneRes("PhiNoneRes", "",  0,  0,  15);

        ksksKpdf.addResonance("NoneResonance", "eta(1475)",
                LineShape::Flat,  DecayType::AB,
                RhoNoneRes, PhiNoneRes,
                RooArgList(spin0), S_WAVE);
    }

    // ksksKpdf.Print();
    // cout << "======== test ========" << endl;
    // ksksKpdf.test();
    // cout << "========= end of test ======" << endl;
    RooRealVar weight("weight", "weight",  -2,  2);

    RooArgSet theSet1, theSet2, theSet3;
    theSet1.add(RooArgSet(v11, v12, v13, v14, v21, v22, v23, v24));
    theSet2.add(RooArgSet(v31, v32, v33, v34));
    theSet3.add(RooArgSet(weight));
    RooArgSet theSet4(theSet1, theSet2, "");
    RooArgSet theSet(theSet4, theSet3, "");

    RooDataSet *data = RooDataSet::read(Datadat, theSet);
    data->Print();

    RooDataSet *datWeight = new RooDataSet(data->GetName(), data->GetTitle(),
                data,*data->get(), 0, weight.GetName());
    datWeight->Print();

    // RooFitResult *result = ksksKpdf.fitTo(*datEtaV2G , Save(kTRUE));
    // return ;
    RooFitResult *result = ksksKpdf.fitTo(*datWeight, Save(kTRUE));
    result->Print();
    ksksKpdf.test();
    TString name = "fitresult";
    ksksKpdf.fitFractions(cout);
    ksksKpdf.project("project.root");
    /*
    result->Write("result");
    ksksKpdf.Write("pdf");

    f.Close();
    char pronm[200];
    sprintf(pronm, "project_%d.root",kk);
    ofstream fout("FCN.dat",ios::app);
    fout<<"FCN:\t"<<result->minNll()<<endl;
    fout.close();
    result->Print();
    ofstream fracout("Fraction.dat",ios::app);
    ksksKpdf.fitFractions(fitparameter,1, fracout);
    PWAPlot plot(pronm, Datadat);
    plot.Draw();
    */
}

void AddRes(TString resonanceList , map<string, int> &AddResonance) {
    AddResonance["a0_980"] = 0; 
    AddResonance["NoneRes"] = 1;
    AddResonance["K*+"] = 1;
    AddResonance["kappa"] = 1;
    if (resonanceList != "") {
        TMVA::Tools::Instance();
        vector<TString> mlist = TMVA::gTools().SplitString(resonanceList, ',');

        for (int i=0; i < mlist.size(); i++) {
            string res(mlist[i]);

            if (AddResonance.find(res) == AddResonance.end()) {
                cout << "The resonance \"" << res
                    << "\" is not defined yed, choose from those" << endl;
                for (map<string, int>::iterator it = AddResonance.begin (); it != AddResonance.end (); it++) {
                    std::cout << it->first << " ";
                }
                std::cout << std::endl;
            }
            else {
                cout << "Add " << res << " successful!!!" << endl;}
            AddResonance[res] = 1;
        }
    }
}
