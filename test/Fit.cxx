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
using namespace RooFit ;
using namespace RooStats ;
using std::cout;
using std::endl;
using std::map;
using std::string;

void AddRes(TString, map<string, int> &);

void Fit(int kk=1, TString resonanceList = "", 
        TString Datadat = "Data.dat",
        TString Fracdat = "Frac.dat",
        TString PHSPdat = "PHSP.dat"
        ){

    gSystem->Load("libPhysics");
    gSystem->Load("../lib/libKsKsK.so");

    map<string,int> AddResonance;
    AddRes(resonanceList, AddResonance);
    // model 
    double high = 1.8865;
    double low = 0-high;
    const double pi = 3.1415926*4;
    RooRealVar v11("v11","v11",low,high);
    RooRealVar v12("v12","v12",low,high);
    RooRealVar v13("v13","v13",low,high);
    RooRealVar v14("v14","v14",low,high);
    RooRealVar v21("v21","v21",low,high);
    RooRealVar v22("v22","v22",low,high);
    RooRealVar v23("v23","v23",low,high);
    RooRealVar v24("v24","v24",low,high);
    RooRealVar v31("v31","v31",low,high);
    RooRealVar v32("v32","v32",low,high);
    RooRealVar v33("v33","v33",low,high);
    RooRealVar v34("v34","v34",low,high);

    MultiDalitPdf ksksKpdf("ksksKpdf","ksksKpdf",   
            v11,v12,v13,v14,
            v21,v22,v23,v24,    
            v31,v32,v33,v34, 
            PHSPdat);
    ksksKpdf.setFracDat(Fracdat);
    ksksKpdf.setIdenticalParticle(12);
    ksksKpdf.setrD(5);
    //-----------------Orbital angular momentum-------------------
    RooRealVar spin0("spin0","spin0", 0);
    RooRealVar spin1("spin1","spin1", 1);
    RooRealVar spin2("spin2","spin2", 2);
    const int S_WAVE(0);
    const int P_WAVE(1);
    spin0.setConstant();
    spin1.setConstant();
    spin2.setConstant();
    RooRealVar rRes("rRes", "", 3);
    rRes.setConstant();
    RooRealVar meta("meta", "meta", 1.405, 1.0, 1.405);
    RooRealVar weta("weta", "weta", 1.405, 1.0, 1.405);

    ksksKpdf.configMother("eta(1405)", "eta(1405)", LineShape::RBW, 
            RooArgList(spin0, rRes, meta, weta));
    //----------------add partial wave----------------------------   
    if(AddResonance["a0_980"]){
        RooRealVar ma0("ma0",    "ma0" , 0.990 , 0.88,   1.08);
        ma0.setConstant();
        RooRealVar gPiEtaa0 ("gPiEtaa0","", 0.348);
        gPiEtaa0.setConstant();
        RooRealVar gKK2a0 ("gKK2a0","", 0.348);
        gKK2a0.setConstant();
        RooRealVar Rhoa0_980("Rhoa0_980","", 1, 0, 10);
        RooRealVar Phia0_980("Phia0_980","", 0, 0, 15);
        ksksKpdf.addResonance("a0_980", "eta(1405)", 
                LineShape::a980_p, DecayType::BC, 
                Rhoa0_980, Phia0_980,
                RooArgList(spin0, rRes, ma0, gPiEtaa0 , gKK2a0), P_WAVE);
    }
    if(AddResonance["rho_1700"]){
        RooRealVar mrho_1700("mrho_1700",    "mrho_1700" , 1.594);
        mrho_1700.setConstant();
        RooRealVar widrho_1700 ("widrho_1700","", 0.259);
        widrho_1700.setConstant();
        RooRealVar r_rho_1700("r_rho_1700","", 3.0);
        r_rho_1700.setConstant();
        RooRealVar Rhorho_1700("Rhorho_1700","", 1, 0, 10);
        RooRealVar Phirho_1700("Phirho_1700","", 0, 0, 15);
        ksksKpdf.addResonance("rho_1700", 
                LineShape::GS, DecayType::BC, 
                Rhorho_1700, Phirho_1700,
                RooArgList(spin1, r_rho_1700, mrho_1700, widrho_1700));
    }
    if(AddResonance["NoneRes"]){
        RooRealVar mNoneRes("mNoneRes",  "mNoneRes" , 0.982 , 0.88,   1.08);
        mNoneRes.setConstant();
        RooRealVar GNoneRes("GNoneRes","GNoneRes", 47.4E-3,  10E-3,  200E-3);
        GNoneRes.setConstant();
        RooRealVar RhoNoneRes("RhoNoneRes","", 1, 0, 10);
        RooRealVar PhiNoneRes("PhiNoneRes","", 0, 0, 15);

        ksksKpdf.addResonance("NoneResonance",  
                LineShape::Flat,  DecayType::AB,
                RhoNoneRes, PhiNoneRes,
                RooArgList(spin0));
    }
    /*
    if(AddResonance["pi1_1400"]){
        // pi1(1400)
        RooRealVar mpi1_1400("mpi1_1400",    "mpi1_1400" , 1.354 , 1.3,   1.4);
        mpi1_1400.setConstant();
        RooRealVar Gpi1_1400("Gpi1_1400",  "Gpi1_1400" , 330E-3,  200E-3,  500E-3);
        Gpi1_1400.setConstant();

        RooRealVar rhopi1_1400("rhopi1_1400","rhopi1_1400", 0.5, 0.0, 10.0);
        RooRealVar phipi1_1400("phipi1_1400","phipi1_1400", 0.0 ,-pi,pi);
        ksksKpdf.addResonance("Dtopi1_1400", "D+ to pi_1(1400)+ Ks", spin1, mpi1_1400,  Gpi1_1400,  
                rhopi1_1400, phipi1_1400, modetype23, r0, r1, g1);
    }
    if(AddResonance["a2_1320"]){
        // a2(1320)
        RooRealVar ma2_1320("ma2_1320", "ma2_1320", 1.3181, 1.2, 1.5);
        ma2_1320.setConstant();
        RooRealVar Ga2_1320("Ga2_1320", "Ga2_1320", 109.8E-3, 100E-3,  150E-3);
        Ga2_1320.setConstant();
        RooRealVar rhoa2_1320("rhoa2_1320","rhoa2_1320", 2.5, 0.0, 10.0);
        RooRealVar phia2_1320("phia2_1320","phia2_1320", 0.0 ,-pi,pi);
        fitparameter.add(rhoa2_1320);
        fitparameter.add(phia2_1320);
        ksksKpdf.addResonance("Dtoa2_1320", "Dto a_2(1320)+ Ks", spin2,
                ma2_1320,  Ga2_1320,  rhoa2_1320, phia2_1320,
                modetype23, r0, r1, g1); 
    }
    if(AddResonance["a0_1450"]){
        // a0(1450)
        RooRealVar ma0_1450("ma0_1450",    "ma0_1450" , 1.474 , 1.2,   1.5);
        ma0_1450.setConstant();
        RooRealVar Ga0_1450("Ga0_1450",  "Ga0_1450" , 265E-3,  150E-3,  350E-3);
        Ga0_1450.setConstant();
        RooRealVar rhoa0_1450("rhoa0_1450","rhoa0_1450", 1.5, 0.0, 10.0);
        RooRealVar phia0_1450("phia0_1450","phia0_1450", 0.0 ,-pi,pi);
        fitparameter.add(rhoa0_1450);
        fitparameter.add(phia0_1450);
        ksksKpdf.addResonance("Dtoa0_1450", "Dto a_0(1450)+ Ks", spin0,
                ma0_1450,  Ga0_1450,  rhoa0_1450, phia0_1450,
                modetype23, r0, r1, g1); 
    }

    if(AddResonance["Pion_Ks_P"]){
        // ( pi+ Ks ) P wave
        RooRealVar mPion_Ks_P("mPion_Ks_P",  "mPion_Ks_P", 1.60 , 1.55,   1.65);
        mPion_Ks_P.setConstant();
        RooRealVar GPion_Ks_P("GPion_Ks_P","GPion_Ks_P", 47.4,  10,  200);
        GPion_Ks_P.setConstant();
        RooRealVar rhoPion_Ks_P("rhoPion_Ks_P","rhoPion_Ks_P", 0.20, 0.0, 30.0);
        RooRealVar phiPion_Ks_P("phiPion_Ks_P","phiPion_Ks_P", -2.22818 ,-pi,pi);

        fitparameter.add(rhoPion_Ks_P);
        fitparameter.add(phiPion_Ks_P);
        ksksKpdf.addResonance("DtoPiPion_Ks_P", "D+ to (pi+ Ks)p eta ", spin1, mPion_Ks_P,  GPion_Ks_P,  
                rhoPion_Ks_P, phiPion_Ks_P, modetype12, r0, r1, g0);
    }
    if(AddResonance["Pion_Ks_D"]){
        // ( pi+ Ks ) P wave
        RooRealVar mPion_Ks_D("mPion_Ks_D",  "mPion_Ks_D", 1.60 , 1.55,   1.65);
        mPion_Ks_D.setConstant();
        RooRealVar GPion_Ks_D("GPion_Ks_D","GPion_Ks_D", 47.4,  10,  200);
        GPion_Ks_D.setConstant();
        RooRealVar rhoPion_Ks_D("rhoPion_Ks_D","rhoPion_Ks_D", 0.20, 0.0, 30.0);
        RooRealVar phiPion_Ks_D("phiPion_Ks_D","phiPion_Ks_D", -2.22818 ,-pi,pi);

        fitparameter.add(rhoPion_Ks_D);
        fitparameter.add(phiPion_Ks_D);
        ksksKpdf.addResonance("DtoPiPion_Ks_D", "D+ to (pi+ Ks)p eta ", spin2, mPion_Ks_D,  GPion_Ks_D,  
                rhoPion_Ks_D, phiPion_Ks_D, modetype12, r0, r1, g0);
    }
    if(AddResonance["Ks_Eta_P"]){
        RooRealVar mKs_Eta_P("mKs_Eta_P",  "mKs_Eta_P", 1 , .20,   34);
        mKs_Eta_P.setConstant();
        RooRealVar GKs_Eta_P("GKs_Eta_P","GKs_Eta_P", 44,  10,  2000);
        GKs_Eta_P.setConstant();
        RooRealVar rhoKs_Eta_P("rhoKs_Eta_P","rhoKs_Eta_P", 0.010, 0.0, 100.0);
        RooRealVar phiKs_Eta_P("phiKs_Eta_P","phiKs_Eta_P", -2.22818 ,-pi,pi);

        fitparameter.add(rhoKs_Eta_P);
        fitparameter.add(phiKs_Eta_P);
        ksksKpdf.addResonance("DtoPiKs_Eta_P", "D+ to (eta Ks)p pi+ ", spin1, mKs_Eta_P,  GKs_Eta_P, 
                rhoKs_Eta_P, phiKs_Eta_P, modetype13, r0, r1, g0); 
    }
    if(AddResonance["Ks_Eta_D"]){
        RooRealVar mKs_Eta_D("mKs_Eta_D",  "mKs_Eta_D", 1 , .20,   34);
        mKs_Eta_D.setConstant();
        RooRealVar GKs_Eta_D("GKs_Eta_D","GKs_Eta_D", 44,  10,  2000);
        GKs_Eta_D.setConstant();
        RooRealVar rhoKs_Eta_D("rhoKs_Eta_D","rhoKs_Eta_D", 0.010, 0.0, 100.0);
        RooRealVar phiKs_Eta_D("phiKs_Eta_D","phiKs_Eta_D", -2.22818 ,-pi,pi);

        fitparameter.add(rhoKs_Eta_D);
        fitparameter.add(phiKs_Eta_D);
        ksksKpdf.addResonance("DtoPiKs_Eta_D", "D+ to (eta Ks)p pi+ ", spin2, mKs_Eta_D,  GKs_Eta_D, 
                rhoKs_Eta_D, phiKs_Eta_D, modetype13, r0, r1, g0); 
    }
    if(AddResonance["Eta_Pion_P"]){
        RooRealVar mEta_Pion_P("mEta_Pion_P",  "mEta_Pion_P", 1 , .20,   34);
        mEta_Pion_P.setConstant();
        RooRealVar GEta_Pion_P("GEta_Pion_P","GEta_Pion_P", 44,  10,  2000);
        GEta_Pion_P.setConstant();
        RooRealVar rhoEta_Pion_P("rhoEta_Pion_P","rhoEta_Pion_P", 1.0, 0.0, 10.0);
        RooRealVar phiEta_Pion_P("phiEta_Pion_P","phiEta_Pion_P", 1.18 ,-pi,pi);

        fitparameter.add(rhoEta_Pion_P);
        fitparameter.add(phiEta_Pion_P);
        ksksKpdf.addResonance("DtoPiEta_Pion_P", "D+ to (pi+ eta)p Ks ", spin1, mEta_Pion_P,  GEta_Pion_P,  
                rhoEta_Pion_P, phiEta_Pion_P, modetype23, r0, r1, g0); 
    }
    if(AddResonance["Eta_Pion_D"]){
        RooRealVar mEta_Pion_D("mEta_Pion_D",  "mEta_Pion_D", 1 , .20,   34);
        mEta_Pion_D.setConstant();
        RooRealVar GEta_Pion_D("GEta_Pion_D","GEta_Pion_D", 44,  10,  2000);
        GEta_Pion_D.setConstant();
        RooRealVar rhoEta_Pion_D("rhoEta_Pion_D","rhoEta_Pion_D", 1.0, 0.0, 10.0);
        RooRealVar phiEta_Pion_D("phiEta_Pion_D","phiEta_Pion_D", 1.18 ,-pi,pi);

        fitparameter.add(rhoEta_Pion_D);
        fitparameter.add(phiEta_Pion_D);
        ksksKpdf.addResonance("DtoPiEta_Pion_D", "D+ to (pi+ eta)p Ks ", spin2, mEta_Pion_D,  GEta_Pion_D,  
                rhoEta_Pion_D, phiEta_Pion_D, modetype23, r0, r1, g0); 
    }
    */

    //ksksKpdf.Print();
    cout<<"========  test ========"<<endl;
    ksksKpdf.test();
    cout<<"========= end of test ======"<<endl;
    exit();
    RooRealVar weight("weight","weight", -2, 2);

    RooArgSet theSet1,theSet2,theSet3;
    theSet1.add(RooArgSet(v11, v12, v13, v14, v21, v22, v23, v24));
    theSet2.add(RooArgSet(v31, v32, v33, v34));
    theSet3.add(RooArgSet(weight));
    RooArgSet theSet4(theSet1,theSet2,"");
    RooArgSet theSet(theSet4,theSet3,"");

    RooDataSet *data = RooDataSet::read(Datadat,theSet);
    data->Print();

    // RooDataSet *datEtaV2G = new
    //     RooDataSet(data->GetName(),data->GetTitle(),
    //             data,*data->get(),0, weight.GetName());
    // datEtaV2G->Print();

    //RooFitResult *result = ksksKpdf.fitTo(*datEtaV2G , Save(kTRUE));
    //return ;
    RooFitResult *result = ksksKpdf.fitTo(*data , Save(kTRUE));
    TString name = "fitresult";
    result->Print();

  //  ksksKpdf.project("project_1.root");
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

void AddRes(TString resonanceList , map<string, int> &AddResonance){
    AddResonance["a0_980"] = 1;  // ***
    AddResonance["NoneRes"] = 0; // **
    AddResonance["a2_1320"] = 0;  // 
    AddResonance["a0_1450"] = 0;  // 
    AddResonance["pi1_1400"] = 0;  // 
    AddResonance["rho_1700"] = 0;  // 
    AddResonance["Pion_Ks_P"] = 0; // 
    AddResonance["Pion_Ks_D"] = 0; // 
    AddResonance["Ks_Eta_P"] = 0;  // 
    AddResonance["Ks_Eta_D"] = 0;  // 
    AddResonance["Eta_Pion_P"] = 0; // *
    AddResonance["Eta_Pion_D"] = 0; // *
    if (resonanceList != "") {
        TMVA::Tools::Instance();
        vector<TString> mlist = TMVA::gTools().SplitString( resonanceList, ',' );

        for (int i=0; i<mlist.size(); i++) {
            string res(mlist[i]);

            if (AddResonance.find(res) == AddResonance.end()) {
                cout << "The resonance \"" << res 
                    << "\" is not defined yed, choose from those" << endl;
                for(map<string,int>::iterator it = AddResonance.begin(); it != AddResonance.end(); it++) {
                    std::cout << it->first << " ";
                }
                std::cout << std::endl;
            }
            else{
                cout<<"Add "<<res<<" successful!!!"<<endl;}
            AddResonance[res] = 1;
        }
    }
}
