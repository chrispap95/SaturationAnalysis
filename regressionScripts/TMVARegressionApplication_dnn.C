/// \file
/// \ingroup tutorial_tmva
/// \notebook -nodraw
/// This macro provides a simple example on how to use the trained regression MVAs
/// within an analysis module
///
///  - Project   : TMVA - a Root-integrated toolkit for multivariate data analysis
///  - Package   : TMVA
///  - Exectuable: TMVARegressionApplication
///
/// \macro_output
/// \macro_code
/// \author Andreas Hoecker

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

using namespace TMVA;

void TMVARegressionApplication_dnn(int energy) {
    //---------------------------------------------------------------
    // This loads the library
    TMVA::Tools::Instance();

    std::cout << std::endl;
    std::cout << "==> Start TMVARegressionApplication" << std::endl;

    // --------------------------------------------------------------------------------------------------
    // --- Create the Reader object

    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

    // Create a set of variables and declare them to the reader
    // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
    Float_t layer, eta, phi, n1, n2, n3, n4, n5, n6, nup, ndown, rechitsum, event;
    Float_t un1, un2, un3, un4, un5, un6;
    Float_t dn1, dn2, dn3, dn4, dn5, dn6;
    Float_t cellType, rechit;
    reader->AddVariable(    "layer",    &layer);
    reader->AddVariable(       "n1",       &n1);
    reader->AddVariable(       "n2",       &n2);
    reader->AddVariable(       "n3",       &n3);
    reader->AddVariable(       "n4",       &n4);
    reader->AddVariable(       "n5",       &n5);
    reader->AddVariable(       "n6",       &n6);
    reader->AddVariable(      "nup",      &nup);
    reader->AddVariable(    "ndown",    &ndown);
    reader->AddVariable(      "un1",      &un1);
    reader->AddVariable(      "un2",      &un2);
    reader->AddVariable(      "un3",      &un3);
    reader->AddVariable(      "un4",      &un4);
    reader->AddVariable(      "un5",      &un5);
    reader->AddVariable(      "un6",      &un6);
    reader->AddVariable(      "dn1",      &dn1);
    reader->AddVariable(      "dn2",      &dn2);
    reader->AddVariable(      "dn3",      &dn3);
    reader->AddVariable(      "dn4",      &dn4);
    reader->AddVariable(      "dn5",      &dn5);
    reader->AddVariable(      "dn6",      &dn6);
    //reader->AddVariable( "cellType", &cellType);

    // Spectator variables declared in the training have to be added to the reader, too
    //Float_t spec1,spec2;
    //reader->AddSpectator( "spec1:=var1*2",  &spec1 );
    //reader->AddSpectator( "spec2:=var1*3",  &spec2 );

    // --- Book the MVA methods

    TString dir    = "dataset/weights/";
    TString prefix = "TMVAReg_saturated_E500to3000Eta1p56Phi0p0_240K_3hl_25nodes";

    // Book method(s)
    TString methodName = "DNN_CPU method";//it->first + " method";
    TString weightfile = dir + prefix + "_" + "DNN_CPU" + ".weights.xml";
    reader->BookMVA( methodName, weightfile );

    // Prepare input tree (this must be replaced by your data source)
    // in this example, there is a toy tree with signal and one with background events
    // we'll later on use only the "signal" events for the test in this example.
    //
    string energy_str = to_string(energy);
    TFile *input(0);
    TString fname = "EvaluationSamples2/out_E"+energy_str+"Eta1p56Phi0p0_0_converted.root";
    if (!gSystem->AccessPathName( fname )) {
        input = TFile::Open( fname ); // check if file in local directory exists
    }
    //else {
    //   TFile::SetCacheFileDir(".");
    //   input = TFile::Open("http://root.cern.ch/files/tmva_reg_example.root", "CACHEREAD"); // if not: download from ROOT server
    // }
    if (!input) {
        std::cout << "ERROR: could not open data file" << std::endl;
        exit(1);
    }
    std::cout << "--- TMVARegressionApp        : Using input file: " << input->GetName() << std::endl;

    // --- Event loop

    // Prepare the tree
    // - here the variable names have to corresponds to your tree
    // - you can use the same variables as above which is slightly faster,
    //   but of course you can use different ones and copy the values inside the event loop
    //
    TTree* theTree = (TTree*)input->Get("t1");
    std::cout << "--- Select signal sample" << std::endl;
    theTree->SetBranchAddress(    "layer",      &layer );
    theTree->SetBranchAddress(       "n1",         &n1 );
    theTree->SetBranchAddress(       "n2",         &n2 );
    theTree->SetBranchAddress(       "n3",         &n3 );
    theTree->SetBranchAddress(       "n4",         &n4 );
    theTree->SetBranchAddress(       "n5",         &n5 );
    theTree->SetBranchAddress(       "n6",         &n6 );
    theTree->SetBranchAddress(       "un1",       &un1 );
    theTree->SetBranchAddress(       "un2",       &un2 );
    theTree->SetBranchAddress(       "un3",       &un3 );
    theTree->SetBranchAddress(       "un4",       &un4 );
    theTree->SetBranchAddress(       "un5",       &un5 );
    theTree->SetBranchAddress(       "un6",       &un6 );
    theTree->SetBranchAddress(       "dn1",       &dn1 );
    theTree->SetBranchAddress(       "dn2",       &dn2 );
    theTree->SetBranchAddress(       "dn3",       &dn3 );
    theTree->SetBranchAddress(       "dn4",       &dn4 );
    theTree->SetBranchAddress(       "dn5",       &dn5 );
    theTree->SetBranchAddress(       "dn6",       &dn6 );
    theTree->SetBranchAddress(       "nup",       &nup );
    theTree->SetBranchAddress(     "ndown",     &ndown );
    theTree->SetBranchAddress(  "cellType",  &cellType ); // 0 for 300um & 1 for 200um
    theTree->SetBranchAddress(     "event",     &event );
    theTree->SetBranchAddress( "rechitsum", &rechitsum );
    theTree->SetBranchAddress(    "rechit",    &rechit ); // Unsaturated rechit estimation

    TString foutname = "RegressionResults2/flatRegressionResult_"+energy_str+"GeV_Eta1p56Phi0p0.root";
    TFile *target  = new TFile( foutname,"RECREATE" );
    TTree* t1 = new TTree("t1","sample");
    Float_t val;
    t1->Branch(      "layer",     &layer,      "layer/F" );
    t1->Branch(         "n1",        &n1,         "n1/F" );
    t1->Branch(         "n2",        &n2,         "n2/F" );
    t1->Branch(         "n3",        &n3,         "n3/F" );
    t1->Branch(         "n4",        &n4,         "n4/F" );
    t1->Branch(         "n5",        &n5,         "n5/F" );
    t1->Branch(         "n6",        &n6,         "n6/F" );
    t1->Branch(        "un1",       &un1,        "un1/F" );
    t1->Branch(        "un2",       &un2,        "un2/F" );
    t1->Branch(        "un3",       &un3,        "un3/F" );
    t1->Branch(        "un4",       &un4,        "un4/F" );
    t1->Branch(        "un5",       &un5,        "un5/F" );
    t1->Branch(        "un6",       &un6,        "un6/F" );
    t1->Branch(        "dn1",       &dn1,        "dn1/F" );
    t1->Branch(        "dn2",       &dn2,        "dn2/F" );
    t1->Branch(        "dn3",       &dn3,        "dn3/F" );
    t1->Branch(        "dn4",       &dn4,        "dn4/F" );
    t1->Branch(        "dn5",       &dn5,        "dn5/F" );
    t1->Branch(        "dn6",       &dn6,        "dn6/F" );
    t1->Branch(        "nup",       &nup,        "nup/F" );
    t1->Branch(      "ndown",     &ndown,      "ndown/F" );
    t1->Branch(   "cellType",  &cellType,   "cellType/F" );
    t1->Branch(      "event",     &event,      "event/F" );
    t1->Branch(  "rechitsum", &rechitsum,  "rechitsum/F" );
    t1->Branch(     "rechit",    &rechit,     "rechit/F" );
    t1->Branch( "regression",       &val, "regression/F" );

    std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
    TStopwatch sw;
    sw.Start();

    float avdevquad = 0;
    int n = 0;
    for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
        if (ievt%50000 == 0) {
            std::cout << "--- ... Processing event: " << ievt << std::endl;
        }

        theTree->GetEntry(ievt);

        // Retrieve the MVA target values (regression outputs) and fill into histograms
        // NOTE: EvaluateRegression(..) returns a vector for multi-target regression
        if (layer > -1) {
            val = (reader->EvaluateRegression( methodName ))[0];
            avdevquad += pow(rechit-val,2);
            n++;
        }else {
            val = 0;
        }
        t1->Fill();
    }
    sw.Stop();
    std::cout << "--- End of event loop: "; sw.Print();

    std::cout << "average quadratic deviation = " << sqrt(avdevquad/(float)n) << std::endl;

    target->cd();
    t1->Write();
    target->Close();

    std::cout << "--- Created root file: \"" << target->GetName()
    << "\" containing the MVA output tree" << std::endl;

    delete reader;

    std::cout << "==> TMVARegressionApplication is done!" << std::endl << std::endl;
}
