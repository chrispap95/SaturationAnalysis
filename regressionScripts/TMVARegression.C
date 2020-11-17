#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Tools.h"
#include "TMVA/Factory.h"

using namespace TMVA;

void TMVARegression( TString fname, TString uniqueID = "testRun", int Nsamples = 1000,
                     int Nlayers = 3, TString nodes = "20")
{
    ROOT::EnableImplicitMT();
    // This loads the library
    TMVA::Tools::Instance();

    TString uniqueid = uniqueID;

    std::cout << std::endl;
    std::cout << "==> Start TMVARegression" << std::endl;

    // Create a new root output file
    TString outfileName( uniqueid+".root" );
    TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

    /*
    Create the factory object. Later you can choose the methods
    whose performance you'd like to investigate. The factory will
    then run the performance analysis for you.
    The first argument is the base of the name of all the
    weightfiles in the directory weight/
    The second argument is the output file for the training results
    All TMVA output can be suppressed by removing the "!" (not) in
    front of the "Silent" argument in the option string
    */

    TMVA::Factory *factory = new TMVA::Factory( uniqueid, outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression" );

    TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");
    /*
    If you wish to modify default settings
    (please check "src/Config.h" to see all available global options)
    */

    //     (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
    //     (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

    /*
    Define the input variables that shall be used for the MVA training
    note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
    [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
    */
    //dataloader->AddVariable( "layer", "layer", "units", 'F' );
    dataloader->AddVariable( "n1", "neighbor 1", "units", 'F' );
    dataloader->AddVariable( "n2", "neighbor 2", "units", 'F' );
    dataloader->AddVariable( "n3", "neighbor 3", "units", 'F' );
    dataloader->AddVariable( "n4", "neighbor 4", "units", 'F' );
    dataloader->AddVariable( "n5", "neighbor 5", "units", 'F' );
    dataloader->AddVariable( "n6", "neighbor 6", "units", 'F' );
    /*
    dataloader->AddVariable( "nup", "neighbor up", "units", 'F' );
    dataloader->AddVariable( "ndown", "neighbor down", "units", 'F' );
    dataloader->AddVariable( "un1", "up 1", "units", 'F' );
    dataloader->AddVariable( "un2", "up 2", "units", 'F' );
    dataloader->AddVariable( "un3", "up 3", "units", 'F' );
    dataloader->AddVariable( "un4", "up 4", "units", 'F' );
    dataloader->AddVariable( "un5", "up 5", "units", 'F' );
    dataloader->AddVariable( "un6", "up 6", "units", 'F' );
    dataloader->AddVariable( "dn1", "down 1", "units", 'F' );
    dataloader->AddVariable( "dn2", "down 2", "units", 'F' );
    dataloader->AddVariable( "dn3", "down 3", "units", 'F' );
    dataloader->AddVariable( "dn4", "down 4", "units", 'F' );
    dataloader->AddVariable( "dn5", "down 5", "units", 'F' );
    dataloader->AddVariable( "dn6", "down 6", "units", 'F' );
    //dataloader->AddVariable( "cellType", "cellType", "units", 'F' );
    */

    /*
    You can add so-called "Spectator variables", which are not used in the MVA training,
    but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
    input variables, the response values of all trained MVAs, and the spectator variables
    */
    //dataloader->AddSpectator( "spec1:=rechitsum",  "Spectator 1", "units", 'F' );
    //dataloader->AddSpectator( "spec2:=var1*3",  "Spectator 2", "units", 'F' );


    // Add the variable carrying the regression target
    dataloader->AddTarget( "rechit" );

    /*It is also possible to declare additional targets for multi-dimensional regression, ie:
          factory->AddTarget( "fvalue2" );
    BUT: this is currently ONLY implemented for MLP
    Read training and test data (see TMVAClassification for reading ASCII files)
    */

    //load the signal and background event samples from ROOT trees
    TFile *input(0);
    if (!gSystem->AccessPathName( fname )) {
        input = TFile::Open( fname ); // check if file in local directory exists
    }
    if (!input) {
        std::cout << "ERROR: could not open data file" << std::endl;
        exit(1);
    }
    std::cout << "--- TMVARegression           : Using input file: " << input->GetName() << std::endl;

    // Register the regression tree

    TTree *regTree = (TTree*)input->Get("tree");

    // global event weights per tree (see below for setting event-wise weights)
    Double_t regWeight  = 1.0;

    // You can add an arbitrary number of regression trees
    dataloader->AddRegressionTree( regTree, regWeight );
    //dataloader->SetWeightExpression("MLdead","Regression");

    /*
    This would set individual event weights (the variables defined in the
    expression need to exist in the original TTree)
    */
    //dataloader->SetWeightExpression( "var1", "Regression" );

    // Apply additional cuts on the signal and background samples (can be different)
    TCut mycut = "";
    // for example: TCut mycut = "abs(var1)<0.5 && abs(var2-0.5)<1";

    // split to training and testing samples
    int nTrain = 0.8*Nsamples;
    int nTest = 0.2*Nsamples;

    // tell the DataLoader to use all remaining events in the trees after training for testing:
    dataloader->PrepareTrainingAndTestTree(mycut,"nTrain_Regression="+to_string(nTrain)+
                                                 ":nTest_Regression="+to_string(nTest)+
                                                 ":SplitMode=Random:NormMode=NumEvents:!V");

    /*
    If no numbers of events are given, half of th e events in the tree are used
    for training, and the other half for testing:
          dataloader->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
    Book MVA methods
    Please lookup the various method configuration options in the corresponding cxx files, eg:
    src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.
    it is possible to preset ranges in the option string in which the cut optimisation should be done:
    "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input
    */

    TString layoutString("Layout=TANH|22");
    for (int i = 0; i < Nlayers; ++i) {
        layoutString += ",TANH|"+nodes;
    }
    layoutString += ",LINEAR";

    TString trainingStrategyString("TrainingStrategy=LearningRate=1e-3,Momentum=0.3,"
                                   "ConvergenceSteps=20,BatchSize=50,TestRepetitions=1,"
                                   "WeightDecay=0.0,Regularization=None,Optimizer=Adam");

    TString nnOptions("!H:V:ErrorStrategy=SUMOFSQUARES:VarTransform=G:WeightInitialization=XAVIERUNIFORM:Architecture=CPU");
    nnOptions.Append(":");
    nnOptions.Append(layoutString);
    nnOptions.Append(":");
    nnOptions.Append(trainingStrategyString);

    factory->BookMethod(dataloader, TMVA::Types::kDL, "DNN_CPU", nnOptions);

    // Now you can tell the factory to train, test, and evaluate the MVAs
    // Train MVAs using the set of training events
    factory->TrainAllMethods();

    // Evaluate all MVAs using the set of test events
    factory->TestAllMethods();

    // Evaluate and compare performance of all configured MVAs
    factory->EvaluateAllMethods();

    // Save the output
    outputFile->Close();

    std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
    std::cout << "==> TMVARegression is done!" << std::endl;

    delete factory;
    delete dataloader;

    // Launch the GUI for the root macros
    if (!gROOT->IsBatch()) TMVA::TMVARegGui( outfileName );
    return;
}
