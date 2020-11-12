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
//#include "TMVA/DataLoader.h"
//#include "TMVA/TMVARegGui.h"

using namespace TMVA;

void TMVARegression()
{
  TMVA::Tools::Instance();

  TString uniqueid = "TMVAReg_saturated_E500to3000Eta1p56Phi0p0_CV_240K_3hl_25nodes";

  std::cout << std::endl;
  std::cout << "==> Start TMVARegression" << std::endl;

  TString outfileName( uniqueid+".root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");

  dataloader->AddVariable(    "layer",         "layer", "units", 'F' );
  dataloader->AddVariable(       "n1",    "neighbor 1", "units", 'F' );
  dataloader->AddVariable(       "n2",    "neighbor 2", "units", 'F' );
  dataloader->AddVariable(       "n3",    "neighbor 3", "units", 'F' );
  dataloader->AddVariable(       "n4",    "neighbor 4", "units", 'F' );
  dataloader->AddVariable(       "n5",    "neighbor 5", "units", 'F' );
  dataloader->AddVariable(       "n6",    "neighbor 6", "units", 'F' );
  dataloader->AddVariable(      "nup",   "neighbor up", "units", 'F' );
  dataloader->AddVariable(    "ndown", "neighbor down", "units", 'F' );
  dataloader->AddVariable(      "un1",          "up 1", "units", 'F' );
  dataloader->AddVariable(      "un2",          "up 2", "units", 'F' );
  dataloader->AddVariable(      "un3",          "up 3", "units", 'F' );
  dataloader->AddVariable(      "un4",          "up 4", "units", 'F' );
  dataloader->AddVariable(      "un5",          "up 5", "units", 'F' );
  dataloader->AddVariable(      "un6",          "up 6", "units", 'F' );
  dataloader->AddVariable(      "dn1",        "down 1", "units", 'F' );
  dataloader->AddVariable(      "dn2",        "down 2", "units", 'F' );
  dataloader->AddVariable(      "dn3",        "down 3", "units", 'F' );
  dataloader->AddVariable(      "dn4",        "down 4", "units", 'F' );
  dataloader->AddVariable(      "dn5",        "down 5", "units", 'F' );
  dataloader->AddVariable(      "dn6",        "down 6", "units", 'F' );
  //dataloader->AddVariable( "cellType",      "cellType", "units", 'F' );

  dataloader->AddTarget( "rechit" );

  //load the signal and background event samples from ROOT trees
  TFile *input(0);
  TString fname = "TrainingSamples.nosync/out_E500to3000Eta1p56Phi0_0_converted.root";
  if (!gSystem->AccessPathName( fname )) {
    input = TFile::Open( fname ); // check if file in local directory exists
  }
  if (!input) {
    std::cout << "ERROR: could not open data file" << std::endl;
    exit(1);
  }
  std::cout << "--- TMVACrossValidationRegression: Using input file: " << input->GetName() << std::endl;

  // Register the regression tree
  TTree *regTree = (TTree*)input->Get("t1");

  Double_t regWeight  = 1.0;
  dataloader->AddRegressionTree( regTree, regWeight );
  TCut mycut = "";
  dataloader->PrepareTrainingAndTestTree(mycut,
    "nTest_Regression=0:"
    "SplitMode=Random:NormMode=NumEvents:!V");

  UInt_t numFolds = 5;
  TString analysisType = "Regression";
  TString splitExpr = "";

  TString cvOptions = Form("!V:!Silent:ModelPersistence"
    ":!FoldFileOutput:AnalysisType=%s:NumFolds=%i:SplitExpr=%s",
    analysisType.Data(), numFolds, splitExpr.Data());

  TMVA::CrossValidation ce{"TMVACrossValidationRegression", dataloader, outputFile, cvOptions};


  TString layoutString("Layout=SYMMRELU|7,Layout=SYMMRELU|10,Layout=SYMMRELU|10,Layout=SYMMRELU|10,"
                       "LINEAR");
  TString training0("LearningRate=5e-4,Momentum=0.5,Repetitions=1,ConvergenceSteps=20,BatchSize=200,"
                    "TestRepetitions=10,WeightDecay=0.01,Regularization=NONE,DropConfig=0.2+0.2+0.2+0.,"
                    "DropRepetitions=2");
  TString training1("LearningRate=5e-4,Momentum=0.7,Repetitions=1,ConvergenceSteps=20,BatchSize=200,"
                    "TestRepetitions=5,WeightDecay=0.01,Regularization=L2,DropConfig=0.1+0.1+0.1,"
                    "DropRepetitions=1");
  TString training2("LearningRate=1e-4,Momentum=0.3,Repetitions=1,ConvergenceSteps=20,BatchSize=200,"
                    "TestRepetitions=5,WeightDecay=0.01,Regularization=NONE");
  TString training3("LearningRate=1e-4,Momentum=0.1,Repetitions=1,ConvergenceSteps=20,BatchSize=200,"
                    "TestRepetitions=5,WeightDecay=0.01,Regularization=NONE");

  TString trainingStrategyString("TrainingStrategy=");
  trainingStrategyString += training0 + "|" + training1 + "|" + training2 + "|" + training3;
  TString nnOptions(
    "!H:V:ErrorStrategy=SUMOFSQUARES:VarTransform=G:WeightInitialization=XAVIERUNIFORM:Architecture=CPU");
  nnOptions.Append(":");
  nnOptions.Append(layoutString);
  nnOptions.Append(":");
  nnOptions.Append(trainingStrategyString);

  ce.BookMethod(TMVA::Types::kDNN, "DNN_CPU", nnOptions);
  ce.Evaluate();

  // Save the output
  outputFile->Close();

  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVARegression is done!" << std::endl;

  delete dataloader;

  // Launch the GUI for the root macros
  if (!gROOT->IsBatch()) TMVA::TMVARegGui( outfileName );
}

int main( int argc, char** argv )
{
   TMVARegression();
   return 0;
}
