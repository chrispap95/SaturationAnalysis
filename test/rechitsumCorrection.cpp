/****************************************************
**  Channelog:
**      Moved to CMSSW output.
****************************************************/

#include<string>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<map>
#include <boost/algorithm/string.hpp>
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "boost/function.hpp"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TMath.h"
#ifdef _DEBUG
#include "debug_new.h"
#endif
#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"

using boost::lexical_cast;
namespace po=boost::program_options;

double DeltaR(double eta1,double phi1,double eta2,double phi2){
    double dr=99999.;
    double deta=fabs(eta1-eta2);
    double dphi=fabs(phi1-phi2);
    if(dphi>TMath::Pi()) dphi=2.*TMath::Pi()-dphi;
    dr=sqrt(deta*deta+dphi*dphi);
    return dr;
}

int main(int argc, char** argv){
    /**********************************
    ** initialize some variables
    **********************************/

    //Input output and config options
    std::string cfg;
    unsigned pNevts;
    std::string outFilePath;
    std::string filePath;
    std::string digifilePath;
    unsigned nRuns;
    std::string recoFileName;
    unsigned debug;
    double deadfrac;
    bool adjacent;
    po::options_description preconfig("Configuration");
    preconfig.add_options()("cfg,c",po::value<std::string>(&cfg)->required());
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(preconfig).allow_unregistered().run(), vm);
    po::notify(vm);
    po::options_description config("Configuration");
    config.add_options()
    //Input output and config options //->required()
    ("pNevts,n",        po::value<unsigned>(&pNevts)->default_value(0))
    ("outFilePath,o",   po::value<std::string>(&outFilePath)->required())
    ("filePath,i",      po::value<std::string>(&filePath)->required())
    ("recoFileName,r",  po::value<std::string>(&recoFileName)->required())
    ("nRuns",           po::value<unsigned>(&nRuns)->default_value(0))
    ("debug,d",         po::value<unsigned>(&debug)->default_value(0))
    ("deadfrac",        po::value<double>(&deadfrac)->default_value(0))
    //Restrict number of adjacent dead cells
    ("adjacent",        po::value<bool>(&adjacent)->default_value(0))
    ;
    po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
    po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
    po::notify(vm);

    std::cout << " -- Input parameters: " << std::endl
    << " -- Input file path: " << filePath << std::endl
    << " -- Output file path: " << outFilePath << std::endl
    << std::endl
    << " -- Processing ";
    if (pNevts == 0) std::cout << "all events." << std::endl;
    else std::cout << pNevts << " events per run." << std::endl;

    TRandom3 lRndm(0);
    std::cout << " -- Random number seed: " << lRndm.GetSeed() << std::endl;

    /**********************************
    ** Input
    **********************************/
    // Get Reco files
    std::ostringstream inputrec;
    if (digifilePath.size()==0)
    inputrec << filePath << "/" << recoFileName;
    else
    inputrec << digifilePath << "/" << recoFileName;

    TChain *lRecTree = 0;

    lRecTree = new TChain("hgcalTupleTree/tree");

    if (nRuns == 0){
        lRecTree->AddFile(inputrec.str().c_str());
    }
    else {
        for (unsigned i(1);i<=nRuns;++i){
            std::ostringstream lstrrec;
            lstrrec << inputrec.str() << "_" << i << ".root";
            lRecTree->AddFile(lstrrec.str().c_str());
        }
    }

    if (!lRecTree){
        std::cout << " -- Error, tree RecoTree cannot be opened. Exiting..." << std::endl;
        return 1;
    }

    /**********************************
    ** Output
    **********************************/
    // Define output file and the histograms contained
    TFile *outputFile = TFile::Open(outFilePath.c_str(),"RECREATE");

    if (!outputFile) {
        std::cout << " -- Error, output file " << outFilePath
        << " cannot be opened. Please create output directory. Exiting..." << std::endl;
        return 1;
    }
    else {
        std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
    }
    outputFile->cd();

    TH1F* h1 = new TH1F("h1","uncorrected",250,600,1100);
    TH1F* h2 = new TH1F("h2","corrected",250,600,1100);

    /**********************************
    **  start event loop
    **********************************/
    const unsigned nEvts = (
        (pNevts > lRecTree->GetEntries() || pNevts==0)
        ? static_cast<unsigned>(lRecTree->GetEntries())
        : pNevts
    );
    std::cout << " -- Processing " << nEvts << " events out of "
    << lRecTree->GetEntries() << std::endl;

    //loop on events
    std::vector<float   > *rechitEnergy = 0;
    std::vector<float   > *rechitEta    = 0;
    std::vector<float   > *rechitPhi    = 0;
    std::vector<float   > *rechitPosx   = 0;
    std::vector<float   > *rechitPosy   = 0;
    std::vector<float   > *rechitPosz   = 0;
    std::vector<unsigned> *rechitLayer  = 0;
    std::vector<unsigned> *rechitIndex  = 0;
    std::vector<int     > *rechitWaferU = 0;
    std::vector<int     > *rechitWaferV = 0;
    std::vector<unsigned> *rechitCellU  = 0;
    std::vector<unsigned> *rechitCellV  = 0;
    std::vector<float   > *simhitEnergy = 0;
    std::vector<float   > *simhitEta    = 0;
    std::vector<float   > *simhitPhi    = 0;
    std::vector<float   > *simhitPosx   = 0;
    std::vector<float   > *simhitPosy   = 0;
    std::vector<float   > *simhitPosz   = 0;
    std::vector<unsigned> *simhitLayer  = 0;
    std::vector<unsigned> *simhitIndex  = 0;
    std::vector<int     > *simhitWaferU = 0;
    std::vector<int     > *simhitWaferV = 0;
    std::vector<unsigned> *simhitCellU  = 0;
    std::vector<unsigned> *simhitCellV  = 0;
    std::vector<float   > *genEta       = 0;
    std::vector<float   > *genPhi       = 0;

    lRecTree->SetBranchAddress("HGCRecHitEnergy" ,&rechitEnergy);
    lRecTree->SetBranchAddress("HGCRecHitEta"    ,&rechitEta);
    lRecTree->SetBranchAddress("HGCRecHitPhi"    ,&rechitPhi);
    lRecTree->SetBranchAddress("HGCRecHitPosx"   ,&rechitPosx);
    lRecTree->SetBranchAddress("HGCRecHitPosy"   ,&rechitPosy);
    lRecTree->SetBranchAddress("HGCRecHitPosz"   ,&rechitPosz);
    lRecTree->SetBranchAddress("HGCRecHitLayer"  ,&rechitLayer);
    lRecTree->SetBranchAddress("HGCRecHitIndex"  ,&rechitIndex);
    lRecTree->SetBranchAddress("HGCRecHitWaferU" ,&rechitWaferU);
    lRecTree->SetBranchAddress("HGCRecHitWaferV" ,&rechitWaferV);
    lRecTree->SetBranchAddress("HGCRecHitCellU"  ,&rechitCellU);
    lRecTree->SetBranchAddress("HGCRecHitCellV"  ,&rechitCellV);
    lRecTree->SetBranchAddress("HGCSimHitsEnergy",&simhitEnergy);
    lRecTree->SetBranchAddress("HGCSimHitsEta"   ,&simhitEta);
    lRecTree->SetBranchAddress("HGCSimHitsPhi"   ,&simhitPhi);
    lRecTree->SetBranchAddress("HGCSimHitsPosx"  ,&simhitPosx);
    lRecTree->SetBranchAddress("HGCSimHitsPosy"  ,&simhitPosy);
    lRecTree->SetBranchAddress("HGCSimHitsPosz"  ,&simhitPosz);
    lRecTree->SetBranchAddress("HGCSimHitsLayer" ,&simhitLayer);
    lRecTree->SetBranchAddress("HGCSimHitsIndex" ,&simhitIndex);
    lRecTree->SetBranchAddress("HGCSimHitsWaferU",&simhitWaferU);
    lRecTree->SetBranchAddress("HGCSimHitsWaferV",&simhitWaferV);
    lRecTree->SetBranchAddress("HGCSimHitsCellU" ,&simhitCellU);
    lRecTree->SetBranchAddress("HGCSimHitsCellV" ,&simhitCellV);
    lRecTree->SetBranchAddress("GenParEta"       ,&genEta);
    lRecTree->SetBranchAddress("GenParPhi"       ,&genPhi);

    unsigned ievtRec = 0;

    // Loop over entries (events)
    for (unsigned ievt(0); ievt<nEvts; ++ievt){
        if (ievtRec>=lRecTree->GetEntries()) continue;
        Long64_t local_entry = lRecTree->LoadTree(ievt);

        if (debug) std::cout << std::endl<<std::endl<<"... Processing entry: " << ievt << std::endl;
        else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;

        if (local_entry < 0) continue;
        if (local_entry == 0) {
            lRecTree->SetBranchAddress("HGCRecHitEnergy" ,&rechitEnergy);
            lRecTree->SetBranchAddress("HGCRecHitEta"    ,&rechitEta);
            lRecTree->SetBranchAddress("HGCRecHitPhi"    ,&rechitPhi);
            lRecTree->SetBranchAddress("HGCRecHitPosx"   ,&rechitPosx);
            lRecTree->SetBranchAddress("HGCRecHitPosy"   ,&rechitPosy);
            lRecTree->SetBranchAddress("HGCRecHitPosz"   ,&rechitPosz);
            lRecTree->SetBranchAddress("HGCRecHitLayer"  ,&rechitLayer);
            lRecTree->SetBranchAddress("HGCRecHitIndex"  ,&rechitIndex);
            lRecTree->SetBranchAddress("HGCRecHitWaferU" ,&rechitWaferU);
            lRecTree->SetBranchAddress("HGCRecHitWaferV" ,&rechitWaferV);
            lRecTree->SetBranchAddress("HGCRecHitCellU"  ,&rechitCellU);
            lRecTree->SetBranchAddress("HGCRecHitCellV"  ,&rechitCellV);
            lRecTree->SetBranchAddress("HGCSimHitsEnergy",&simhitEnergy);
            lRecTree->SetBranchAddress("HGCSimHitsEta"   ,&simhitEta);
            lRecTree->SetBranchAddress("HGCSimHitsPhi"   ,&simhitPhi);
            lRecTree->SetBranchAddress("HGCSimHitsPosx"  ,&simhitPosx);
            lRecTree->SetBranchAddress("HGCSimHitsPosy"  ,&simhitPosy);
            lRecTree->SetBranchAddress("HGCSimHitsPosz"  ,&simhitPosz);
            lRecTree->SetBranchAddress("HGCSimHitsLayer" ,&simhitLayer);
            lRecTree->SetBranchAddress("HGCSimHitsIndex" ,&simhitIndex);
            lRecTree->SetBranchAddress("HGCSimHitsWaferU",&simhitWaferU);
            lRecTree->SetBranchAddress("HGCSimHitsWaferV",&simhitWaferV);
            lRecTree->SetBranchAddress("HGCSimHitsCellU" ,&simhitCellU);
            lRecTree->SetBranchAddress("HGCSimHitsCellV" ,&simhitCellV);
            lRecTree->SetBranchAddress("GenParEta"       ,&genEta);
            lRecTree->SetBranchAddress("GenParPhi"       ,&genPhi);
        }

        lRecTree->GetEntry(ievtRec);

        double etagen   = 99999.;
        double phigen   = 99999.;
        //double thetagen = -1.;
        if((*genEta).size()>0) {
            etagen   = (*genEta)[0];
            phigen   = (*genPhi)[0];
        }

        if (debug) std::cout << " - Event contains " << (*rechitEnergy).size()
        << " rechits." << std::endl;
        double coneSize = 0.3;

        std::set<std::tuple<int, int, int, int, int>> saturatedList1;
        std::set<std::tuple<int, int, int, int, int>> saturatedList2;
        float rechitsum = 0;
        float rechitsum_corr = 0;

        // Loop over rechits of event
        for (unsigned iH(0); iH<(*rechitEnergy).size(); ++iH){
            unsigned layer   = (*rechitLayer)[iH];
            double   zh      = (*rechitPosz)[iH];
            double   lenergy = (*rechitEnergy)[iH];
            double   leta    = (*rechitEta)[iH];
            double   lphi    = (*rechitPhi)[iH];
            double   dR      = DeltaR(etagen,phigen,leta,lphi);

            int waferU  = (*rechitWaferU)[iH];
            int waferV  = (*rechitWaferV)[iH];
            unsigned cellU   = (*rechitCellU)[iH];
            unsigned cellV   = (*rechitCellV)[iH];
            unsigned index   = (*rechitIndex)[iH];

            /* Select hits that are:
            **     - in CE-E
            **     - within DeltaR < 0.3 wrt gen particle
            **     - in positive endcap
            */
            if(!index && zh > 0 && dR < coneSize) {
                if(lenergy>27.7 && lenergy<27.85){
                    // Format: (layer, waferU, waferV, cellU, cellV)
                    std::tuple<int, int, int, int, int> saturatedCell;
                    std::get<0>(saturatedCell) = layer;
                    std::get<1>(saturatedCell) = waferU;
                    std::get<2>(saturatedCell) = waferV;
                    std::get<3>(saturatedCell) = cellU;
                    std::get<4>(saturatedCell) = cellV;
                    saturatedList1.insert(saturatedCell);
                }
                else if(lenergy>41.3 && lenergy<41.45){
                    // Format: (layer, waferU, waferV, cellU, cellV)
                    std::tuple<int, int, int, int, int> saturatedCell;
                    std::get<0>(saturatedCell) = layer;
                    std::get<1>(saturatedCell) = waferU;
                    std::get<2>(saturatedCell) = waferV;
                    std::get<3>(saturatedCell) = cellU;
                    std::get<4>(saturatedCell) = cellV;
                    saturatedList2.insert(saturatedCell);
                }
                else {
                    rechitsum += lenergy;
                    rechitsum_corr += lenergy;
                }
            }
        }

        // Loop over simhits of event
        for (unsigned iH(0); iH<(*simhitEnergy).size(); ++iH){
            unsigned layer   = (*simhitLayer)[iH];
            double   zh      = (*simhitPosz)[iH];
            double   lenergy = (*simhitEnergy)[iH];
            double   leta    = (*simhitEta)[iH];
            double   lphi    = (*simhitPhi)[iH];
            double   dR      = DeltaR(etagen,phigen,leta,lphi);

            int waferU  = (*simhitWaferU)[iH];
            int waferV  = (*simhitWaferV)[iH];
            unsigned cellU   = (*simhitCellU)[iH];
            unsigned cellV   = (*simhitCellV)[iH];
            unsigned index   = (*simhitIndex)[iH];

            /* Select hits that are:
            **     - in CE-E
            **     - within DeltaR < 0.3 wrt gen particle
            **     - in positive endcap
            */
            if(!index && zh > 0 && dR < coneSize) {
                std::tuple<int, int, int, int, int> tempCell;
                std::get<0>(tempCell) = layer;
                std::get<1>(tempCell) = waferU;
                std::get<2>(tempCell) = waferV;
                std::get<3>(tempCell) = cellU;
                std::get<4>(tempCell) = cellV;
                std::set<std::tuple<int, int, int, int, int>>::iterator ibc1 = saturatedList1.find(tempCell);
                std::set<std::tuple<int, int, int, int, int>>::iterator ibc2 = saturatedList2.find(tempCell);
                if(ibc1 != saturatedList1.end()){
                    rechitsum_corr += lenergy*122.251;
                }
                if(ibc2 != saturatedList2.end()){
                    rechitsum_corr += lenergy*182.138;
                }
            }
        }

        h1->Fill(rechitsum);
        h2->Fill(rechitsum_corr);
        ievtRec++;
    }

    if(debug) std::cout << "Writing files ..." << std::endl;
    outputFile->cd();
    outputFile->Write();
    outputFile->Close();

    return 0;
}
