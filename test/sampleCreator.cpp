/****************************************************
**  Channelog:
**      Moved to CMSSW output.
**      Editor: Christos Papageorgakis
**
**      This code uses ntuples to produce rechitsums
**      for a given saturated Si cell fraction.
**      The output also contains an ntuple of all the
**      saturated cells that can be used to train a DNN to
**      estimate the lost energy.
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

/* Function that gives a vector of tuples that describe the neighbors
** of the tuple that is given as the input.
*/
std::vector<std::tuple<int, int, int, int, int, int>> getNeighbors(
    std::tuple<int, int, int, int, int, int> saturatedCell)
{
    std::vector<std::tuple<int, int, int, int, int, int>> neighbors;
    // Find same-layer neighboring cells
    // cell ( 0,-1) wrt given
    std::tuple<int, int, int, int, int, int> n1(saturatedCell);
    std::get<4>(n1) -= 1;
    // cell (-1,-1) wrt given
    std::tuple<int, int, int, int, int, int> n2(saturatedCell);
    std::get<3>(n2) -= 1;
    std::get<4>(n2) -= 1;
    // cell (-1, 0) wrt given
    std::tuple<int, int, int, int, int, int> n3(saturatedCell);
    std::get<3>(n3) -= 1;
    // cell ( 0,+1) wrt given
    std::tuple<int, int, int, int, int, int> n4(saturatedCell);
    std::get<4>(n4) += 1;
    // cell (+1, 0) wrt given
    std::tuple<int, int, int, int, int, int> n5(saturatedCell);
    std::get<3>(n5) += 1;
    std::get<4>(n5) += 1;
    // cell (+1,+1) wrt given
    std::tuple<int, int, int, int, int, int> n6(saturatedCell);
    std::get<3>(n6) += 1;

    // Check boundary conditions and make transitions between wafers when on the edge
    // For n1
    if (std::get<3>(n1) > -1 && std::get<3>(n1) < 8 && std::get<4>(n1) == -1){
        std::get<1>(n1) += 1;
        std::get<3>(n1) += 8;
        std::get<4>(n1) = 15;
    }else if (std::get<3>(n1)-std::get<4>(n1) == 9){
        std::get<1>(n1) += 1;
        std::get<2>(n1) += 1;
        std::get<3>(n1) -= 8;
        std::get<4>(n1) += 8;
    }
    // For n2
    if (std::get<3>(n2) > -1 && std::get<3>(n2) < 8 && std::get<4>(n2) == -1){
        std::get<1>(n2) += 1;
        std::get<3>(n2) += 8;
        std::get<4>(n2) = 15;
    }else if (std::get<4>(n2) > -1 && std::get<4>(n2) < 8 && std::get<3>(n2) == -1){
        std::get<2>(n2) -= 1;
        std::get<3>(n2) = 15;
        std::get<4>(n2) += 8;
    }
    // For n3
    if (std::get<4>(n3) > -1 && std::get<4>(n3) < 8 && std::get<3>(n3) == -1){
        std::get<2>(n3) -= 1;
        std::get<3>(n3) = 15;
        std::get<4>(n3) += 8;
    }else if (std::get<4>(n3)-std::get<3>(n3) == 8){
        std::get<1>(n3) -= 1;
        std::get<2>(n3) -= 1;
        std::get<3>(n3) += 8;
        std::get<4>(n3) -= 8;
    }
    // For n4
    if (std::get<3>(n4) > 7 && std::get<3>(n4) < 16 && std::get<4>(n4) == 16){
        std::get<1>(n4) -= 1;
        std::get<3>(n4) -= 8;
        std::get<4>(n4) = 0;
    }else if (std::get<4>(n4)-std::get<3>(n4) == 8){
        std::get<1>(n4) -= 1;
        std::get<2>(n4) -= 1;
        std::get<3>(n4) += 8;
        std::get<4>(n4) -= 8;
    }
    // For n5
    if (std::get<4>(n5) > 7 && std::get<4>(n5) < 16 && std::get<3>(n5) == 16){
        std::get<2>(n5) += 1;
        std::get<3>(n5) = 0;
        std::get<4>(n5) -= 8;
    }else if (std::get<3>(n5) > 7 && std::get<3>(n5) < 16 && std::get<4>(n5) == 16){
        std::get<1>(n5) -= 1;
        std::get<3>(n5) -= 8;
        std::get<4>(n5) = 0;
    }
    // For n6
    if (std::get<4>(n6) > 7 && std::get<4>(n6) < 16 && std::get<3>(n6) == 16){
        std::get<2>(n6) += 1;
        std::get<3>(n6) = 0;
        std::get<4>(n6) -= 8;
    }else if (std::get<3>(n6)-std::get<4>(n6) == 9){
        std::get<1>(n6) += 1;
        std::get<2>(n6) += 1;
        std::get<3>(n6) -= 8;
        std::get<4>(n6) += 8;
    }

    neighbors.push_back(n1);
    neighbors.push_back(n2);
    neighbors.push_back(n3);
    neighbors.push_back(n4);
    neighbors.push_back(n5);
    neighbors.push_back(n6);
    return neighbors;
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
    unsigned nRuns,firstRun;
    std::string recoFileName;
    unsigned debug;
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
    ("firstRun",        po::value<unsigned>(&firstRun)->default_value(1))
    ("debug,d",         po::value<unsigned>(&debug)->default_value(0))
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
        for (unsigned i = firstRun;i<=nRuns+firstRun;++i){
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

    /**********************************
    ** ML Study output section
    **     - MLlayer is saturated cell layer
    **     - MLeta is gen eta
    **     - MLphi is gen phi
    **     - MLni is ith saturated cell neighbor
    **     - MLuni is 6 neighbors at layer+1
    **     - MLdni is 6 neighbors at layer-1
    **     - MLsaturated is saturated cell rechit
    ** We also need to create a ?set? container to store the values before writing to TTree
    **********************************/
    float MLlayer, MLeta, MLphi, MLsaturated, MLnup, MLndown, MLevent;
    float MLwaferU, MLwaferV, MLcellU, MLcellV;
    float MLn1, MLn2, MLn3, MLn4, MLn5, MLn6;
    float MLdn1, MLdn2, MLdn3, MLdn4, MLdn5, MLdn6;
    float MLun1, MLun2, MLun3, MLun4, MLun5, MLun6;
    float MLrechitsum, MLsimHits, cellType;
    TTree* t1 = new TTree("t1","sample");
    t1->Branch("MLlayer"     ,&MLlayer     ,"MLlayer/F"     );
    t1->Branch("MLwaferU"    ,&MLwaferU    ,"MLwaferU/F"    );
    t1->Branch("MLwaferV"    ,&MLwaferV    ,"MLwaferV/F"    );
    t1->Branch("MLcellU"     ,&MLcellU     ,"MLcellU/F"     );
    t1->Branch("MLcellV"     ,&MLcellV     ,"MLcellV/F"     );
    t1->Branch("MLeta"       ,&MLeta       ,"MLeta/F"       );
    t1->Branch("MLphi"       ,&MLphi       ,"MLphi/F"       );
    t1->Branch("MLn1"        ,&MLn1        ,"MLn1/F"        );
    t1->Branch("MLn2"        ,&MLn2        ,"MLn2/F"        );
    t1->Branch("MLn3"        ,&MLn3        ,"MLn3/F"        );
    t1->Branch("MLn4"        ,&MLn4        ,"MLn4/F"        );
    t1->Branch("MLn5"        ,&MLn5        ,"MLn5/F"        );
    t1->Branch("MLn6"        ,&MLn6        ,"MLn6/F"        );
    t1->Branch("MLsaturated" ,&MLsaturated ,"MLsaturated/F" );
    t1->Branch("MLnup"       ,&MLnup       ,"MLnup/F"       );
    t1->Branch("MLndown"     ,&MLndown     ,"MLndown/F"     );
    t1->Branch("MLun1"       ,&MLun1       ,"MLun1/F"       );
    t1->Branch("MLun2"       ,&MLun2       ,"MLun2/F"       );
    t1->Branch("MLun3"       ,&MLun3       ,"MLun3/F"       );
    t1->Branch("MLun4"       ,&MLun4       ,"MLun4/F"       );
    t1->Branch("MLun5"       ,&MLun5       ,"MLun5/F"       );
    t1->Branch("MLun6"       ,&MLun6       ,"MLun6/F"       );
    t1->Branch("MLdn1"       ,&MLdn1       ,"MLdn1/F"       );
    t1->Branch("MLdn2"       ,&MLdn2       ,"MLdn2/F"       );
    t1->Branch("MLdn3"       ,&MLdn3       ,"MLdn3/F"       );
    t1->Branch("MLdn4"       ,&MLdn4       ,"MLdn4/F"       );
    t1->Branch("MLdn5"       ,&MLdn5       ,"MLdn5/F"       );
    t1->Branch("MLdn6"       ,&MLdn6       ,"MLdn6/F"       );
    t1->Branch("MLevent"     ,&MLevent     ,"MLevent/F"     );
    t1->Branch("MLrechitsum" ,&MLrechitsum ,"MLrechitsum/F" );
    t1->Branch("MLsimHits"   ,&MLsimHits   ,"MLsimHits/F"   );
    t1->Branch("cellType"    ,&cellType    ,"cellType/F"    );

    // Format:
    // <layer, waferU, waferV, cellU, cellV, cellType>
    // cellType is 0 for 300um and 1 for 200um
    std::set<std::tuple<int, int, int, int, int, int>> saturatedList;

    // Define average energy in layers plus and minus 1
    std::set<std::tuple<int, int, int, int, int, int>> adj_to_saturated;
    std::set<std::tuple<int, int, int, int, int, int>> adj_to_saturated_inlay;

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
        /*
        ** Define a vector of the array:
        ** {saturated cell:
        **      layer, waferU, waferV, cellU, cellV,
        **      eta, phi,
        **      MLn1, MLn2, MLn3, MLn4, MLn5, MLn6,
        **      recHit,
        **      MLup, MLdown,
        **      MLun1, MLun2, MLun3, MLun4, MLun5, MLun6,
        **      MLdn1, MLdn2, MLdn3, MLdn4, MLdn5, MLdn6,
        **      MLevent,
        **      MLrechitsum,
        **      MLsimHits,
        **      cellType
        ** }
        */
        std::vector<std::array<float, 32>> MLvectorev;

        if (ievtRec>=lRecTree->GetEntries()) continue;
        Long64_t local_entry = lRecTree->LoadTree(ievt);

        if (debug) std::cout << std::endl<<std::endl << "... Processing entry: " << ievt << std::endl;
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
        if((*genEta).size()>0) {
            etagen   = (*genEta)[0];
            phigen   = (*genPhi)[0];
        }

        if (debug) std::cout << " - Event contains " << (*rechitEnergy).size()
        << " rechits." << std::endl;
        double coneSize = 0.3;

        // First loop over rechits of event
        for (unsigned iH(0); iH<(*rechitEnergy).size(); ++iH){
            int      layer   = (*rechitLayer)[iH];
            float   zh      = (*rechitPosz)[iH];
            float   lenergy = (*rechitEnergy)[iH];
            float   leta    = (*rechitEta)[iH];
            float   lphi    = (*rechitPhi)[iH];
            float   dR      = DeltaR(etagen,phigen,leta,lphi);

            int waferU  = (*rechitWaferU)[iH];
            int waferV  = (*rechitWaferV)[iH];
            int cellU   = (*rechitCellU)[iH];
            int cellV   = (*rechitCellV)[iH];
            int index   = (*rechitIndex)[iH];

            /* Select hits that are:
            **     - in CE-E
            **     - within DeltaR < 0.3 wrt gen particle
            **     - in positive endcap
            */
            if(!index && zh > 0 && dR < coneSize) {
                if(lenergy>27.7 && lenergy<27.85){
                    // Format: (layer, waferU, waferV, cellU, cellV, cellType)
                    std::tuple<int, int, int, int, int, int> saturatedCell;
                    std::get<0>(saturatedCell) = layer;
                    std::get<1>(saturatedCell) = waferU;
                    std::get<2>(saturatedCell) = waferV;
                    std::get<3>(saturatedCell) = cellU;
                    std::get<4>(saturatedCell) = cellV;
                    std::get<5>(saturatedCell) = 0;
                    saturatedList.insert(saturatedCell);
                    std::array<float, 32> tempArr = {
                        (float)layer, (float)waferU, (float)waferV,
                        (float)cellU, (float)cellV, (float)leta, (float)lphi,
                        0, 0, 0, 0, 0, 0, // n1, n2, n3, n4, n5, n6
                        lenergy,
                        0, 0,             // nup, ndown
                        0, 0, 0, 0, 0, 0, // un1, un2, un3, un4, un5, un6
                        0, 0, 0, 0, 0, 0, // dn1, dn2, dn3, dn4, dn5, dn6
                        (float)ievt,
                        0, 0,             // recHitsum, simhits
                        0                 // cellType
                    };
                    MLvectorev.push_back(tempArr);
                }else if(lenergy>41.3 && lenergy<41.45){
                    // Format: (layer, waferU, waferV, cellU, cellV, cellType)
                    std::tuple<int, int, int, int, int, int> saturatedCell;
                    std::get<0>(saturatedCell) = layer;
                    std::get<1>(saturatedCell) = waferU;
                    std::get<2>(saturatedCell) = waferV;
                    std::get<3>(saturatedCell) = cellU;
                    std::get<4>(saturatedCell) = cellV;
                    std::get<5>(saturatedCell) = 1;
                    saturatedList.insert(saturatedCell);
                    std::array<float, 32> tempArr = {
                        (float)layer, (float)waferU, (float)waferV,
                        (float)cellU, (float)cellV, (float)leta, (float)lphi,
                        0, 0, 0, 0, 0, 0, // n1, n2, n3, n4, n5, n6
                        lenergy,
                        0, 0,             // nup, ndown
                        0, 0, 0, 0, 0, 0, // un1, un2, un3, un4, un5, un6
                        0, 0, 0, 0, 0, 0, // dn1, dn2, dn3, dn4, dn5, dn6
                        (float)ievt,
                        0, 0,             // recHitsum, simhits
                        1                 // cellType
                    };
                    MLvectorev.push_back(tempArr);
                }
            }
        }

        // Insert loops for adj_to_saturated and adj_to_saturated_inlay population here
        for (auto itr = saturatedList.begin(); itr != saturatedList.end(); ++itr) {
            adj_to_saturated.insert(
                {
                    0, //corresponds to cell bellow
                    std::get<0>(*itr)-1,
                    std::get<1>(*itr),
                    std::get<2>(*itr),
                    std::get<3>(*itr),
                    std::get<4>(*itr)
                }
            );
            adj_to_saturated.insert(
                {
                    1, //corresponds to cell above
                    std::get<0>(*itr)+1,
                    std::get<1>(*itr),
                    std::get<2>(*itr),
                    std::get<3>(*itr),
                    std::get<4>(*itr)
                }
            );

            std::vector<std::tuple<int,int,int,int,int,int>> inLayerNeighbors;
            inLayerNeighbors = getNeighbors(*itr);
            int iN = 0;
            for(auto itr2 = inLayerNeighbors.begin(); itr2!=inLayerNeighbors.end(); ++itr2){
                adj_to_saturated_inlay.insert(
                    {
                        iN,
                        std::get<0>(*itr2),
                        std::get<1>(*itr2),
                        std::get<2>(*itr2),
                        std::get<3>(*itr2),
                        std::get<4>(*itr2)
                    }
                );
                ++iN;
            }
        }

        MLrechitsum = 0;
        // Second loop over rechits of event
        for (unsigned iH(0); iH<(*rechitEnergy).size(); ++iH){
            int      layer   = (*rechitLayer)[iH];
            double   zh      = (*rechitPosz)[iH];
            double   lenergy = (*rechitEnergy)[iH];
            double   leta    = (*rechitEta)[iH];
            double   lphi    = (*rechitPhi)[iH];
            double   dR      = DeltaR(etagen,phigen,leta,lphi);

            int waferU  = (*rechitWaferU)[iH];
            int waferV  = (*rechitWaferV)[iH];
            int cellU   = (*rechitCellU)[iH];
            int cellV   = (*rechitCellV)[iH];
            int index   = (*rechitIndex)[iH];

            /* Select hits that are:
            **     - in CE-E
            **     - within DeltaR < 0.3 wrt gen particle
            **     - in positive endcap
            */
            if(!index && zh > 0 && dR < coneSize) {
                std::tuple<int, int, int, int, int, int> tempsi1(layer,waferU,waferV,cellU,cellV,0);
                std::tuple<int, int, int, int, int, int> tempsi2(layer,waferU,waferV,cellU,cellV,1);
                std::set<std::tuple<int, int, int, int, int, int>>::iterator ibc1=saturatedList.find(tempsi1);
                std::set<std::tuple<int, int, int, int, int, int>>::iterator ibc2=saturatedList.find(tempsi2);

                // Calculate energy without saturated Si cells
                if(ibc1 == saturatedList.end() && ibc2 == saturatedList.end()) {
                    MLrechitsum += lenergy;
                }

                /* Perform Simple Average
                ** First, check if the cell is in a neighbors list
                */
                std::tuple<int, int, int, int, int, int> tempsiU(
                    1,layer,waferU,waferV,cellU,cellV
                );
                std::tuple<int, int, int, int, int, int> tempsiD(
                    0,layer,waferU,waferV,cellU,cellV
                );
                std::set<std::tuple<int, int, int, int, int, int>>::iterator itrU=adj_to_saturated.find(tempsiU);
                std::set<std::tuple<int, int, int, int, int, int>>::iterator itrD=adj_to_saturated.find(tempsiD);

                if(itrU!=adj_to_saturated.end()) {
                    for(auto itr = MLvectorev.begin(); itr != MLvectorev.end(); itr++) {
                        if( (*itr)[0] == layer-1 &&
                            (*itr)[1] == waferU && (*itr)[2] == waferV &&
                            (*itr)[3] == cellU  && (*itr)[4] == cellV
                        ){
                            (*itr)[14] = lenergy;
                        }
                    }
                }
                if(itrD!=adj_to_saturated.end()) {
                    for(auto itr = MLvectorev.begin(); itr != MLvectorev.end(); itr++) {
                        if( (*itr)[0] == layer+1 &&
                            (*itr)[1] == waferU && (*itr)[2] == waferV &&
                            (*itr)[3] == cellU  && (*itr)[4] == cellV
                        ){
                            (*itr)[15] = lenergy;
                        }
                    }
                }

                // Get rechits of saturated cells' neighbors
                for(int n = 0; n < 6; ++n){
                    // Same layer neighbors
                    std::tuple<int, int, int, int, int, int> tempsiNn(
                        n,layer,waferU,waferV,cellU,cellV
                    );
                    //Check if cell is an nth neighbor of some saturated cell
                    std::set<std::tuple<int, int, int, int, int, int>>::iterator itrNn=adj_to_saturated_inlay.find(tempsiNn);
                    if(itrNn!=adj_to_saturated_inlay.end()) {
                        std::vector<std::tuple<int,int,int,int,int, int>> sameLayerNeighbors;
                        sameLayerNeighbors = getNeighbors(tempsi1);
                        // Get neighbor number
                        int nn = (std::get<0>(*itrNn)+3)%6;
                        std::tuple<int, int, int, int> saturatedCell;
                        std::get<0>(saturatedCell) = std::get<1>(sameLayerNeighbors[nn]);
                        std::get<1>(saturatedCell) = std::get<2>(sameLayerNeighbors[nn]);
                        std::get<2>(saturatedCell) = std::get<3>(sameLayerNeighbors[nn]);
                        std::get<3>(saturatedCell) = std::get<4>(sameLayerNeighbors[nn]);
                        for(auto itr = MLvectorev.begin(); itr != MLvectorev.end(); itr++) {
                            if( (*itr)[0] == layer &&
                                (*itr)[1] == std::get<0>(saturatedCell) && (*itr)[2] == std::get<1>(saturatedCell) &&
                                (*itr)[3] == std::get<2>(saturatedCell) && (*itr)[4] == std::get<3>(saturatedCell)
                            ){
                                (*itr)[n+7] = lenergy;
                            }
                        }
                    }

                    // Next layer neighbors
                    std::tuple<int, int, int, int, int, int> tempsiUNn(
                        n,layer-1,waferU,waferV,cellU,cellV
                    );
                    std::set<std::tuple<int, int, int, int, int, int>>::iterator itrUNn=adj_to_saturated_inlay.find(tempsiUNn);
                    if(itrUNn!=adj_to_saturated_inlay.end()) {
                        std::vector<std::tuple<int,int,int,int,int, int>> nextLayerNeighbors;
                        nextLayerNeighbors = getNeighbors(tempsi1);
                        // Get neighbor number
                        int nn = (std::get<0>(*itrUNn)+3)%6;
                        std::tuple<int, int, int, int> saturatedCell;
                        std::get<0>(saturatedCell) = std::get<1>(nextLayerNeighbors[nn]);
                        std::get<1>(saturatedCell) = std::get<2>(nextLayerNeighbors[nn]);
                        std::get<2>(saturatedCell) = std::get<3>(nextLayerNeighbors[nn]);
                        std::get<3>(saturatedCell) = std::get<4>(nextLayerNeighbors[nn]);
                        for(auto itr = MLvectorev.begin(); itr != MLvectorev.end(); itr++) {
                            if( (*itr)[0] == layer-1 &&
                            (*itr)[1] == std::get<0>(saturatedCell) && (*itr)[2] == std::get<1>(saturatedCell) &&
                            (*itr)[3] == std::get<2>(saturatedCell) && (*itr)[4] == std::get<3>(saturatedCell)
                            ){
                                (*itr)[n+16] = lenergy;
                            }
                        }
                    }

                    // Previous layer neighbors
                    std::tuple<int, int, int, int, int, int> tempsiDNn(
                        n,layer+1,waferU,waferV,cellU,cellV
                    );
                    std::set<std::tuple<int, int, int, int, int, int>>::iterator itrDNn=adj_to_saturated_inlay.find(tempsiDNn);
                    if(itrDNn!=adj_to_saturated_inlay.end()) {
                        std::vector<std::tuple<int,int,int,int,int, int>> prevLayerNeighbors;
                        prevLayerNeighbors = getNeighbors(tempsi1);
                        // Get neighbor number
                        int nn = (std::get<0>(*itrDNn)+3)%6;
                        std::tuple<int, int, int, int> saturatedCell;
                        std::get<0>(saturatedCell) = std::get<1>(prevLayerNeighbors[nn]);
                        std::get<1>(saturatedCell) = std::get<2>(prevLayerNeighbors[nn]);
                        std::get<2>(saturatedCell) = std::get<3>(prevLayerNeighbors[nn]);
                        std::get<3>(saturatedCell) = std::get<4>(prevLayerNeighbors[nn]);
                        for(auto itr = MLvectorev.begin(); itr != MLvectorev.end(); itr++) {
                            if( (*itr)[0] == layer+1 &&
                            (*itr)[1] == std::get<0>(saturatedCell) && (*itr)[2] == std::get<1>(saturatedCell) &&
                            (*itr)[3] == std::get<2>(saturatedCell) && (*itr)[4] == std::get<3>(saturatedCell)
                            ){
                                (*itr)[n+22] = lenergy;
                            }
                        }
                    }
                }
            }
        }

        // Loop over simhits of event
        for (unsigned iH(0); iH<(*simhitEnergy).size(); ++iH){
            int   layer   = (*simhitLayer)[iH];
            float zh      = (*simhitPosz)[iH];
            float lenergy = (*simhitEnergy)[iH];
            float leta    = (*simhitEta)[iH];
            float lphi    = (*simhitPhi)[iH];
            float dR      = DeltaR(etagen,phigen,leta,lphi);

            int waferU = (*simhitWaferU)[iH];
            int waferV = (*simhitWaferV)[iH];
            int cellU  = (*simhitCellU)[iH];
            int cellV  = (*simhitCellV)[iH];
            int index  = (*simhitIndex)[iH];

            /* Select hits that are:
            **     - in CE-E
            **     - within DeltaR < 0.3 wrt gen particle
            **     - in positive endcap
            */
            if(!index && zh > 0 && dR < coneSize) {
                for(auto itr = MLvectorev.begin(); itr != MLvectorev.end(); ++itr) {
                    if(
                        (*itr)[0] == layer  && (*itr)[1] == waferU &&
                        (*itr)[2] == waferV && (*itr)[3] == cellU  &&
                        (*itr)[4] == cellV
                    ){
                        (*itr)[30] += lenergy;
                    }
                }
            }
        }

        //Export the ML dataset values to the TTree
        for(auto itr = MLvectorev.begin(); itr != MLvectorev.end(); ++itr) {
            if ((*itr)[5] > 0 || (*itr)[0]==-1) {
                /* This condition is necessary to ensure the cell was within
                ** the cone.
                */
                MLlayer     = (*itr)[0];
                MLwaferU    = (*itr)[1];
                MLwaferV    = (*itr)[2];
                MLcellU     = (*itr)[3];
                MLcellV     = (*itr)[4];
                MLeta       = (*itr)[5];
                MLphi       = (*itr)[6];
                MLn1        = (*itr)[7];
                MLn2        = (*itr)[8];
                MLn3        = (*itr)[9];
                MLn4        = (*itr)[10];
                MLn5        = (*itr)[11];
                MLn6        = (*itr)[12];
                MLsaturated = (*itr)[13];
                MLnup       = (*itr)[14];
                MLndown     = (*itr)[15];
                MLun1       = (*itr)[16];
                MLun2       = (*itr)[17];
                MLun3       = (*itr)[18];
                MLun4       = (*itr)[19];
                MLun5       = (*itr)[20];
                MLun6       = (*itr)[21];
                MLdn1       = (*itr)[22];
                MLdn2       = (*itr)[23];
                MLdn3       = (*itr)[24];
                MLdn4       = (*itr)[25];
                MLdn5       = (*itr)[26];
                MLdn6       = (*itr)[27];
                MLevent     = (float)ievt;
                //MLrechitsum = rechitsumsaturated_Si;
                MLsimHits   = (*itr)[30];
                cellType    = (*itr)[31];
                t1->Fill();
            }
        }
        ievtRec++;
    }

    if(debug) std::cout << "Writing files ..." << std::endl;
    outputFile->cd();
    outputFile->Write();
    outputFile->Close();

    return 0;
}
