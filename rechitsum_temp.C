#include "RootStyle.cc"

void rechitsum_temp(){
    set_root_style();
    //set_tdr_style();
    //set_vasu_style();
    //setTDRStyle();

    TString filename = "RegressionResults2/flatRegressionResult_2800GeV.root";
    TFile* fin = TFile::Open(filename);
    TTree* t1 = dynamic_cast< TTree* >(fin->Get("t1"));
    Float_t rechitsum, truth, regression, event, nup, ndown, cellType, layer;
    Float_t  n1,  n2,  n3,  n4,  n5,  n6;
    Float_t un1, un2, un3, un4, un5, un6;
    Float_t dn1, dn2, dn3, dn4, dn5, dn6;
    t1->SetBranchAddress(         "n1",         &n1 );
    t1->SetBranchAddress(         "n2",         &n2 );
    t1->SetBranchAddress(         "n3",         &n3 );
    t1->SetBranchAddress(         "n4",         &n4 );
    t1->SetBranchAddress(         "n5",         &n5 );
    t1->SetBranchAddress(         "n6",         &n6 );
    t1->SetBranchAddress(        "un1",        &un1 );
    t1->SetBranchAddress(        "un2",        &un2 );
    t1->SetBranchAddress(        "un3",        &un3 );
    t1->SetBranchAddress(        "un4",        &un4 );
    t1->SetBranchAddress(        "un5",        &un5 );
    t1->SetBranchAddress(        "un6",        &un6 );
    t1->SetBranchAddress(        "dn1",        &dn1 );
    t1->SetBranchAddress(        "dn2",        &dn2 );
    t1->SetBranchAddress(        "dn3",        &dn3 );
    t1->SetBranchAddress(        "dn4",        &dn4 );
    t1->SetBranchAddress(        "dn5",        &dn5 );
    t1->SetBranchAddress(        "dn6",        &dn6 );
    t1->SetBranchAddress(        "nup",        &nup );
    t1->SetBranchAddress(      "ndown",      &ndown );
    t1->SetBranchAddress(  "rechitsum",  &rechitsum );
    t1->SetBranchAddress(     "rechit",      &truth );
    t1->SetBranchAddress(      "event",      &event );
    t1->SetBranchAddress( "regression", &regression );
    t1->SetBranchAddress(   "cellType",   &cellType );
    t1->SetBranchAddress(      "layer",      &layer );


    TString histname = "single gamma 2800 GeV;recHitSum [GeV];Entries";
    TH1F* h1 = new TH1F("h1",histname,300,1800,3000);
    TH1F* h2 = new TH1F("h2",histname,300,1800,3000);
    TH1F* h3 = new TH1F("h3",histname,300,1800,3000);
    int n = t1->GetEntries();
    Float_t event_tmp;
    Float_t rechitsum_nocorr = 0;
    Float_t rechitsum_truth  = 0;
    Float_t rechitsum_MLregr = 0;
    for(int i = 0; i <= n; ++i){
        t1->GetEntry(i);
        if(i == 0) {
            event_tmp = event;
            rechitsum_nocorr = rechitsum;
            rechitsum_truth  = rechitsum;
            rechitsum_MLregr = rechitsum;
        }

        if(event_tmp != event) {
            h1->Fill(rechitsum_nocorr);
            h2->Fill(rechitsum_truth);
            h3->Fill(rechitsum_MLregr);
            rechitsum_nocorr =  rechitsum;
            if (layer > 0 && cellType < 0.5) rechitsum_nocorr  += 27.77;
            if (layer > 0 && cellType > 0.5) rechitsum_nocorr  += 41.37;
            rechitsum_truth  =  rechitsum;
            rechitsum_truth  += truth;
            rechitsum_MLregr =  rechitsum;
            rechitsum_MLregr += regression;
            event_tmp = event;
        }else{
            rechitsum_truth  += truth;
            rechitsum_MLregr += regression;
            if (layer > 0 && cellType < 0.5) rechitsum_nocorr  += 27.77;
            if (layer > 0 && cellType > 0.5) rechitsum_nocorr  += 41.37;
        }
    }

    THStack* hs = new THStack("hs","hs");
    hs->Add(h1);
    hs->Add(h2);
    hs->Add(h3);
    h2->SetLineColor(kRed);
    h3->SetLineColor(kGreen);
    hs->Draw("nostack");
    hs->SetTitle(";Energy [GeV]; Events");

    auto *legend = new TLegend(0.15,0.6,0.35,0.88);
    legend->SetBorderSize(0);
    legend->AddEntry(h1,"no correction","l");
    legend->AddEntry(h2,"truth estimation","l");
    legend->AddEntry(h3,"ML correction","l");
    legend->Draw();
}
