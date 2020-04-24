#include "RootStyle.cc"

bool isSaturated(float rechit);

std::vector<float> averageCalculator(int En);

void averageSaturatedNeighbors(){
    Float_t En[]   = {550,750,1000,1400,2000,2800};
    Float_t nn6[]  = { 0., 0.,  0.,  0.,  0.,  0.};
    Float_t nn2[]  = { 0., 0.,  0.,  0.,  0.,  0.};
    Float_t nn8[]  = { 0., 0.,  0.,  0.,  0.,  0.};
    Float_t nn20[] = { 0., 0.,  0.,  0.,  0.,  0.};

    for (int i = 0; i < 6; ++i){
        std::vector<float> tempVec = averageCalculator(En[i]);
        nn6[i]  = tempVec[0];
        nn2[i]  = tempVec[1];
        nn8[i]  = tempVec[2];
        nn20[i] = tempVec[3];
    }

    set_root_style();
    //set_tdr_style();
    //set_vasu_style();
    //setTDRStyle();

    TCanvas* c = new TCanvas("c","c",1);
    TGraph* gr6  = new TGraph(6, En,  nn6);
    TGraph* gr2  = new TGraph(6, En,  nn2);
    TGraph* gr8  = new TGraph(6, En,  nn8);
    TGraph* gr20 = new TGraph(6, En, nn20);
    gr6->SetMarkerColor(kBlue);
    gr6->SetMarkerStyle(21);
    gr2->SetMarkerColor(kRed);
    gr2->SetMarkerStyle(21);
    gr8->SetMarkerColor(kGreen);
    gr8->SetMarkerStyle(21);
    gr20->SetMarkerColor(kBlack);
    gr20->SetMarkerStyle(21);

    TMultiGraph* gr = new TMultiGraph();
    gr->Add(gr6);
    gr->Add(gr2);
    gr->Add(gr8);
    gr->Add(gr20);
    gr->SetTitle(";Energy [GeV]; <N>");
    gr->Draw("AP");

    TLegend* legend = new TLegend(0.12,0.6,0.6,0.88);
    legend->SetBorderSize(0);
    legend->AddEntry(gr6,"6 closest neighbors (same layer)","p");
    legend->AddEntry(gr2,"2 closest neighbors (adjacent layers)","p");
    legend->AddEntry(gr8,"8 closest neighbors (any layer)","p");
    legend->AddEntry(gr20,"20 closest neighbors (any layer)","p");
    legend->Draw();
}

std::vector<float> averageCalculator(int En){
    TString filename = "RegressionResults2/flatRegressionResult_"+to_string(En)+"GeV.root";
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

    int n = t1->GetEntries();
    int Ntot = 0;
    int nn = 0;
    int dn = 0;
    int un = 0;
    int pnn = 0;
    for(int i = 0; i <= n; ++i){
        t1->GetEntry(i);
        if(layer > 0) Ntot++;
        if(isSaturated(n1)) nn++;
        if(isSaturated(n2)) nn++;
        if(isSaturated(n3)) nn++;
        if(isSaturated(n4)) nn++;
        if(isSaturated(n5)) nn++;
        if(isSaturated(n6)) nn++;
        if(isSaturated(dn1)) dn++;
        if(isSaturated(dn2)) dn++;
        if(isSaturated(dn3)) dn++;
        if(isSaturated(dn4)) dn++;
        if(isSaturated(dn5)) dn++;
        if(isSaturated(dn6)) dn++;
        if(isSaturated(un1)) un++;
        if(isSaturated(un2)) un++;
        if(isSaturated(un3)) un++;
        if(isSaturated(un4)) un++;
        if(isSaturated(un5)) un++;
        if(isSaturated(un6)) un++;
        if(isSaturated(ndown)) pnn++;
        if(isSaturated(nup))   pnn++;
    }

    /*
    ** output vector contains:
    **     - <N> for 6 closest neighbors (same layer)
    **     - <N> for 2 closest neighbors (other layers)
    **     - <N> for 8 closest neighbors (any layer)
    **     - <N> for 20 closest neighbors (any layer)
    */
    float an6, an2, an8, an20;
    an6  = (float)nn/(float)Ntot;
    an2  = (float)pnn/(float)Ntot;
    an8  = an2+an6;
    an20 = an8+((float)un+(float)dn)/(float)Ntot;

    std::vector<float> output;
    output.push_back(an6);
    output.push_back(an2);
    output.push_back(an8);
    output.push_back(an20);

    cout << "<N_samelayer> \t=\t" << (float)nn/(float)Ntot  << endl;
    cout << "<N_pncells> \t=\t"   << (float)pnn/(float)Ntot << endl;
    cout << "<N_nextlayer> \t=\t" << (float)un/(float)Ntot  << endl;
    cout << "<N_prevlayer> \t=\t" << (float)dn/(float)Ntot  << endl;
    cout << "<N> \t\t=\t" << ((float)nn+(float)un+(float)dn+(float)pnn)/(float)Ntot << endl;

    return output;
}

bool isSaturated(float rechit) {
    if(rechit>27.7 && rechit<27.85){
        return 1;
    }
    else if(rechit>41.3 && rechit<41.45){
        return 1;
    }
    return 0;
}
