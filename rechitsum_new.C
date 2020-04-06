std::vector<double> rechitsum_new(int En, int bins, int range, double fit_cut){
    TString filename = "RegressionResults/flatRegressionResult_"+to_string(En)+"GeV.root";
    TFile* fin = TFile::Open(filename);
    TTree* t1 = dynamic_cast< TTree* >(fin->Get("t1"));
    Float_t rechitsum, truth, regression, event, nup, ndown;
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


    TString histname = "single gamma "+to_string(En)+"GeV;recHitSum [GeV];Entries";
    TH1F* h1 = new TH1F("h1",histname,bins,En-range,En+range);
    int n = t1->GetEntries();
    Float_t event_tmp;
    Float_t rechitsum_truth  = 0;
    Float_t rechitsum_MLregr = 0;
    for(int i = 0; i <= n; ++i){
        t1->GetEntry(i);
        if(i == 0) {
            event_tmp = event;
            rechitsum_truth  = rechitsum;
            rechitsum_MLregr = rechitsum;
        }
        if(event_tmp != event) {
            h1->Fill(rechitsum);
            rechitsum_truth  =  rechitsum;
            rechitsum_truth  += truth;
            rechitsum_MLregr =  rechitsum;
            rechitsum_MLregr += regression;
            event_tmp = event;
        }else{
            rechitsum_truth  += truth;
            rechitsum_MLregr += regression;
        }
    }

    h1->Draw();
    TFitResultPtr r = h1->Fit("gaus","Sq","");
    r = h1->Fit("gaus","S","",r->Parameter(1)-fit_cut*r->Parameter(2),r->Parameter(1)+3*r->Parameter(2));

    std::vector<double> output_vector;
    double a0  = r->Parameter(0);
    double a1  = r->Parameter(1);
    double a2  = r->Parameter(2);
    double a0e = r->ParError(0);
    double a1e = r->ParError(1);
    double a2e = r->ParError(2);
    output_vector.push_back(a0);
    output_vector.push_back(a1);
    output_vector.push_back(a2);
    output_vector.push_back(a0e);
    output_vector.push_back(a1e);
    output_vector.push_back(a2e);

    return output_vector;
}
