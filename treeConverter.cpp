/*
** Script:
**    - Renames branch names
**    - Adds rechit truth estimation
**    - Skims the data for cells that don't qualify for training set
**        ~ cells that are on the wafer edges
*/

void treeConverter(TString input){
    TFile* f = TFile::Open(input+".root");
    TTree* t = dynamic_cast< TTree* >(f->Get("t1"));
    Float_t n1, n2, n3, n4, n5, n6, nup, ndown, rechit, layer, rechitsum;
    Float_t un1, un2, un3, un4, un5, un6, dn1, dn2, dn3, dn4, dn5, dn6;
    Float_t saturated, simhits, cellType, event;
    t->SetBranchAddress(     "MLlayer",     &layer );
    t->SetBranchAddress(        "MLn1",        &n1 );
    t->SetBranchAddress(        "MLn2",        &n2 );
    t->SetBranchAddress(        "MLn3",        &n3 );
    t->SetBranchAddress(        "MLn4",        &n4 );
    t->SetBranchAddress(        "MLn5",        &n5 );
    t->SetBranchAddress(        "MLn6",        &n6 );
    t->SetBranchAddress(       "MLun1",       &un1 );
    t->SetBranchAddress(       "MLun2",       &un2 );
    t->SetBranchAddress(       "MLun3",       &un3 );
    t->SetBranchAddress(       "MLun4",       &un4 );
    t->SetBranchAddress(       "MLun5",       &un5 );
    t->SetBranchAddress(       "MLun6",       &un6 );
    t->SetBranchAddress(       "MLdn1",       &dn1 );
    t->SetBranchAddress(       "MLdn2",       &dn2 );
    t->SetBranchAddress(       "MLdn3",       &dn3 );
    t->SetBranchAddress(       "MLdn4",       &dn4 );
    t->SetBranchAddress(       "MLdn5",       &dn5 );
    t->SetBranchAddress(       "MLdn6",       &dn6 );
    t->SetBranchAddress(       "MLnup",       &nup );
    t->SetBranchAddress(     "MLndown",     &ndown );
    t->SetBranchAddress(    "cellType",  &cellType ); // 0 for 300um & 1 for 200um
    t->SetBranchAddress(     "MLevent",     &event );
    t->SetBranchAddress(   "MLsimHits",   &simhits );
    t->SetBranchAddress( "MLrechitsum", &rechitsum );
    t->SetBranchAddress( "MLsaturated", &saturated ); // saturated rechit value

    TFile* fout = new TFile(input+"_converted.root","RECREATE");
    TTree* t1 = new TTree("t1","sample");
    t1->Branch(     "layer",     &layer,     "layer/F" );
    t1->Branch(        "n1",        &n1,        "n1/F" );
    t1->Branch(        "n2",        &n2,        "n2/F" );
    t1->Branch(        "n3",        &n3,        "n3/F" );
    t1->Branch(        "n4",        &n4,        "n4/F" );
    t1->Branch(        "n5",        &n5,        "n5/F" );
    t1->Branch(        "n6",        &n6,        "n6/F" );
    t1->Branch(       "dn1",       &dn1,       "dn1/F" );
    t1->Branch(       "dn2",       &dn2,       "dn2/F" );
    t1->Branch(       "dn3",       &dn3,       "dn3/F" );
    t1->Branch(       "dn4",       &dn4,       "dn4/F" );
    t1->Branch(       "dn5",       &dn5,       "dn5/F" );
    t1->Branch(       "dn6",       &dn6,       "dn6/F" );
    t1->Branch(       "un1",       &un1,       "un1/F" );
    t1->Branch(       "un2",       &un2,       "un2/F" );
    t1->Branch(       "un3",       &un3,       "un3/F" );
    t1->Branch(       "un4",       &un4,       "un4/F" );
    t1->Branch(       "un5",       &un5,       "un5/F" );
    t1->Branch(       "un6",       &un6,       "un6/F" );
    t1->Branch(       "nup",       &nup,       "nup/F" );
    t1->Branch(     "ndown",     &ndown,     "ndown/F" );
    t1->Branch(  "cellType",  &cellType,  "cellType/F" );
    t1->Branch(     "event",     &event,     "event/F" );
    t1->Branch(   "simhits",   &simhits,  "simmhits/F" );
    t1->Branch( "rechitsum", &rechitsum, "rechitsum/F" );
    t1->Branch(    "rechit",    &rechit,    "rechit/F" ); // Unsaturated rechit estimation
    t1->Branch( "saturated", &saturated, "saturated/F" ); // Saturated rechit value

    for (int i = 0; i < t->GetEntries(); ++i) {
        t->GetEntry(i);
        /*
        ** Fit results:
        ** mean  122.092 \pm 0.00015
        ** sigma 0.013
        */
        //if(!cellType) rechit = 122.24305*simhits; // 300um cells
        if(!cellType && layer>0 && simhits>0) rechit = 122.092*simhits; // 300um cells
        //else rechit = 182.13023*simhits;          // 200um cells
        t1->Fill();
    }
    fout->cd();
    fout->Write();
    fout->Close();
}
