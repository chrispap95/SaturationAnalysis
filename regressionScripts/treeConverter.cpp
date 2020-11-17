/*
** Script:
**    - Renames branch names
**    - Adds rechit truth estimation
**    - Skims the data for cells that don't qualify for training set
**        ~ cells that are on the wafer edges
*/

void treeConverter(TString input, bool isTraining){
    TFile* f = TFile::Open(input+".root");
    TTree* t = dynamic_cast< TTree* >(f->Get("t1"));
    Float_t n1, n2, n3, n4, n5, n6, nup, ndown, rechit, layer, rechitsum;
    Float_t un1, un2, un3, un4, un5, un6, dn1, dn2, dn3, dn4, dn5, dn6;
    Float_t saturated, simhits, thickness, event, cellU, cellV;
    TString prefix = "ML"; // or for older code prefix = "ML"
    t->SetBranchAddress( prefix+"layer", &layer );
    t->SetBranchAddress( prefix+"cellU", &cellU );
    t->SetBranchAddress( prefix+"cellV", &cellV );
    t->SetBranchAddress( prefix+"n1", &n1 );
    t->SetBranchAddress( prefix+"n2", &n2 );
    t->SetBranchAddress( prefix+"n3", &n3 );
    t->SetBranchAddress( prefix+"n4", &n4 );
    t->SetBranchAddress( prefix+"n5", &n5 );
    t->SetBranchAddress( prefix+"n6", &n6 );
    t->SetBranchAddress( prefix+"un1", &un1 );
    t->SetBranchAddress( prefix+"un2", &un2 );
    t->SetBranchAddress( prefix+"un3", &un3 );
    t->SetBranchAddress( prefix+"un4", &un4 );
    t->SetBranchAddress( prefix+"un5", &un5 );
    t->SetBranchAddress( prefix+"un6", &un6 );
    t->SetBranchAddress( prefix+"dn1", &dn1 );
    t->SetBranchAddress( prefix+"dn2", &dn2 );
    t->SetBranchAddress( prefix+"dn3", &dn3 );
    t->SetBranchAddress( prefix+"dn4", &dn4 );
    t->SetBranchAddress( prefix+"dn5", &dn5 );
    t->SetBranchAddress( prefix+"dn6", &dn6 );
    t->SetBranchAddress( prefix+"nup", &nup );
    t->SetBranchAddress( prefix+"ndown", &ndown );
    t->SetBranchAddress( prefix+"thickness", &thickness );
    t->SetBranchAddress( prefix+"event", &event );
    t->SetBranchAddress( prefix+"simHits", &simhits );
    t->SetBranchAddress( prefix+"rechitsum", &rechitsum );
    t->SetBranchAddress( prefix+"saturated", &saturated );

    TFile* fout = new TFile(input+"_converted.root","RECREATE");
    TTree* tree = new TTree("tree","sample");
    tree->Branch( "layer", &layer, "layer/F" );
    tree->Branch( "n1", &n1, "n1/F" );
    tree->Branch( "n2", &n2, "n2/F" );
    tree->Branch( "n3", &n3, "n3/F" );
    tree->Branch( "n4", &n4, "n4/F" );
    tree->Branch( "n5", &n5, "n5/F" );
    tree->Branch( "n6", &n6, "n6/F" );
    tree->Branch( "dn1", &dn1, "dn1/F" );
    tree->Branch( "dn2", &dn2, "dn2/F" );
    tree->Branch( "dn3", &dn3, "dn3/F" );
    tree->Branch( "dn4", &dn4, "dn4/F" );
    tree->Branch( "dn5", &dn5, "dn5/F" );
    tree->Branch( "dn6", &dn6, "dn6/F" );
    tree->Branch( "un1", &un1, "un1/F" );
    tree->Branch( "un2", &un2, "un2/F" );
    tree->Branch( "un3", &un3, "un3/F" );
    tree->Branch( "un4", &un4, "un4/F" );
    tree->Branch( "un5", &un5, "un5/F" );
    tree->Branch( "un6", &un6, "un6/F" );
    tree->Branch( "nup", &nup, "nup/F" );
    tree->Branch( "ndown", &ndown, "ndown/F" );
    tree->Branch( "thickness", &thickness, "thickness/F" );
    tree->Branch( "event", &event, "event/F" );
    tree->Branch( "simhits", &simhits, "simmhits/F" );
    tree->Branch( "rechitsum", &rechitsum, "rechitsum/F" );
    tree->Branch( "rechit", &rechit, "rechit/F" ); // Unsaturated rechit estimation
    tree->Branch( "saturated", &saturated, "saturated/F" ); // Saturated rechit value

    for (int i = 0; i < t->GetEntries(); ++i) {
        t->GetEntry(i);
        /*
        ** Fit results:
        ** mean  122.092 \pm 0.00015
        ** sigma 0.013
        */
        if (thickness == 300) rechit = 122.092*simhits; // 300um cells
        else rechit = -1;

        bool isHalfCell = (cellV==0 || cellV==15 || cellU==0 || cellU==15) ? 1 : 0;
        bool hasSatNeighbor = (n1>27.7 || n2>27.7 || n3>27.7 || n4>27.7 || n5>27.7 || n6>27.7) ? 1 : 0;
        if (!isTraining) tree->Fill();
        else if (isHalfCell || hasSatNeighbor || thickness!=300 || layer<=0 || layer!=12) continue;
        else tree->Fill();
    }
    fout->cd();
    fout->Write();
    fout->Close();
}
