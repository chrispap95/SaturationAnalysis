{
    TCanvas* c = new TCanvas("c","c",1650,600);
    gStyle->SetOptStat(0);
    gPad->SetLeftMargin(0.04);
    Float_t n1 = 9.545;
    Float_t n2 = 9.456;
    Float_t n3 = 10.68;
    Float_t n4 = 8.675;
    Float_t n5 = 8.723;
    Float_t n6 = 10.68;
    Float_t dn1 = 18;
    Float_t dn2 = 17.97;
    Float_t dn3 = 6.939;
    Float_t dn4 = 2.563;
    Float_t dn5 = 2.613;
    Float_t dn6 = 7.066;
    Float_t un1 = 3.122;
    Float_t un2 = 3.096;
    Float_t un3 = 7.742;
    Float_t un4 = 17.42;
    Float_t un5 = 17.77;
    Float_t un6 = 7.838;
    Float_t down = 17.99;
    Float_t up = 19.38;
    //Float_t saturated = 27.74;


    // Make histo
    TH2Poly* h_honeycomb = new TH2Poly();
    h_honeycomb->Honeycomb(0,0,1,13,5);
    h_honeycomb->SetTitle("");
    h_honeycomb->GetXaxis()->SetLabelOffset(999);
    h_honeycomb->GetXaxis()->SetLabelSize(0);
    h_honeycomb->GetXaxis()->SetTickLength(0);
    h_honeycomb->GetYaxis()->SetLabelOffset(999);
    h_honeycomb->GetYaxis()->SetLabelSize(0);
    h_honeycomb->GetYaxis()->SetTickLength(0);

    float a = 0;
    h_honeycomb->Fill(3.5+a,2.5,dn4);
    h_honeycomb->Fill(2.5+a,4,dn3);
    h_honeycomb->Fill(3.5+a,5.5,dn1);
    h_honeycomb->Fill(5+a,5.5,dn2);
    h_honeycomb->Fill(6+a,4,dn6);
    h_honeycomb->Fill(5+a,2.5,dn5);
    h_honeycomb->Fill(5+a,4,down);
    a = 7;
    h_honeycomb->Fill(3.5+a,2.5,n4);
    h_honeycomb->Fill(2.5+a,4,n3);
    h_honeycomb->Fill(3.5+a,5.5,n1);
    h_honeycomb->Fill(5+a,5.5,n2);
    h_honeycomb->Fill(6+a,4,n6);
    h_honeycomb->Fill(5+a,2.5,n5);
    //h_honeycomb->Fill(5+a,4,saturated);
    a = 14;
    h_honeycomb->Fill(3.5+a,2.5,un4);
    h_honeycomb->Fill(2.5+a,4,un3);
    h_honeycomb->Fill(3.5+a,5.5,un1);
    h_honeycomb->Fill(5+a,5.5,un2);
    h_honeycomb->Fill(6+a,4,un6);
    h_honeycomb->Fill(5+a,2.5,un5);
    h_honeycomb->Fill(5+a,4,up);

    h_honeycomb->GetZaxis()->SetTitle("recHit [GeV]");
    h_honeycomb->GetZaxis()->SetTitleOffset(0.54);
    h_honeycomb->GetZaxis()->SetTitleSize(0.04);
    h_honeycomb->SetMarkerSize(1.5);
    h_honeycomb->Draw("colz 0 text");

    // Add necessary text
    TLatex ltx;
    ltx.DrawLatex(3,7,"Previous layer");
    ltx.DrawLatex(9.5,7,"Saturated cell layer");
    ltx.DrawLatex(17.25,7,"Next layer");

    // Add separating lines
    TLine l1(7.8,0,7.8,8);
    l1.Draw();
    TLine l2(14.7,0,14.7,8);
    l2.Draw();
    //TMarker mr1(4.35,4,3);
    //mr1.Draw();
    TMarker mr2(11.3,4,3);
    mr2.Draw();
    //TMarker mr3(18.23,4,3);
    //mr3.Draw();
}
