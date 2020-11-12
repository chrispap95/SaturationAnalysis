#include "RootStyle.cc"

void resComparison(){
    set_root_style();
    //set_tdr_style();
    //set_vasu_style();
    //setTDRStyle();

    TCanvas* c = new TCanvas("c","canvas",1);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    TFile* f0 = TFile::Open("outputFiles/out_truth.root");
    TFile* f1 = TFile::Open("outputFiles/out_mlCorr.root");
    TFile* f2 = TFile::Open("outputFiles/out_noCorr.root");

    TGraphErrors* gr0 = (TGraphErrors*)f0->Get("Graph");
    TGraphErrors* gr1 = (TGraphErrors*)f1->Get("Graph");
    TGraphErrors* gr2 = (TGraphErrors*)f2->Get("Graph");
    TMultiGraph* g = new TMultiGraph();

    gr0->SetTitle("truth estimation");
    gr1->SetTitle("ML correction");
    gr2->SetTitle("no correction");

    gr0->SetMarkerColor(kBlack);
    gr1->SetMarkerColor(kGreen);
    gr2->SetMarkerColor(kBlue);
    gr0->SetLineColor(kBlack);
    gr1->SetLineColor(kGreen);
    gr2->SetLineColor(kBlue);

    gr0->GetFunction("f2")->SetLineColor(kBlack);
    /*
    gr1->GetFunction("f2")->SetLineColor(kGreen);
    gr2->GetFunction("f2")->SetLineColor(kBlue);
    */

    g->Add(gr0,"P");
    g->Add(gr1,"P");
    g->Add(gr2,"P");
    g->SetTitle(";E [GeV];#frac{width}{mean}");
    g->Draw("AP");

    TLegend* lg = new TLegend(0.5,0.5,0.88,0.88);
    lg->SetBorderSize(0);
    lg->AddEntry(gr0,"truth estimation","lp");
    lg->AddEntry(gr1,"ML correction","lp");
    lg->AddEntry(gr2,"no correction","lp");
    lg->Draw();
    gPad->Update();
    //c->Print("presentation/deadfrac_comparison05_noadj.eps")
}
