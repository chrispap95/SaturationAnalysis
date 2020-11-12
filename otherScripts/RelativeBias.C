{
  TFile* file = TFile::Open("~/profile.root");
  TCanvas* c1 = (TCanvas*)file->Get("canvas1");
  TProfile* prf = (TProfile*)c1->GetPrimitive("MVA_DNN_CPUtest_reg_tgt0_rtgt0_pfx");
  int nbins = prf->GetNbinsX();
  int bin0 = prf->GetXaxis()->GetBinUpEdge(0);
  int bin50 = prf->GetXaxis()->GetBinUpEdge(nbins);

  TH1F* href = new TH1F("href","href",nbins,bin0,bin50);
  for (int i = 1; i <= nbins; ++i) {
    href->SetBinContent(i,prf->GetBinContent(i)/prf->GetBinCenter(i));
    href->SetBinError(i,prf->GetBinError(i)/prf->GetBinCenter(i));
  }
  TCanvas c2;
  href->Draw();
}
