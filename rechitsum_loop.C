#include "rechitsum_new.C"

void rechitsum_loop(){
    gStyle->SetOptStat(0);
    gStyle->SetOptFit();

    TCanvas* c1 = new TCanvas("c1","c1",1000,800);
    TCanvas* c2 = new TCanvas("c2","c2",1000,800);
    /*TCanvas* c3 = new TCanvas("c3","c3",1000,800);
    TCanvas* c4 = new TCanvas("c4","c4",1000,800);
    TCanvas* c5 = new TCanvas("c5","c5",1000,800);*/
    c1->Divide(2,2);
    c2->Divide(2,2);
    /*c3->Divide(2,2);
    c4->Divide(2,2);
    c5->Divide(2,2);*/

    /*double energies[]  = {  5, 10, 15, 20, 30, 40, 60, 80,100,140,200,280,400,550,750,1000,1400,2000,2800};
    int bins[]          = {100,100,100,100,100,100,100,100,100,100,100,100,100,100,100, 100, 100, 100, 100};
    int range[]         = {  3,  5,  5, 10, 10, 10, 15, 15, 15, 20, 25, 30, 40, 50,100, 150, 400, 800,1100};
    double fit_cut[]    = {1.6,1.6,1.6,1.6,1.6,1.6,1.6,1.6,1.6,1.6,1.6,1.6,1.1,1.0,1.2, 1.3, 1.3, 1.3, 1.3};

    double scemean[]   = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double scemeane[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double sceres[]    = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double scerese[]   = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double energiese[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};*/

    /*double energies[]  = {750,1000,1400,2000,2800};
    int bins[]         = {100, 100, 100, 100, 200};
    int range[]        = {100, 100, 100, 100, 200};
    double fit_cut[]   = {1.0, 0.7, 0.7, 0.9, 0.9};*/

    double energies[]  = {750,1000,1400,2000,2800};
    int bins[]         = {100, 100, 100, 100, 100};
    int range[]        = {100, 150, 200, 300, 400};
    double fit_cut[]   = {1.0, 1.0, 1.0, 1.0, 1.0};

    double scemean[]   = {0.,0.,0.,0.,0.};
    double scemeane[]  = {0.,0.,0.,0.,0.};
    double sceres[]    = {0.,0.,0.,0.,0.};
    double scerese[]   = {0.,0.,0.,0.,0.};
    double energiese[] = {0.,0.,0.,0.,0.};

    for(int j = 0; j < 5; ++j){
        if(j < 4) c1->cd(j%4+1);
        else if(j < 8) c2->cd(j%4+1);
        //else if(j < 12) c3->cd(j%4+1);
        //else if(j < 16) c4->cd(j%4+1);
        //else c5->cd(j%4+1);
        std::vector<double> temp = rechitsum_new(energies[j],bins[j],range[j],fit_cut[j]);
        scemean[j]  = temp[1];
        scemeane[j] = temp[4];
        //sceres[j]   = temp[2]/energies[j];
        sceres[j]   = temp[2]/temp[1];
        scerese[j]  = sceres[j]*sqrt(pow(temp[5]/temp[2],2)+pow(temp[4]/temp[1],2));
        std::cout   << " fit results for " << energies[j]
        << " mean " << scemean[j] <<  "+-" << scemeane[j]
        <<  " res " <<  sceres[j] <<  "+-" <<  scerese[j] << std::endl;
    }

    TCanvas* c_res = new TCanvas("c_res","c_res",1);
    TGraphErrors *gr = new TGraphErrors(5,energies,sceres,energiese,scerese);
    gr->SetTitle("gamma resolution versus energy;E [GeV];width/mean");
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    TF1  *f2 = new TF1("f2","sqrt(([0]/sqrt(x))**2+([1]/x)**2+([2])**2)");
    f2->SetNpx(1000);
    gr->Fit("f2");
    gr->Draw("AP");

    /*TString outname = "outputFiles/out0"+to_string(df)+".root";
    TFile* out = new TFile(outname,"RECREATE");
    gr->Write();
    out->Close();

    //Print PDFs
    TString cname1  = "outputFiles2/canvas1_df0"+to_string(df)+"_LSaver.pdf";
    TString cname2  = "outputFiles2/canvas2_df0"+to_string(df)+"_LSaver.pdf";
    TString cname3  = "outputFiles2/canvas3_df0"+to_string(df)+"_LSaver.pdf";
    TString cname4  = "outputFiles2/canvas4_df0"+to_string(df)+"_LSaver.pdf";
    TString cname5  = "outputFiles2/canvas5_df0"+to_string(df)+"_LSaver.pdf";
    TString cnamegr = "outputFiles2/resplot_df0"+to_string(df)+"_LSaver.pdf";
    c1->Print(cname1);
    c2->Print(cname2);
    c3->Print(cname3);
    c4->Print(cname4);
    //c5->Print(cname5);
    c_res->Print(cnamegr);*/
}
