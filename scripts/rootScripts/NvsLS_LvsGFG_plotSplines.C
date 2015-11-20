// Plots the four splines for LS_GFG,LS_LFG,N_GFG,N_LFG
void plotSplines(TString target,TString process = "qel_cc_n"){
  TString types[4] = {"LS_GFG","LS_LFG","N_GFG","N_LFG"};
  TString files[4];
  for(int i=0;i<4;i++) files[i] = "rootfile_" + types[i] + ".root";
  Color_t colors[4] = {kBlue,kCyan,kRed,kMagenta};

  TGraph * graphs[4];
  TMultiGraph * mg = new TMultiGraph();

  int i = 0;
  TFile * tempFile = new TFile(files[i]);
  TDirectory * tempDir = (TDirectory*) tempFile->Get("nu_mu_" + target);
  graphs[i] = (TGraph*) tempDir->Get("qel_cc_n");
  graphs[i]->SetLineColor(colors[i]);
  graphs[i]->SetTitle(types[i]);
  graphs[i]->Draw("AL");
  c1->SaveAs(types[i] + "_" + target + ".png");
  mg->Add(graphs[i], "AL");

  i = 1;
  TFile * tempFile = new TFile(files[i]);
  TDirectory * tempDir = (TDirectory*) tempFile->Get("nu_mu_" + target);
  graphs[i] = (TGraph*) tempDir->Get("qel_cc_n");
  graphs[i]->SetLineColor(colors[i]);
  graphs[i]->SetTitle(types[i]);
  graphs[i]->Draw("AL");
  c1->SaveAs(types[i] + "_" + target + ".png");
  mg->Add(graphs[i], "AL");

  i = 2;
  TFile * tempFile = new TFile(files[i]);
  TDirectory * tempDir = (TDirectory*) tempFile->Get("nu_mu_" + target);
  graphs[i] = (TGraph*) tempDir->Get("qel_cc_n");
  graphs[i]->SetLineColor(colors[i]);
  graphs[i]->SetTitle(types[i]);
  graphs[i]->Draw("AL");
  c1->SaveAs(types[i] + "_" + target + ".png");
  mg->Add(graphs[i], "AL");

  i = 3;
  TFile * tempFile = new TFile(files[i]);
  TDirectory * tempDir = (TDirectory*) tempFile->Get("nu_mu_" + target);
  graphs[i] = (TGraph*) tempDir->Get("qel_cc_n");
  graphs[i]->SetLineColor(colors[i]);
  graphs[i]->SetTitle(types[i]);
  graphs[i]->Draw("AL");
  c1->SaveAs(types[i] + "_" + target + ".png");
  mg->Add(graphs[i], "AL");
  
  mg->Draw();
  c1->BuildLegend();
  //c1->SetTitle("XSec versus Energy for a " + target + " target, " + process + " events");
}
