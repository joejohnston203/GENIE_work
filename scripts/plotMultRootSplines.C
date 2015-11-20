// Root script to plot the graph of a spline from a .root file
// file is .root file
// target is in form C12, Pb208, etc
// process is of form qel_cc_n. Possible processes can be shown in
// root by calling .ls when the root file is open
void plotSpline(TString files[4], int numFiles,
		TString target,
		TString process = "qel_cc_n",
		TString title = "Cross Section vs Energy"){

  TGraph * graphs[numFiles];
  TMultiGraph * mg = new TMultiGraph();

  Color_t colors[1] = {kBlue};

  for(int i=0; i<numFiles; i++){
    TFile * tempFile = new TFile(files[i]);
    TDirectory * tempDir = (TDirectory*) tempFile->Get("nu_mu_" + target);
    graphs[i] = (TGraph*) tempDir->Get(process);
    graphs[i]->SetLineColor(colors[i]);
    graphs[i]->SetTitle(files[i]);
    graphs[i]->Draw("AL");
    c1->SaveAs(target + "_graph_" + i + ".png");
    mg->Add(graph, "AL");
  }
  mg->Draw();
  c1->BuildLegend();
}
