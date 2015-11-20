// Root script to plot the graph of a spline from a .root file
// file is .root file
// target is in form C12, Pb208, etc
// process is of form qel_cc_n. Possible processes can be shown in
// root by calling .ls when the root file is open
void plotSpline(TString file, TString target,
		TString process = "qel_cc_n",
		TString title = "Cross Section vs Energy"){

  TGraph * graph;
  TMultiGraph * mg = new TMultiGraph();

  Color_t color = kBlue;

  TFile * tempFile = new TFile(file);
  TDirectory * tempDir = (TDirectory*) tempFile->Get("nu_mu_" + target);
  graph = (TGraph*) tempDir->Get(process);
  graph->SetLineColor(color);
  graph->SetTitle(title);
  graph->Draw("A*");
  c1->SaveAs(target + ".png");
  //mg->Add(graph, "AL");

  //mg->Draw();
  //c1->BuildLegend();
}
