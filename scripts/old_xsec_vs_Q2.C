/*
-file1 is the filename of the first root file to be graphed
-type1 is one of N_RFG, N_LFG, LS_RFG, LS_LFG for file 1
-target1 is of the form C12, Pb208, etc for file 1
-E1 is the neutrino beam energy in GeV for file 1
-file2 is the filename of the second root file to be graph
-type2, target2, and E2 are for the second file
-line1,line2 are the integers corresponding to the colors of the lines
*/

int q2Graph(TString file1,TString type1,TString target1,double E1,
	    TString file2,TString type2,TString target2,double E2,
	    double q2min,double q2max,
	    int color1=kRed,int color2=kBlue,int numBins=100)
{
  TString splineName1 = "/usr/GENIEdevel/src/LlewellynSmith/splines/rootfile_" + type1 + ".root";
  TFile * spline = new TFile(splineName1); //Get root spline to get total xsec
  TDirectory * dir1 = (TDirectory*) spline->Get("nu_mu_" + target1);
  TGraph * graph1 = (TGraph*) dir1->Get("qel_cc_n");
  Double_t totalXSec1 = graph1->Eval(E1);

  // create histogram of Q2 data
  TChain * chain1 = new TChain("gst");
  chain1->Add(file1);

  // Parameters for TH1D are name, title, number bins, lower bound, upper bound
  TString title1 = type1 + " XSec, " + target1 + " target";
  TH1D* hst_Q2_1 = new TH1D("hst_Q2_1",title1,numBins,q2min,q2max);
  double totalEvents1 = (double)(chain1->Draw("Q2>>hst_Q2_1","qel"));
  double binWidth1 = (q2max-q2min)/double(numBins);
  // Scale so the y axis is cross section, not event count
  hst_Q2_1->Scale(totalXSec1/totalEvents1/binWidth1);
  hst_Q2_1->SetLineColor(color1);

  //Create Graph for second file
  TString splineName2 = "/usr/GENIEdevel/src/LlewellynSmith/splines/rootfile_" + type2 + ".root";
  TFile * spline = new TFile(splineName2); //Get root spline to get total xsec
  TDirectory * dir2 = (TDirectory*) spline->Get("nu_mu_" + target2);
  TGraph * graph2 = (TGraph*) dir2->Get("qel_cc_n");
  Double_t totalXSec2 = graph2->Eval(E2);

  // create histogram of Q2 data
  TChain * chain2 = new TChain("gst");
  chain2->Add(file2);

  // Parameters for TH1D are name, title, number bins, lower bound, upper bound
  TString title2 = type2 + " XSec, " + target2 + " target";
  TH1D* hst_Q2_2 = new TH1D("hst_Q2_2",title2,numBins,q2min,q2max);
  double totalEvents2 = (double)(chain2->Draw("Q2>>hst_Q2_2","qel","same"));
  double binWidth2 = (q2max-q2min)/double(numBins);
  // Scale so the y axis is cross section, not event count
  hst_Q2_2->Scale(totalXSec2/totalEvents2/binWidth2);
  hst_Q2_2->SetLineColor(color2);

  return 0;
}

/*
TString file1 = "N_RFG_numu_C12.gst.root";
TString type1 = "N_RFG";
TString target1 = "C12";
double E1 = 1.0;

TString file2 = "LS_RFG_numu_C12.gst.root";
TString type2 = "LS_RFG";
TString target2 = "C12";
double E2 = 1.0;
double q2min = 0.0;
double q2max = 2.0;
*/
