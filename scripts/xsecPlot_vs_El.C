
/*
Plot the differential cross section dSigma/dEl vs Tl, the kinetic energy of the lepton
-file is the filename of the first root file to be graphed
-type is one of N_RFG, N_LFG, LS_RFG, LS_LFG for file1
-target is of the form C12, Pb208, etc for file1
-Enu is the neutrino beam energy in GeV for file1
-line1 is the color of the line
*/
int xsec_vs_El(TString file,TString type,
	       TString splineName1,
	       TString target,double Enu,
	    double Elmin,double Elmax,
	    bool same = false,int color = kBlack,int numBins = 100,
	       TString title = "")
{
  // Type must be one of the four options or else there won't be a spline
  if(!(type=="LS_RFG"||type=="LS_LFG"||type=="N_RFG"||type=="N_LFG")){
    std::cout << "Invalid option for type" << std::endl;
    return 0;
  }

  //TString splineName1 = "/usr/GENIEdevel/src/LlewellynSmith/splines/rootfile_" + type + ".root";
  TFile * spline = new TFile(splineName1); //Get root spline to get total xsec
  TDirectory * dir1 = (TDirectory*) spline->Get("nu_mu_" + target);
  TGraph * graph1 = (TGraph*) dir1->Get("qel_cc_n");
  Double_t totalXSec1 = graph1->Eval(Enu);

  // create histogram of Q2 data
  TChain * chain1 = new TChain("gst");
  chain1->Add(file);

  // Put togther title for curve
  TString title1;
  if(title == ""){
    TString xsecStr,fgmStr,tgtStr,eStr;
    if(type == "LS_RFG"){
      xsecStr = "No RPA";
      fgmStr = "Relativistc FG";
    }else if(type == "N_RFG"){
      xsecStr = "RPA";
      fgmStr = "Relativistc FG";
    }else if(type == "LS_LFG"){
      xsecStr = "No RPA";
      fgmStr = "Local FG";
    }else if(type == "N_LFG"){
      xsecStr = "RPA";
      fgmStr = "Local FG";
    }
    tgtStr = target + " tgt";
    eStr += Enu;
    eStr.ReplaceAll(" ","");
    eStr = "Enu = " + eStr + " GeV";
    title1 = xsecStr + ", " + fgmStr;
    //title1 = xsecStr + ", " + fgmStr + ", " + tgtStr + ", " + eStr;
  }else{
    title1 = title;
  }

  // Parameters for TH1D are name, title, number bins, lower bound, upper bound
  TH1D* hst_El = new TH1D("hst_El",title1,numBins,Elmin,Elmax);
  if(same){
    double totalEvents1 = (double)(chain1->Draw("El-0.106>>hst_El","qel","same"));
  }else{
   double totalEvents1 = (double)(chain1->Draw("El-0.106>>hst_El","qel"));
  }

  double binWidth1 = (Elmax-Elmin)/double(numBins);
  // Scale so the y axis is cross section, not event count
  hst_El->Scale(totalXSec1/totalEvents1/binWidth1);

  // Format the histogram
  hst_El->SetLineColor(color);
  TAxis* xax = hst_El->GetXaxis();
  xax->SetTitle("Tmu (GeV)");
  TAxis* yax = hst_El->GetYaxis();
  yax->SetTitle("Differential XSec, dSigma/dEl");

  return 0;
}


/*
Plot the differential XSec dSigma/dQ2
-file is the filename of the first root file to be graphed
-type is one of N_RFG, N_LFG, LS_RFG, LS_LFG for file1
-target is of the form C12, Pb208, etc for file1
-Enu is the neutrino beam energy in GeV for file1
-line1 is the color of the line
*/
int xsec_vs_q2(TString file,TString type, TString splineName,
	       TString target,double Enu,
	       double q2min,double q2max,
	       bool same = false,int color = kBlack,int numBins = 100,
	       TString title = "")
{
  // Type must be one of the four options or else there won't be a spline
  if(!(type=="LS_RFG"||type=="LS_LFG"||type=="N_RFG"||type=="N_LFG")){
    std::cout << "Invalid option for type" << std::endl;
    return 0;
  }

  //TString splineName1 = "/usr/GENIEdevel/src/LlewellynSmith/splines/rootfile_" + type + ".root";
  TFile * spline = new TFile(splineName); //Get root spline to get total xsec
  TDirectory * dir1 = (TDirectory*) spline->Get("nu_mu_" + target);
  TGraph * graph1 = (TGraph*) dir1->Get("qel_cc_n");
  Double_t totalXSec1 = graph1->Eval(Enu);

  // create histogram of Q2 data
  TChain * chain1 = new TChain("gst");
  chain1->Add(file);

  // Put togther title for curve
  TString title1;
  if(title == ""){
    TString xsecStr,fgmStr,tgtStr,eStr;
    if(type == "LS_RFG"){
      xsecStr = "No RPA";
      fgmStr = "Relativistc FG";
    }else if(type == "N_RFG"){
      xsecStr = "RPA";
      fgmStr = "Relativistc FG";
    }else if(type == "LS_LFG"){
      xsecStr = "No RPA";
      fgmStr = "Local FG";
    }else if(type == "N_LFG"){
      xsecStr = "RPA";
      fgmStr = "Local FG";
    }
    tgtStr = target + " tgt";
    eStr += Enu;
    eStr.ReplaceAll(" ","");
    eStr = "Enu = " + eStr + " GeV";
    TString title1 = xsecStr + ", " + fgmStr;
    //TString title1 = xsecStr + ", " + fgmStr + ", " + tgtStr + ", " + eStr;
  }else{
    title1 = title;
  }

  // Parameters for TH1D are name, title, number bins, lower bound, upper bound
  TH1D* hst_Q2 = new TH1D("hst_Q2",title1,numBins,q2min,q2max);
  if(same){
    double totalEvents1 = (double)(chain1->Draw("Q2>>hst_Q2","qel","same"));
  }else{
    double totalEvents1 = (double)(chain1->Draw("Q2>>hst_Q2","qel"));
  }
  double binWidth1 = (q2max-q2min)/double(numBins);
  // Scale so the y axis is cross section, not event count
  hst_Q2->Scale(totalXSec1/totalEvents1/binWidth1);
  
  // Format the histogram
  hst_Q2->SetLineColor(color);
  TAxis* xax = hst_Q2->GetXaxis();
  xax->SetTitle("Q2 (GeV)");
  TAxis* yax = hst_Q2->GetYaxis();
  yax->SetTitle("Differential XSec, dSigma/dQ2");
  
  return 0;
}
