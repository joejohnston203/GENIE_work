#include <TStyle.h>
#include <TLatex.h>

/*
To create a validation plot, first run the modified version of Nieves' fortran code
qe-gen-integral-dOmega.o, which calculates dSigma/dEl as a function of Tl. (Note that
this will require setting nuclear parameters in fort.4). Use the version that allows 
you to set bounds on values of costheta to integrate over. The output will be stored
in fort.60. Next, get a root file which has gevgen data created with the Nieves' cross
section, preferably with the LFG model. Then run root, load this file (.L validation.C),
and call validation plot with the proper parameters. If you want a legend on the graph,
right click on the edge and select build legend. The range of the y-axis may also need
to be edited if the output from Nieves is larger than that from the root file.
*/

/*
Make a validation plot using one file that contains the output from Nieves' code and another
that uses the results of a gevgen run. The Nieves text file should be generated from the
code in qe-gen-integral-dOmega.f. The gevgen run should use Nieves' cross section.
-rootFile is the name of the gevgen output file, stored as a .root.gst file
-spline is the .root file where the spline is stored that was used when running gevgen
-target is of the form C12, Pb208, etc to describe the target
-Enu is the neutrino beam energy in GeV
-nievesData is the text file containing Tmu in the first column and dSigma/dEl in the second
-title, xlabel, and ylabel are used to format the graph
-colorRoot will set the colors of the line from the root file. The line from the Nieves file
  will automatically be black, which can be changed after the graph is made
-numBins can be used to change how the histogram that is used for the root file is generated
*/

// Example call:
// validation_plot("numu_C12_run1_noFSI.gst.root","data_noR_noC_S_spline.txt","C12",1.0,-0.05,0.05,"/data/jpj13/validation_Summer_2015/validation/C12_dir/C12_Enu_1.0_RPA0_coul0/C12_1.0GeV_ctmin-0.05_ctmax0.05",false)

void validation_plot(double ymax, TString rootFile,TString spline,
			TString target,double Enu,
		    double cosThetaMin = -1.0,double cosThetaMax = 1.0,
		    TString nievesData = "fort.60",
		    bool same = false,
		    TString title = "",
		    TString xlabel = "T_lepton(GeV)",
		    TString ylabel = "Dif XSec dSigma/dEl, 10^-38 cm^2/GeV",
		    int colorRoot = kRed,
		    int numBins = 100)
{
  set_root_env();

  // Y will be T_lepton = El-0.106
  // Tmin will be 0 and Tmax will be El
  xsec_vs_Y(ymax,rootFile,spline,target,Enu,"El-0.106",0.0,Enu,cosThetaMin,cosThetaMax,same,colorRoot,
	    "GENIE Code",xlabel,ylabel,numBins);

  TGraph *nieves = new TGraph(nievesData);
  
  nieves->GetHistogram()->SetMinimum(0.0);
  nieves->GetHistogram()->SetMaximum(ymax);
  nieves->Draw("C");

  // Place nicely formatted title
  gStyle->SetOptTitle(0);
  add_plot_label(title, 0.5, 0.96);
}

/*
Plot the differential cross section dSigma/dY where Y is some variable stored in the root file
-file is the filename of the first root file to be graphed
-type is one of N_RFG, N_LFG, LS_RFG, LS_LFG for file1
-spline is the root file where from the spline is stored that was used when running gevgen, or 
 a text file containing the data as formatted by xmlGetData.sh
-target is of the form C12, Pb208, etc for file1
-Enu is the neutrino beam energy in GeV for file1
-cmin and cmax are bounds on cos(theta), with theta the outoging lepton angle
-color is the color of the line
*/
int xsec_vs_Y(double ymax,TString file,TString splineName,TString target,double Enu,
	      TString Y,double Ymin,double Ymax,
	      double cmin = -1.0,double cmax = 1.0,
	      bool same = false,int color = kBlack,
	      TString title = "",TString xlabel = "",TString ylabel = "",
	      int numBins = 100)
{
  // Get the total cross section at the current energy
  Double_t totalXSec;
  if(splineName.EndsWith(".root")){
    cout << "Spline is a root file" << endl;
    TFile * spline = new TFile(splineName); //Get root spline to get total xsec
    TDirectory * dir1 = (TDirectory*) spline->Get("nu_mu_" + target);
    TGraph * graph1 = (TGraph*) dir1->Get("qel_cc_n");
    totalXSec = graph1->Eval(Enu);
  }else{
    cout << "Spline is a text file. Attempting to create spline from given points." << endl;
    // Assume data is formatted as required to make a TGraph
    TGraph * tempGraph = new TGraph(splineName);
    double factor = 3.90e10; // Splines from xml file are smaller than the root file by factor
    totalXSec = factor*tempGraph->Eval(Enu);
    cout << "Total XSec = " << totalXSec << endl;
  }

  // create histogram binned according to Y
  TChain * chain1 = new TChain("gst");
  chain1->Add(file);

  // Put togther title for curve
  TString title1;
  if(title == ""){
    TString tgtStr,eStr;
    tgtStr = target + " target";
    eStr += Enu;
    eStr.ReplaceAll(" ","");
    eStr = "Enu = " + eStr + " GeV";
    title1 = "dSigma/d(" + Y + ") for " + tgtStr + ", " + eStr;
  }else{
    title1 = title;
  }

  // Parameters for TH1D are name, title, number bins, lower bound, upper bound
  TH1D* hst_Y = new TH1D("hst_Y",title1,numBins,Ymin,Ymax);

  TString cminstr = "", cmaxstr = "";
  cminstr += cmin;
  cmaxstr += cmax;
  if(same){
    chain1->Draw(Y+">>hst_Y","qel&&cthl>="+cminstr+"&&cthl<="+cmaxstr,"same");
  }else{
    chain1->Draw(Y+">>hst_Y","qel&&cthl>="+cminstr+"&&cthl<="+cmaxstr);
  }

  // Scale so the y axis is cross section, not event count
  double binWidth1 = (Ymax-Ymin)/double(numBins);
  double totalEvents = chain1->GetEntries("qel");
  hst_Y->Scale(totalXSec/totalEvents/binWidth1);

  // Format the histogram
  hst_Y->GetYaxis()->SetRangeUser(0.0,ymax);  
  hst_Y->SetLineColor(color);

  TAxis* xax = hst_Y->GetXaxis();
  if(xlabel == "")
    xax->SetTitle(Y);
  else
    xax->SetTitle(xlabel);

  TAxis* yax = hst_Y->GetYaxis();
  if(ylabel == "")
    yax->SetTitle("Differential XSec, dSigma/d(" + Y + ")");
  else
    yax->SetTitle(ylabel);

  return 0;
}

// code to nicely format the plots
void set_root_env(){

  //TStyle* genieStyle = new TStyle("genieStyle", "GENIE Style");

  //set the background color to white
 gStyle->SetFillColor(10);
 gStyle->SetFrameFillColor(10);
 gStyle->SetCanvasColor(10);
 gStyle->SetPadColor(10);
 gStyle->SetTitleFillColor(0);
 gStyle->SetStatColor(10);

//dont put a colored frame around the plots
 gStyle->SetFrameBorderMode(0);
 gStyle->SetCanvasBorderMode(0);
 gStyle->SetPadBorderMode(0);
 gStyle->SetLegendBorderSize(3);

//use the primary color palette
 gStyle->SetPalette(1,0);

//set the default line color for a histogram to be black
 gStyle->SetHistLineColor(kBlack);

//set the default line color for a fit function to be red
 gStyle->SetFuncColor(kRed);

//make the axis labels black
 gStyle->SetLabelColor(kBlack,"xyz");

//set the default title color to be black
 gStyle->SetTitleColor(kBlack);
 
//set the margins
 gStyle->SetPadBottomMargin(0.18);
 gStyle->SetPadTopMargin(0.08);
 gStyle->SetPadRightMargin(0.08);
 gStyle->SetPadLeftMargin(0.17);

//set axis label and title text sizes
 gStyle->SetLabelFont(42,"xyz");
 gStyle->SetLabelSize(0.06,"xyz");
 gStyle->SetLabelOffset(0.015,"xyz");
 gStyle->SetTitleFont(42,"xyz");
 gStyle->SetTitleSize(0.05,"xyz");
 gStyle->SetTitleOffset(1.4,"y");
 gStyle->SetTitleOffset(1.3,"x");
 gStyle->SetStatFont(42);
 gStyle->SetStatFontSize(0.07);
 gStyle->SetTitleBorderSize(1);
 gStyle->SetStatBorderSize(0);
 gStyle->SetTextFont(42);
 gStyle->SetTitleW(0.5);
 gStyle->SetTitleH(0.1);

//set line widths
 gStyle->SetFrameLineWidth(2);
 gStyle->SetFuncWidth(2);
 gStyle->SetHistLineWidth(2);

//set the number of divisions to show
 gStyle->SetNdivisions(506, "xy");
 //gStyle->SetPadTickX(-50202);

//turn off xy grids
 gStyle->SetPadGridX(0);
 gStyle->SetPadGridY(0);

//set the tick mark style
 gStyle->SetPadTickX(1);
 gStyle->SetPadTickY(1);

//turn off stats
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(0);

//marker/line settings
//gStyle->SetMarkerStyle(20);
 gStyle->SetMarkerSize(.95);//0.7
 gStyle->SetLineWidth(2); 
 gStyle->SetErrorX(0);
 gStyle->SetHistLineStyle(0); //It was 3 for a dotted line

//done
 gStyle->cd();
 gROOT->ForceStyle();
 gStyle->ls();

}

// Code that can optionally be used to add a title
void add_plot_label( char* label, double x, double y, double size = 0.05, int color = 1, int font = 62, int align = 22 ){

  TLatex *latex = new TLatex( x, y, label );
  latex->SetNDC();
  latex->SetTextSize(size);
  latex->SetTextColor(color);
  latex->SetTextFont(font);
  latex->SetTextAlign(align);
  latex->Draw();

}
