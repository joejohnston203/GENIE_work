#include <TStyle.h>
#include <TLatex.h>

/*
To create a validation plot, check the dsigmadTl_C12 and dsigmadTl_Pb208
folders to ensure that run_fortran_C12.sh and run_fortran_Pb208.sh have
been run. These will run Nieves fortran code for different energies and
angles, and store the output so the methods here can access them.

Next, use GENIE to generate .xml files with the desired splines. Put the
information for just the CCQE spline in a file, then run the following
bash commands to extract just the energy and cross section into a text file:
*/
  // sed -n 's/.*nknots="\([0-9]*\).*/#\1/p' spline.xml > spline.txt
  // echo -e "#E\txsec" >> spline.txt
  // sed -n 's/\t<knot> <E> *\([0-9.]*\) <\/E> <xsec> *\([0-9.e-]*\).*/\1\t\2/p' spline.xml  >> spline.txt
/*

The use GENIE to generate .gst.root files for the CCQE xsec model and 
settings you are interested in comparing to Nieves' code.

Finally, update the settings in dsigmadTl_make_plots() to link to your spline
and .gst.root files, and execute this script to generate the plots
*/

void dsigmadTl_make_plots() {
  // Set paramters for plotting here
  // do not inclue .png extension in output_file_name
  TString output_file_name = "N_LFG_200MeV_C12_RPA_Coul";
  // location of .gst.root file from gevgen
  TString rootFile = "~/GENIE_work/runs/2016_05_21/2016_05_21_C12_0.2GeV_N_RPA_Coul_LFG/numu_C12_CCQE_run0.gst.root";
  // 2 column text file (as above) or .root spline file
  TString spline = "~/GENIE_work/runs/2016_05_19_C12_N_RPA_Coul_1GeV/data_C12_N_LFG.xml";
  //TString spline = "/home/joe/GENIE_work/runs/lfgcomp_2016_03_17/2016_03_17_C12_LS_LFG_1GeV/data_2016_03_17_C12_LS_LFG.txt";
  //TString spline = "~/GENIE_work/runs/my_af_runs/2016_05_02_af2_C12_LS_RFG_1GeV/data_2016_05_02_af2_C12_LS_RFG_1GeV.txt";
  TString target = "C12"; // C12 or Pb208
  TString Enu = "0.2"; // in GeV- 0.2, 1.0, or 5.0 
  TString RPA = "1";   // 0 or 1
  TString Coul = "1";  // 0 or 1
  double ymax1 = 25.0; // y axis max for plot where all solid angles are integrated
  double ymax2 = 2.0; // y axis max for plot with specific ctl values
  // "" gives a default title
  //TString plotTitle = "C12, Enu = 1.0 GeV, Nieves noRPA, RFG"; 
  TString plotTitle = output_file_name;
  make_plots(output_file_name, rootFile, spline, target,
	     Enu, RPA, Coul, ymax1, ymax2,plotTitle);
  exit();
}

void make_plots(TString graph_title, TString rootFile,
		TString spline, TString target, TString Enu,
		TString RPA, TString Coul, double ymax1, double ymax2,
		TString plotTitle = "") {
  // Make a plot of dSigma/dTl, where all solid angles have been integrated
  TString ctlminmaxList = "-1.0 1.0";
  save_plot(graph_title + "_total.png",ymax1,rootFile,spline,
	    target,Enu,RPA,Coul,ctlminmaxList,plotTitle);
  
  // Make dSigma/dTl plots for specific cos(theta) angles
  //  TString ctlminmax2= "-1.0 -0.9 -0.05 0.05 0.4 0.5 0.7 0.8 0.9 0.95";
  TString ctlminmax2= "-1.0 -0.9 0.4 0.5 0.9 0.95";
  save_plot(graph_title + ".png",ymax2,rootFile,spline,target,
	    Enu,RPA,Coul,ctlminmax2,plotTitle);
}

/*
Creates and saves a plot for the given settings
  graph_title - name of the created file (without .png extension)
  ymax        - maximum y axis value
  rootFile    - .gst.root file from a gevgen run
  spline      - text file with Energy and XSec in two columns, or a .root
                file created with gspl2root
  target      - C12 or Pb208
  Enu         - String with neutrino energy in GeV, 0.2, 1.0, or 5.0
  RPA         - 0 (off) or 1 (on)
  Coul        - 0 (off) or 1 (on)
  ctlList     - TString with the desired ctl bounds separated by spaces.
                Integration bounds come in pairs, so for example, giving
		"-1.0 1.0" as the argument would create a graph where all
		ctl values are integrated over.
  plotTitle   - plot title
 */
void save_plot(TString graph_title, double ymax, TString rootFile,
	       TString spline, TString target, TString Enu,
	       TString RPA, TString Coul, TString ctlList,
	       TString plotTitle = ""){
  if(plotTitle.EqualTo(""))
    plotTitle = "dSigma/dTl, " + target + 
      ", Enu = " + Enu + " GeV, gevgen vs Nieves' code";
  
  TCanvas * c1 = new TCanvas("c1","c1",900,700);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.15);
  
  // Add a legend
  TLegend* leg1 = new TLegend(.2,.7,.5,.9,"");
  // Get the angles
  TObjArray * ctlminmax = ctlList.Tokenize(TString(" "));
  TString ctlmin, ctlmax;
  TH1D* hst;
  int colors[] = {kRed,kBlue,kGreen,kViolet,kCyan,kYellow,kPink,kBlack};// size 8
  int colorsSize = 8;
  for(int i=0;i+1<ctlminmax->GetEntries();i+=2){
    ctlmin = ((TObjString*)ctlminmax->At(i))->String();
    ctlmax = ((TObjString*)ctlminmax->At(i+1))->String();
    TString nievesName = target + "_" + Enu + "_R" + RPA + "_C" + Coul ;
    TString nievesData = "dsigmadTl_" + target + "/" 
      + nievesName + "/" + nievesName + "_" + ctlmin + "_" + ctlmax + ".txt";

    int color;
    if(i/2<colorsSize) color = colors[i/2];
    else               color = kBlack;

    if(i==0)
      hst = validation_plot(ymax,rootFile,spline,target,Enu,
			    ctlmin,ctlmax,nievesData,false,
			    "T_lepton(GeV)","dSigma/dEl, 10^-38 cm^2/GeV",
			    color);
    else
      hst = validation_plot(ymax,rootFile,spline,target,Enu,
			    ctlmin,ctlmax,nievesData,true,
			    "T_lepton(GeV)","dSigma/dEl, 10^-38 cm^2/GeV",
			    color);
    leg1->AddEntry(hst,"ctl=" + ctlmin + " to " + ctlmax,"L");
  }
  leg1->Draw();

  // Add a title
  gStyle->SetOptTitle(0);
  add_plot_label(plotTitle, 0.5,0.96);
  c1->SaveAs(graph_title);
}

/*
This makes a validation plot using one file that contains the output from Nieves' code 
and another
that uses the results of a gevgen run. The Nieves text file should be generated from the
code in qe-gen-integral-dOmega.f. The gevgen run should use Nieves' cross section.
-ymax is the maximum y axis value
-rootFile is the name of the gevgen output file, stored as a .root.gst file
-spline should be a text file with the energy and xsec extracted from a .xml spline file
 (can also be a root file created with gspl2root)
-target is of the form C12, Pb208, etc to describe the target
-Enu is the neutrino beam energy in GeV
-nievesData is the text file containing Tmu in the first column and dSigma/dEl in the second
-title, xlabel, and ylabel are used to format the graph
-colorRoot will set the colors of the line from the root file. The line from the Nieves file
  will automatically be black, which can be changed after the graph is made
-numBins can be used to change how the histogram that is used for the root file is generated

The TH1D* pointer is returned so the curve can be added to a legend
*/

// Example call:
// validation_plot("numu_C12_run1_noFSI.gst.root","data_noR_noC_S_spline.txt","C12",1.0,-0.05,0.05,"/data/jpj13/validation_Summer_2015/validation/C12_dir/C12_Enu_1.0_RPA0_coul0/C12_1.0GeV_ctmin-0.05_ctmax0.05",false)

TH1D* validation_plot(double ymax, TString rootFile,TString spline,
			TString target,TString Enu,
		    TString cosThetaMin = "-1.0",TString cosThetaMax = "1.0",
		    TString nievesData = "fort.60",
		    bool same = false,
		    TString xlabel = "T_lepton(GeV)",
		    TString ylabel = "dSigma/dEl, 10^-38 cm^2/GeV",
		    int colorRoot = kRed,
		    int numBins = 100)
{
  set_root_env();

  // Y will be T_lepton = El-0.106
  // Tmin will be 0 and Tmax will be El
  TH1D* hst_Y = xsec_vs_Y(ymax,rootFile,spline,target,Enu,"El-0.106",0.0,Enu.Atof(),
			  cosThetaMin,cosThetaMax,same,colorRoot,
			  xlabel,ylabel,numBins);

  TGraph *nieves = new TGraph(nievesData);
  
  nieves->GetHistogram()->SetMinimum(0.0);
  nieves->GetHistogram()->SetMaximum(ymax);
  nieves->Draw("C");

  return hst_Y;
}

/*
Plot the differential cross section dSigma/dY where Y is some variable stored in the root file
-file is the filename of the first root file to be graphed
-type is one of N_RFG, N_LFG, LS_RFG, LS_LFG for file1
-spline should be a text file with the energy and xsec extracted from a .xml spline file
 (can also be a root file created with gspl2root)
-target is of the form C12, Pb208, etc for file1
-Enu is the neutrino beam energy in GeV for file1
-cmin and cmax are bounds on cos(theta), with theta the outoging lepton angle
-color is the color of the line
*/
TH1D* xsec_vs_Y(double ymax,TString file,TString splineName,TString target,TString Enu,
	      TString Y,double Ymin,double Ymax,
	      TString cmin = -1.0,TString cmax = 1.0,
	      bool same = false,int color = kRed,
	      TString xlabel = "",TString ylabel = "",
	      int numBins = 100)
{
  // Get the total cross section at the current energy
  Double_t totalXSec;
  if(splineName.EndsWith(".root")){
    cout << "Spline is a root file" << endl;
    TFile * spline = new TFile(splineName); //Get root spline to get total xsec
    TDirectory * dir1 = (TDirectory*) spline->Get("nu_mu_" + target);
    TGraph * graph1 = (TGraph*) dir1->Get("qel_cc_n");
    totalXSec = graph1->Eval(Enu.Atof());
  }else{
    cout << "Spline is a text file. Attempting to create spline from given points." << endl;
    // Assume data is formatted as required to make a TGraph
    TGraph * tempGraph = new TGraph(splineName);
    double factor = 3.90e10; // Splines from xml file are smaller than the root file by factor
    totalXSec = factor*tempGraph->Eval(Enu.Atof());
    cout << "Total XSec = " << totalXSec << endl;
  }

  // create histogram binned according to Y
  TChain * chain1 = new TChain("gst");
  chain1->Add(file);

  // Put togther title for curve
  TString title1 = "dSigma/d(" + Y + ") for " + target + " target, " + Enu + " GEV";

  // Parameters for TH1D are name, title, number bins, lower bound, upper bound
  TH1D* hst_Y = new TH1D("hst_Y",title1,numBins,Ymin,Ymax);

  if(same){
    chain1->Draw(Y+">>hst_Y","qel&&cthl>="+cmin+"&&cthl<="+cmax,"same");
  }else{
    chain1->Draw(Y+">>hst_Y","qel&&cthl>="+cmin+"&&cthl<="+cmax);
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

  return hst_Y;
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
