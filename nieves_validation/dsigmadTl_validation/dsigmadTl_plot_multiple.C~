#include <TStyle.h>
#include <TLatex.h>
#include "dsigmadTl_make_plots.C"

void dsigmadTl_plot_multiple() {
  // LS, LFG, 1 GeV, C12
  // do not inclue .png extension in output_file_name
  TString output_file_name = "LS_LFG_1GeV_C12";
  // location of .gst.root file from gevgen
  TString rootFile = "~/GENIE_work/runs/2016_03_17/2016_03_17_C12_LS_LFG_1GeV/numu_C12_CCQE_run1.gst.root";
  // 2 column text file (as above) or .root spline file
  TString spline = "~/GENIE_work/runs/2016_03_17/2016_03_17_C12_LS_LFG_1GeV/data_2016_03_17_C12_LS_LFG.txt";
  TString target = "C12"; // C12 or Pb208
  TString Enu = "1.0"; // in GeV- 0.2, 1.0, or 5.0 
  TString RPA = "0";   // 0 or 1
  TString Coul = "0";  // 0 or 1
  double ymax1 = 25.0; // y axis max for plot where all solid angles are integrated
  double ymax2 = 12.0; // y axis max for plot with specific ctl values
  // "" gives a default title
  TString plotTitle = "C12, Enu = 1.0 GeV, LwlynSmith, LFG"; 
  make_plots(output_file_name, rootFile, spline, target,
	     Enu, RPA, Coul, ymax1, ymax2,plotTitle);

  // LS, LFG, 200 MeV, C12
  // do not inclue .png extension in output_file_name
  TString output_file_name = "LS_LFG_200MeV_C12";
  // location of .gst.root file from gevgen
  TString rootFile = "~/GENIE_work/runs/2016_03_17/2016_03_21_C12_LS_LFG_200MeV/numu_C12_CCQE_run0.gst.root";
  // 2 column text file (as above) or .root spline file
  TString spline = "~/GENIE_work/runs/2016_03_17/2016_03_21_C12_LS_LFG_200MeV/data_2016_03_21_C12_LS_LFG_200MeV.txt";
  TString target = "C12"; // C12 or Pb208
  TString Enu = "0.2"; // in GeV- 0.2, 1.0, or 5.0 
  TString RPA = "0";   // 0 or 1
  TString Coul = "0";  // 0 or 1
  double ymax1 = 30.0; // y axis max for plot where all solid angles are integrated
  double ymax2 = 2.0; // y axis max for plot with specific ctl values
  // "" gives a default title
  TString plotTitle = "C12, Enu = 200 MeV, LwlynSmith, LFG"; 
  make_plots(output_file_name, rootFile, spline, target,
	     Enu, RPA, Coul, ymax1, ymax2,plotTitle);
  
  exit();
}
