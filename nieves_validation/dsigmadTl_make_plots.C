#include <TStyle.h>
#include <TLatex.h>
#include "dsigmadTl_script.C"

void dsigmadTl_make_plots() {
  make_plots();
}

void make_plots() {
  // Set paramters for plotting here
  TString graph_title = "LS_LFG";
  double ymax = 25.0;
  TString rootFile = "LS_LFG_C12_E1GeV.gst.root";
  TString spline = "LS_LFG_C12_E1GeV.data.txt";
  TString target = "C12";
  double Enu = 1.0;
  TString RPA = "0";
  TString Coul = "0";
  
  double ctlminmax[]={-1.0,1.0};
  int ctlsize = 2;
  save_plot(graph_title + "_all_anlges.png",ymax,rootFile,spline,target,Enu,RPA,Coul,
	    ctlminmax,ctlsize);

  double ctlminmax2[]={-1.0,-0.9,-0.05,0.05,0.4,0.5,0.7,0.8,0.9,0.95};
  int ctlsize2 = 10;
  save_plot(graph_title + ".png",ymax,rootFile,spline,target,Enu,RPA,Coul,
	    ctlminmax2,ctlsize2);
}

