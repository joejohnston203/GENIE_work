//
// NuINT11 Conference, GENIE-MiniBooNE comparisons
//
// CCQE
// d2Sigma/dTmudCosmu
// dSigma/dQ2
//
// Brandon Eberly, University of Pittsburgh
// Costas Andreopoulos, STFC / Rutherford Appleton Laboratory
//
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

//
// consts
//
const int kNSamples       = 1;
const int kNRunsPerCase   = 5;
const int kNEvtPerRun     = 100000;

const char * kLabel[kNSamples] = 
{  
 // available samples
 /* 0 */ "Tables/nu_mu_CH2"
};
const int kRunNuCCQE[kNSamples][kNRunsPerCase] = {1,2,3,4,5};

int kNeutrons[kNSamples] = 
{  
 // Number of neutrons in molecule for each sample
 /* 0, CH2 */  6
};

void nuint11_CCQElike_2D(int isample)
{
  if(isample<0 || isample >= kNSamples) return;

  const char * label = kLabel[isample];

  int nNeutrons = kNeutrons[isample];

  // get cross section graphs
  TFile fsig("./Splines/xsec_C12_H1.root","read");
  TGraph * sig_graph_ccC12 = (TGraph*) fsig.Get("nu_mu_C12/tot_cc");
  TGraph * sig_graph_ccH1 = (TGraph*) fsig.Get("nu_mu_H1/tot_cc");
  
  //make 3rd order splines from cross section graphs
  TSpline3* spl_ccC12 = new TSpline3("spl_ccC12",sig_graph_ccC12);
  TSpline3* spl_ccH1  = new TSpline3("spl_ccH1",sig_graph_ccH1);

  // range & spacing
  const int    nTMu   = 18;
  const double TMuMin =  0.2;
  const double TMuMax =  2.0;
  
  const int nCosMu = 20;
  const double CosMuMin = -1.0;
  const double CosMuMax = 1.0;
  
  const int nQ2 = 40;
  const double Q2Min = 0.0;
  const double Q2Max = 2.0;
  
  const int nEv = 100;
  const double  EvMin = 0.0;
  const double EvMax = 10.0;

  // load event data
  TChain * chain = new TChain("gst");
  for(int ir=0; ir<kNRunsPerCase; ir++) {
     // build filename
     ostringstream filename;
     int run_number = kRunNuCCQE[isample][ir];
     filename << "Sim/miniboone_gst" << run_number << ".root";
     // add to chain
     cout << "Adding " << filename.str() << " to event chain" << endl;
     chain->Add(filename.str().c_str());
  }
  
  // book histograms
  TH2D* hst_2D_ccqelike = new TH2D("hst_2D_ccqelike","Tmu v CosMu CCQE-like", nCosMu, CosMuMin, CosMuMax, nTMu, TMuMin, TMuMax);
  TH1D* hst_Q2_ccqelike = new TH1D("hst_Q2_ccqelike","Q2 CCQE-like", nQ2, Q2Min, Q2Max);
  TH1D* hst_2D_cctot = new TH1D("hst_2D_cctot","Ev CCtotal", nEv, EvMin, EvMax);

  // fill histograms
  chain->Draw("El - 0.1057:pzl/sqrt(El*El - 0.1057*0.1057)>>hst_2D_ccqelike","cc&&neu==14&&nfpip==0&&nfpim==0&&nfpi0==0&&nfem==0","GOFF");
  chain->Draw("Q2>>hst_Q2_ccqelike","cc&&neu==14&&nfpip==0&&nfpim==0&&nfpi0==0&&nfem==0","GOFF");
  long int nCCTot = chain->Draw("Ev>>hst_2D_cctot","cc&&neu==14","GOFF");
  
  //integrate flux for normalization
  ifstream fluxfile("MBData/minibooneFluxNuMu.txt");
  if (!fluxfile) {
    cout <<"Could not open miniboone flux file!"<<endl;
  }
  float fluxE, fluxV;
  float tot_flux = 0.0;
  float tot_fluxSig = 0.0;
  while (fluxfile) {
    fluxfile >> fluxE;
    fluxfile >> fluxV;
    tot_flux += fluxV;
    tot_fluxSig += fluxV*(spl_ccC12->Eval(fluxE+0.025) + 2.0*(spl_ccH1->Eval(fluxE+0.025)) );
  }
  tot_flux *= 0.05;  //50 MeV
  tot_fluxSig *= 0.05;
  cout <<"Total miniboone flux is "<<tot_flux<<endl;
  fluxfile.close();
  
  //bin widths
  const float binWTMu = (TMuMax - TMuMin)/nTMu;
  const float binWCosMu = (CosMuMax - CosMuMin)/nCosMu;
  const float binWQ2 = (Q2Max - Q2Min)/nQ2;
    
  // Normalization factor to number of nucleons in molecule
  double normalize = tot_fluxSig/(nNeutrons*binWTMu*binWCosMu*tot_flux*nCCTot); //units of 10^-38 cm^2/GeV
  double normalizeQ2 = tot_fluxSig/(nNeutrons*binWQ2*tot_flux*nCCTot);
  
  //To get CCQE-like cross section, loop over bins in hst_2D_ccqelike, select desired angle, multiply contents by norm factors
  double mu_T_0deg[nTMu];
  double xsec_0deg[nTMu]; 
  double mu_T_e_0deg[nTMu];
  double xsec_e_0deg[nTMu];
  
  double mu_T_30deg[nTMu];
  double xsec_30deg[nTMu]; 
  double mu_T_e_30deg[nTMu];
  double xsec_e_30deg[nTMu];
  
  //histogram is y vx x = TMu vs CosMu
  for(int i = 1; i <= hst_2D_ccqelike->GetNbinsX(); i++) {
    double currentCos  = hst_2D_ccqelike->GetXaxis()->GetBinCenter(i); 
    if (currentCos > 0.9 && currentCos <= 1.0) {
      for(int j = 1; j <= hst_2D_ccqelike->GetNbinsY(); j++) {

	double currentT  = hst_2D_ccqelike->GetYaxis()->GetBinCenter(j);
	double currentN = hst_2D_ccqelike->GetBinContent(i,j);

        mu_T_0deg[j-1] = currentT;
	xsec_0deg[j-1] = currentN*normalize;
	
	xsec_e_0deg[j-1] = normalize*(hst_2D_ccqelike->GetBinError(i,j));
	mu_T_e_0deg[j-1] = 0.5*binWTMu;   
      }
    }
    if (currentCos > 0.8 && currentCos <= 0.9) {
      for(int j = 1; j <= hst_2D_ccqelike->GetNbinsY(); j++) {

	double currentT  = hst_2D_ccqelike->GetYaxis()->GetBinCenter(j);
	double currentN = hst_2D_ccqelike->GetBinContent(i,j);

        mu_T_30deg[j-1] = currentT;
	xsec_30deg[j-1] = currentN*normalize;
	
	xsec_e_30deg[j-1] = normalize*(hst_2D_ccqelike->GetBinError(i,j));
	mu_T_e_30deg[j-1] = 0.5*binWTMu;   
      }
    }
  }
  
  TH2F* hst_2D_xsec = new TH2F("hst_2D_xsec","Tmu v CosMu CCQE-like", nCosMu, CosMuMin, CosMuMax, nTMu, TMuMin, TMuMax);
  for (int i=1; i<= hst_2D_ccqelike->GetNbinsX(); i++) {
    double currentCos = hst_2D_ccqelike->GetXaxis()->GetBinCenter(i);
    for (int j=1; j <= hst_2D_ccqelike->GetNbinsY(); j++) {
      double currentT = hst_2D_ccqelike->GetYaxis()->GetBinCenter(j);
      double currentN = hst_2D_ccqelike->GetBinContent(i,j);
      
      currentN *= normalize;
      hst_2D_xsec->Fill(currentCos, currentT, currentN);
    }
  }

  
  //Now for Q2, loop over bins in hst_Q2_ccqelike, multiply contents by norm factors
  double Q2[nQ2];
  double xsec_Q2[nQ2];
  double Q2_e[nQ2];
  double xsec_Q2_e[nQ2];
  
  for (int i = 1; i<= hst_Q2_ccqelike->GetNbinsX(); i++) {
    Q2[i-1] = hst_Q2_ccqelike->GetBinCenter(i);
    xsec_Q2[i-1] = hst_Q2_ccqelike->GetBinContent(i);
    
    xsec_Q2[i-1] *= normalizeQ2;
    
    Q2_e[i-1] = 0.5*binWQ2;
    xsec_Q2_e[i-1] = normalizeQ2*(hst_Q2_ccqelike->GetBinError(i));
  }
  
  //build genie graphs
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  TGraphErrors *gr_0deg = new TGraphErrors(nTMu, mu_T_0deg, xsec_0deg, mu_T_e_0deg, xsec_e_0deg);
  
  gr_0deg->SetTitle("CCQE-like d2(sig)/d(Tmu)d(cosThetaMu) Flux-Integrated Cross Section per Neutron at cosThetaMu > 0.9");
  gr_0deg->SetMarkerStyle(21);
  gr_0deg->SetMarkerColor(kRed);
  gr_0deg->SetMarkerSize(1.5);
  gr_0deg->GetXaxis()->SetTitle("Muon Kinetic Energy (GeV)");
  gr_0deg->GetXaxis()->CenterTitle();
  gr_0deg->GetYaxis()->SetTitle("d2(sig)/d(Tmu)d(cosThetaMu) (10^-38 cm^2/GeV)");
  gr_0deg->GetYaxis()->CenterTitle();  

  gr_0deg->Draw("AP");
  
  //build graph of miniboone data
  //data file has columns Tmu center, Tmu width, xsec, xsec_err, xsec_modeldep
  //xsec in units of 10-41 cm^2/GeV, Tmu in units of GeV
  double mb_Tmu_0deg[18], mb_Tmu_err_0deg[18], mb_xsec_0deg[18], mb_xsec_err_0deg[18];
  double mb_modeldep = 0.0;
  ifstream mb_0deg("MBData/minibooneCCQElike_2D_cosTheta0.9-1.0_data.txt");
  if (!mb_0deg) {
    cout <<"Could not open miniboone data file"<<endl;
  }
  int imb=0;
  while (mb_0deg && imb < 18) {
    mb_0deg >> mb_Tmu_0deg[imb];
    mb_0deg >> mb_Tmu_err_0deg[imb];
    mb_0deg >> mb_xsec_0deg[imb];
    mb_0deg >> mb_xsec_err_0deg[imb];
    mb_0deg >> mb_modeldep;
    
    //add 10.7% normalization error to xsec error
    mb_xsec_err_0deg[imb] = sqrt(mb_xsec_err_0deg[imb]*mb_xsec_err_0deg[imb] + 0.107*0.107*mb_xsec_0deg[imb]*mb_xsec_0deg[imb]);
    
    //add model dependence, convert to units of 10-38 cm^2/GeV
    mb_xsec_0deg[imb] += mb_modeldep;
    mb_xsec_0deg[imb] *= 0.001;
    mb_xsec_err_0deg[imb] *= 0.001;
    
    imb++;
  }
  
  TGraphErrors *mb_gr_0deg = new TGraphErrors(18, mb_Tmu_0deg, mb_xsec_0deg, mb_Tmu_err_0deg, mb_xsec_err_0deg);  
  mb_gr_0deg->SetTitle("CCQE-like d2(sig)/d(Tmu)d(cosThetaMu) Flux-Integrated Cross Section per Neutron at cosThetaMu > 0.9");
  mb_gr_0deg->SetMarkerStyle(22);
  mb_gr_0deg->SetMarkerColor(kBlue);
  mb_gr_0deg->SetMarkerSize(1.5);
  mb_gr_0deg->GetXaxis()->SetTitle("Muon Kinetic Energy (GeV)");
  mb_gr_0deg->GetXaxis()->CenterTitle();
  mb_gr_0deg->GetYaxis()->SetTitle("d2(sig)/d(Tmu)d(cosThetaMu) (10^-38 cm^2/GeV)");
  mb_gr_0deg->GetYaxis()->CenterTitle(); 
  
  mb_gr_0deg->Draw("P"); 
  
  //==================================================
  // 0.8 < cosTheta < 0.9
  //==================================================

  //build genie graphs
  TCanvas *c2 = new TCanvas("c2","c2",800,600);
  TGraphErrors *gr_30deg = new TGraphErrors(nTMu, mu_T_30deg, xsec_30deg, mu_T_e_30deg, xsec_e_30deg);
  
  gr_30deg->SetTitle("CCQE-like d2(sig)/d(Tmu)d(cosThetaMu) Flux-Integrated Cross Section per Neutron at 0.8 < cosThetaMu < 0.9");
  gr_30deg->SetMarkerStyle(21);
  gr_30deg->SetMarkerColor(kRed);
  gr_30deg->SetMarkerSize(1.5);
  gr_30deg->GetXaxis()->SetTitle("Muon Kinetic Energy (GeV)");
  gr_30deg->GetXaxis()->CenterTitle();
  gr_30deg->GetYaxis()->SetTitle("d2(sig)/d(Tmu)d(cosThetaMu) (10^-38 cm^2/GeV)");
  gr_30deg->GetYaxis()->CenterTitle();  

  gr_30deg->Draw("AP");
  
  //build graph of miniboone data
  //data file has columns Tmu center, Tmu width, xsec, xsec_err, xsec_modeldep
  //xsec in units of 10-41 cm^2, Tmu in units of GeV
  double mb_Tmu_30deg[16], mb_Tmu_err_30deg[16], mb_xsec_30deg[16], mb_xsec_err_30deg[16];
  mb_modeldep = 0.0;
  ifstream mb_30deg("MBData/minibooneCCQElike_2DcosTheta0.8-0.9_data.txt");
  if (!mb_30deg) {
    cout <<"Could not open miniboone data file"<<endl;
  }
  imb=0;
  while (mb_30deg && imb < 16) {
    mb_30deg >> mb_Tmu_30deg[imb];
    mb_30deg >> mb_Tmu_err_30deg[imb];
    mb_30deg >> mb_xsec_30deg[imb];
    mb_30deg >> mb_xsec_err_30deg[imb];
    mb_30deg >> mb_modeldep;
    
    //add 10.7% normalization error to xsec error
    mb_xsec_err_30deg[imb] = sqrt(mb_xsec_err_30deg[imb]*mb_xsec_err_30deg[imb] + 0.107*0.107*mb_xsec_30deg[imb]*mb_xsec_30deg[imb]);
    
    //add model dependence, convert to units of 10-38 cm^2
    mb_xsec_30deg[imb] += mb_modeldep;
    mb_xsec_30deg[imb] *= 0.001;
    mb_xsec_err_30deg[imb] *= 0.001;
    
    imb++;
  }
  
  TGraphErrors *mb_gr_30deg = new TGraphErrors(16, mb_Tmu_30deg, mb_xsec_30deg, mb_Tmu_err_30deg, mb_xsec_err_30deg);  
  mb_gr_30deg->SetTitle("CCQE-like d2(sig)/d(Tmu)d(cosThetaMu) Flux-Integrated Cross Section per Neutron at 0.8 < cosThetaMu < 0.9");
  mb_gr_30deg->SetMarkerStyle(22);
  mb_gr_30deg->SetMarkerColor(kBlue);
  mb_gr_30deg->SetMarkerSize(1.5);
  mb_gr_30deg->GetXaxis()->SetTitle("Muon Kinetic Energy (GeV)");
  mb_gr_30deg->GetXaxis()->CenterTitle();
  mb_gr_30deg->GetYaxis()->SetTitle("d2(sig)/d(Tmu)d(cosThetaMu) (10^-38 cm^2/GeV)");
  mb_gr_30deg->GetYaxis()->CenterTitle(); 
  
  mb_gr_30deg->Draw("P"); 

  //============================================================================
  // Q2
  //============================================================================
  TCanvas *c3 = new TCanvas("c3","c3",800,600);
  TGraphErrors *gr_Q2 = new TGraphErrors(nQ2, Q2, xsec_Q2, Q2_e, xsec_Q2_e);
  
  gr_Q2->SetTitle("CCQE-like d(sig)/d(Q^2) Flux-Integrated Cross Section per Neutron");
  gr_Q2->SetMarkerStyle(21);
  gr_Q2->SetMarkerColor(kRed);
  gr_Q2->SetMarkerSize(1.5);
  gr_Q2->GetXaxis()->SetTitle("Q^2 (GeV^2)");
  gr_Q2->GetXaxis()->CenterTitle();
  gr_Q2->GetYaxis()->SetTitle("d(sig)/d(Q^2) (10^-38 cm^2/GeV^2)");
  gr_Q2->GetYaxis()->CenterTitle();  

  gr_Q2->Draw("AP");
  
  //build graph of miniboone data
  //data file has columns Q^2 center, Q^2 width, xsec, xsec_err, xsec_modeldep
  //xsec in units of  cm^2, Q^2 in units of GeV
  double mb_Q2[17], mb_Q2_err[17], mb_xsec_Q2[17], mb_xsec_Q2_err[17];
  mb_modeldep = 0.0;
  ifstream mb_Q2_file("MBData/minibooneCCQElike_Q2_data.txt");
  if (!mb_Q2_file) {
    cout <<"Could not open miniboone data file"<<endl;
  }
  imb=0;
  while (mb_Q2_file && imb < 17) {
    mb_Q2_file >> mb_Q2[imb];
    mb_Q2_file >> mb_Q2_err[imb];
    mb_Q2_file >> mb_xsec_Q2[imb];
    mb_Q2_file >> mb_xsec_Q2_err[imb];
    mb_Q2_file >> mb_modeldep;
    
    //add 10.7% normalization error to xsec error
    mb_xsec_Q2_err[imb] = sqrt(mb_xsec_Q2_err[imb]*mb_xsec_Q2_err[imb] + 0.107*0.107*mb_xsec_Q2[imb]*mb_xsec_Q2[imb]);
    
    //add model dependence, convert to units of 10-38 cm^2
    mb_xsec_Q2[imb] += mb_modeldep;
    mb_xsec_Q2[imb] /= 1.0e-38;
    mb_xsec_Q2_err[imb] /= 1.0e-38;
    
    imb++;
  }
  
  TGraphErrors *mb_gr_Q2 = new TGraphErrors(17, mb_Q2, mb_xsec_Q2, mb_Q2_err, mb_xsec_Q2_err);  
  mb_gr_Q2->SetTitle("CCQE-like d(sig)/d(Q^2) Flux-Integrated Cross Section per Neutron");
  mb_gr_Q2->SetMarkerStyle(22);
  mb_gr_Q2->SetMarkerColor(kBlue);
  mb_gr_Q2->SetMarkerSize(1.5);
  mb_gr_Q2->GetXaxis()->SetTitle("Q^2 (GeV^2)");
  mb_gr_Q2->GetXaxis()->CenterTitle();
  mb_gr_Q2->GetYaxis()->SetTitle("d(sig)/d(Q^2) (10^-38 cm^2/GeV^2)");
  mb_gr_Q2->GetYaxis()->CenterTitle(); 
  
  mb_gr_Q2->Draw("P");
  
  //make 2d ratio plot
  TCanvas *c4 = new TCanvas("c4","c4",800,600);
  TH2F* hst_2D_MBxsec = new TH2F("hst_2D_MBxsec","Tmu v CosMu CCQE-like", nCosMu, CosMuMin, CosMuMax, nTMu, TMuMin, TMuMax);
  double mb_Tmu_bin[nTMu]={0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25,1.35,1.45,1.55,1.65,1.75,1.85,1.95};
  double mb_Cmu_bin[nCosMu]={0.95,0.85,0.75,0.65,0.55,0.45,0.35,0.25,0.15,0.05,-0.05,-0.15,-0.25,-0.35,-0.45,-0.55,-0.65,-0.75,-0.85,-0.95};
  int nline = 0;
  int ncol = 0;
  string line1, line2, line3;
  ifstream mb_all("MBData/minibooneCCQElike_2D_allxsec.txt");
  ifstream mb_bkg("MBData/minibooneCCQElike_2D_allbg.txt");
  ifstream mb_err("MBData/minibooneCCQElike_2D_allerr.txt");
  
  double mb_integral = 0.0;
  while (mb_all && mb_bkg && mb_err && nline < nCosMu) {
    getline(mb_all, line1);
    getline(mb_bkg, line2);
    getline(mb_err, line3);
    
    double xsec, bkg, err;
    ncol = 0;
    istringstream ostr1(line1);
    istringstream ostr2(line2);
    istringstream ostr3(line3);
    while (ncol < nTMu) {
      ostr1 >> xsec;
      ostr2 >> bkg;
      ostr3 >> err;

      xsec += bkg;
      xsec *= 0.001;
      err  *= 0.0001;
      
      mb_integral += xsec*binWTMu*binWCosMu;
      
      if (xsec!= 0.0 && err/xsec >= 0.224) {
        xsec = 0.0;
      }
      
      hst_2D_MBxsec->Fill(mb_Cmu_bin[nline], mb_Tmu_bin[ncol], xsec);
      ncol++;
    }
    nline++;
  }
  
  //get normalization of GENIE to MB
  double area_normfactor = hst_2D_xsec->Integral("width")/mb_integral;
  
  //now divide (why doesn't root have this in the TH2 class?)
  TH2F* ratio_hist = new TH2F("ratio_hist","CCQELike d(sig)/d(Tmu)d(cosThetaMu) Flux-Integrated Cross Section Ratio:  GENIE/MiniBooNE, Absolute Normalization",nTMu, TMuMin, TMuMax, nCosMu, CosMuMin, CosMuMax);
  TH2F* area_ratio_hist = new TH2F("area_ratio_hist","CCQELike d(sig)/d(Tmu)d(cosThetaMu) Flux-Integrated Cross Section Ratio:  GENIE/MiniBooNE, Area Normalization",nTMu, TMuMin, TMuMax, nCosMu, CosMuMin, CosMuMax);
  for (int i=1; i<= hst_2D_xsec->GetNbinsX(); i++) {
    double currentCos = hst_2D_xsec->GetXaxis()->GetBinCenter(i);
    for (int j=1; j <= hst_2D_xsec->GetNbinsY(); j++) {
      double currentT = hst_2D_xsec->GetYaxis()->GetBinCenter(j);
      double genie_xsec = hst_2D_xsec->GetBinContent(i,j);
      double mb_xsec = hst_2D_MBxsec->GetBinContent(i,j);
      
      double ratio = mb_xsec !=0.0 ? genie_xsec/mb_xsec : 0.0;
      double a_ratio = mb_xsec !=0.0 ? genie_xsec/(mb_xsec*area_normfactor) : 0.0;
      //cout<<"Cos: "<<currentCos<<"T: "<<currentT<<" genie_xsec: "<<genie_xsec<<" mb_xsec: "<<mb_xsec<<" ratio "<<ratio<<endl;

      ratio_hist->Fill(currentT, currentCos, ratio);
      area_ratio_hist->Fill(currentT, currentCos, a_ratio);
    }
  }  
  
  gStyle->SetPalette(1);
  ratio_hist->GetZaxis()->SetRangeUser(0.0,1.6);
  ratio_hist->GetXaxis()->SetTitle("Muon Kinetic Energy (GeV)");
  ratio_hist->GetYaxis()->SetTitle("Cos(ThetaMu)");
  ratio_hist->DrawCopy("COLZ");  
  
  TCanvas *c5 = new TCanvas("c5","c5",800,600);
  area_ratio_hist->GetZaxis()->SetRangeUser(0.0,1.6);
  area_ratio_hist->GetXaxis()->SetTitle("Muon Kinetic Energy (GeV)");
  area_ratio_hist->GetYaxis()->SetTitle("Cos(ThetaMu)");
  area_ratio_hist->DrawCopy("COLZ");
  

  fsig.Close();
  delete chain;
}
