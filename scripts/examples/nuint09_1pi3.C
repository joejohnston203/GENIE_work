//
// NuINT09 Conference, Benchmark Calculations (GENIE contribution)
//
// 1PI.3:
// d^Sigma / dOmega_pi+ dKE_pi+ at E_nu= 1.0 and 1.5 GeV
//
// Costas Andreopoulos, STFC / Rutherford Appleton Laboratory
//
#include <iomanip>

//
// consts
//
const int kNSamples       = 3;
const int kNWCur          = 1;
const int kNEnergies      = 2;
const int kNRunsPerCase   = 5;
const int kNEvtPerRun     = 100000;

const char * kLabel[kNSamples] = 
{  
 // available samples
 /* 0 */ "nu_mu_C12",
 /* 1 */ "nu_mu_O16",
 /* 2 */ "nu_mu_Fe56"
};
const int kRunNu1PI3[kNSamples][kNWCur][kNEnergies][kNRunsPerCase] =
{
 /* indices : sample ; cc/nc ; energy */
 {
  /* 0,0,0 (nu_mu C12,  CC, 1.0 GeV) */ { {900200, 900201, 900202, 900203, 900204},
  /* 0,0,1 (nu_mu C12,  CC, 1.5 GeV) */   {900300, 900301, 900302, 900303, 900304} }
 },
 {
  /* 1,0,0 (nu_mu O16,  CC, 1.0 GeV) */ { {910200, 910201, 910202, 910203, 910204},
  /* 1,0,1 (nu_mu O16,  CC, 1.5 GeV) */   {910300, 910301, 910302, 910303, 910304} }
 },
 {
  /* 2,0,0 (nu_mu Fe56, CC, 1.0 GeV) */ { {920200, 920201, 920202, 920203, 920204},
  /* 2,0,1 (nu_mu Fe56, CC, 1.5 GeV) */   {920300, 920301, 920302, 920303, 920304} }
 }
};
int kA[kNSamples] = 
{  
 // A for nuclear target at each sample
 /* 0 */  12,
 /* 1 */  16,
 /* 2 */  56
};

void nuint09_1pi3(int isample)
{
  if(isample<0 || isample >= kNSamples) return;

  const char * label = kLabel[isample];

  int A = kA[isample];

  // get cross section graphs

  TFile fsig("./cc1pip_tmp.root","read"); // generated at [1PI.1]
  TGraph * sig_graph_cc1pip = (TGraph*) fsig.Get("CC1pip");

  // range & spacing

  const int    nKEpip   = 60;
  const double KEpipmin =  0.00;
  const double KEpipmax =  1.50;

  const int    ncosth   = 30;
  const double costhmin = -1;
  const double costhmax = +1;

  // create output stream

  ostringstream out_filename;
  out_filename << label << ".1pi_3.d2sig1pi_dKEpidOmega.data";
  ofstream out_stream(out_filename.str().c_str(), ios::out);

  // write out txt file

  out_stream << "# [" << label << "]" << endl;
  out_stream << "#  " << endl;
  out_stream << "# [1PI.3]:" << endl;
  out_stream << "#  d^Sigma / dOmega_pi+ dKE_pi+ at E_nu= 1.0 and 1.5 GeV" << endl;
  out_stream << "#  " << endl;
  out_stream << "#  Note:" << endl;
  out_stream << "#   - pi+ kinetic energy KE in GeV, linear spacing between KEmin = " << KEpipmin << " GeV, KEmax = " << KEpipmax << " GeV "  << endl;
  out_stream << "#   - cross sections in 1E-38 cm^2 " << endl;
  out_stream << "#   - quoted cross section is nuclear cross section divided with number of nucleons A" << endl;
  out_stream << "#  Columns:" << endl;
  out_stream << "#  |  KE(pi+)   |   cos(theta_pi+)   |  dsig(numu A -> mu- 1pi+ X; Enu = 1.0 GeV)   |  dsig(numu A -> mu- 1pi+ X; Enu = 1.5 GeV)  | "  << endl;

  out_stream << setiosflags(ios::fixed) << setprecision(6);

  //
  // load event data
  // 
    
  TChain * chain = new TChain("gst");
  
  // loop over CC/NC cases
  for(int iwkcur=0; iwkcur<kNWCur; iwkcur++) {
    // loop over energies
    for(int ie=0; ie<kNEnergies; ie++) {
       // loop over runs for current case
       for(int ir=0; ir<kNRunsPerCase; ir++) {
          // build filename
          ostringstream filename;
          int run_number = kRunNu1PI3[isample][iwkcur][ie][ir];
          filename << "../gst/gntp." << run_number << ".gst.root";
          // add to chain
          cout << "Adding " << filename.str() << " to event chain" << endl;
          chain->Add(filename.str().c_str());
       }
    }
  } 

  // 
  // get CC1pi+ cross sections at given energies for normalization purposes
  //
  double sig_cc1pip_1000MeV = 1. / A; //HERE
  double sig_cc1pip_1500MeV = 1. / A; //HERE
  
  //
  // book histograms
  //
  TH2D * hst_d2sig_dKEpipdOmg_1000MeV = new TH2D("hst_d2sig_dKEpipdOmg_1000MeV",
           "d2sig / dKEpi+ dOmega, numu A -> mu- 1pi+ X, Enu=1.0 GeV", nKEpip, KEpipmin, KEpipmax, ncosth, costhmin, costhmax);
  TH2D * hst_d2sig_dKEpipdOmg_1500MeV = new TH2D("hst_d2sig_dKEpipdOmg_1500MeV",
           "d2sig / dKEpi+ dOmega, numu A -> mu- 1pi+ X, Enu=1.5 GeV", nKEpip, KEpipmin, KEpipmax, ncosth, costhmin, costhmax);

  //
  // fill histograms
  //
  chain->Draw("pzi/sqrt(pzi*pzi+pyi*pyi+pxi*pxi):(Ei-0.139)>>hst_d2sig_dKEpipdOmg_1000MeV","cc&&Ev>0.99&&Ev<1.01&&pdgi==211&&nipip==1&&nipim==0&&nipi0==0","GOFF");
  chain->Draw("pzi/sqrt(pzi*pzi+pyi*pyi+pxi*pxi):(Ei-0.139)>>hst_d2sig_dKEpipdOmg_1500MeV","cc&&Ev>1.49&&Ev<1.51&&pdgi==211&&nipip==1&&nipim==0&&nipi0==0","GOFF");
                
  //
  // normalize
  //
  double norm_cc1pip_1000MeV = hst_d2sig_dKEpipdOmg_1000MeV -> Integral("width") / (2*TMath::Pi() * sig_cc1pip_1000MeV);
  double norm_cc1pip_1500MeV = hst_d2sig_dKEpipdOmg_1500MeV -> Integral("width") / (2*TMath::Pi() * sig_cc1pip_1500MeV);
  
  if (norm_cc1pip_1000MeV > 0) hst_d2sig_dKEpipdOmg_1000MeV -> Scale(1./norm_cc1pip_1000MeV);
  if (norm_cc1pip_1500MeV > 0) hst_d2sig_dKEpipdOmg_1500MeV -> Scale(1./norm_cc1pip_1500MeV);

  for(int i = 1; i <= hst_d2sig_dKEpipdOmg_1500MeV->GetNbinsX(); i++) {
     for(int j = 1; j <= hst_d2sig_dKEpipdOmg_1500MeV->GetNbinsY(); j++) {
          
        double KEpip    = hst_d2sig_dKEpipdOmg_1500MeV -> GetXaxis() -> GetBinCenter(i);
        double costhpip = hst_d2sig_dKEpipdOmg_1500MeV -> GetYaxis() -> GetBinCenter(j);

        double d2sig_dKEpipdOmg_1000MeV = hst_d2sig_dKEpipdOmg_1000MeV -> GetBinContent(i,j);
        double d2sig_dKEpipdOmg_1500MeV = hst_d2sig_dKEpipdOmg_1500MeV -> GetBinContent(i,j);

        d2sig_dKEpipdOmg_1000MeV = TMath::Max(0., d2sig_dKEpipdOmg_1000MeV);
        d2sig_dKEpipdOmg_1500MeV = TMath::Max(0., d2sig_dKEpipdOmg_1500MeV);

        out_stream << setw(15) << KEpip 
                   << setw(15) << costhpip
                   << setw(15) << d2sig_dKEpipdOmg_1000MeV
                   << setw(15) << d2sig_dKEpipdOmg_1500MeV
                   << endl;
     }
  }

  out_stream.close();

  // visual inspection
  //TCanvas * c1 = new TCanvas("c1","",20,20,500,500);
  //...
  //c1->Update();

  fsig.Close();
  delete chain;
}
