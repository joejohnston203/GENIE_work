#include<iostream>
#include<iomanip>
#include<stdio>
#include<fstream>
using namespace std;

void Drawplots_lfg(){
  //files for all CC processes, Eneu<10 GeV
  //string file1name = "June_14_12_numu_C_le_CC_v264.gst.root";
  string file1name = "Mar_31_16_numu_Pb_200_CCQE_v210_RFG_hA.gst.root";
  //  string file2name = "June_14_12_numu_C_me_CC_v264.gst.root";
  //  string file2name = "Jan_17_15_numu_C_2_v290_hAo.gst.root";
  string file2name = "Mar_31_16_numu_Pb_200_CCQE_v210_LFG_hA.gst.root";
  //string file2name = "Jun_11_12_numu_C_me_CC_v271.gst.root";
  
  //legend from minerva_1pi.C
  //  gROOT->SetStyle("T2k");
  //  gStyle->SetOptTitle(0);

  TLegend* leg1 = new TLegend(.7,.8, 1.0,1.0,"");
  TH1D * fake1 = new TH1D();
  TH1D * fake2 = new TH1D();
  TH1D * fake3 = new TH1D();
  //  fake1 -> SetMarkerStyle(20);
  //  fake1 -> SetMarkerSize(1);
  //  fake1 -> SetMarkerColor(kBlack);
  //  fake1 -> SetLineStyle(1);
  //  fake1 -> SetLineSize(2);
  //  fake1 -> SetLineColor(1);
  fake1 -> SetLineStyle(1);
  fake1 -> SetLineWidth(2);
  fake1 -> SetLineColor(2);
  fake2 -> SetLineStyle(1);
  fake2 -> SetLineWidth(2);
  fake2 -> SetLineColor(4);
  leg1->AddEntry(fake1,"GENIE 2.10+RFG","L");
  leg1->AddEntry(fake2,"GENIE 2.10+LFG","L");
  //  leg1->AddEntry(fake3,"GENIE 2.9 BS MB FF","L");

 TH1I * Q21qel = new TH1I("Q21","Q2 LE qel",100,0,.5);
 TH1I * Q22qel = new TH1I("Q22","Q2 ME qel",100,0,.5);
 TH1I * W1qel = new TH1I("W1","W qel",90,0,2);
 TH1I * W2qel = new TH1I("W2","W qel",90,0,2);
 TH1I * Enu1qel = new TH1I("Enu1","Enu LE qel",100,0,2);
 TH1I * Enu2qel = new TH1I("Enu2","Enu ME qel",100,0,2);
 TH1I * Q21res = new TH1I("Q21","Q2 LE res",100,0,10);
 TH1I * Q22res = new TH1I("Q22","Q2 ME res",100,0,10);
 TH1I * W1res = new TH1I("W1","W LE res",60,0,4);
 TH1I * W2res = new TH1I("W2","W ME res",60,0,4);
 TH1I * Enu1res = new TH1I("Enu1","Ev LE res",250,0,25);
 TH1I * Enu2res = new TH1I("Enu2","Ev ME res",250,0,25);
 TH1I * muonKEqel1 = new TH1I("muonKE1","Muon KE qel",200,0,.2);
 TH1I * muonKEqel2 = new TH1I("muonKE2","Muon KE qel",200,0,.2);
 TH1I * muonMomqel1 = new TH1I("muonMom1","Muon Mom qel",200,0,.2);
 TH1I * muonMomqel2 = new TH1I("muonMom2","Muon Mom qel",200,0,.2);
 TH1I * FermiMomqel1 = new TH1I("FermiMomqel1","Nucleon inital mom qel",200,0,.5);
 TH1I * FermiMomqel2 = new TH1I("FermiMomqel2","Nucleon inital mom qel",200,0,.5);
 TH1I * protonKEqel1 = new TH1I("protonKE1","Proton KE qel",100,0,.2);
 TH1I * protonKEqel2 = new TH1I("protonKE2","Proton KE qel",100,0,.2);
 TH1I * protonMomqel1 = new TH1I("protonMom1","Proton Mom qel",100,0,.6);
 TH1I * protonMomqel2 = new TH1I("protonMom2","Proton Mom qel",100,0,.6);
 TH1I * pionKEres1 = new TH1I("pionKE1","Pion KE res LE",35,0,3.5);
 TH1I * nfqel1 = new TH1I("nf1","nf qel LE",60,0,6);
 TH1I * nfres1 = new TH1I("nf1","nf res LE",100,0,13);
 TH1I * pionKEres2 = new TH1I("pionKE2","Pion KE res ME",100,0,10);
 TH1I * nfqel2 = new TH1I("nf2","nf qel ME",100,0,50);
 TH1I * nfres2 = new TH1I("nf2","nf res ME",100,0,50);
 /*
 TH1F * test = new TH1F("test","4-vector of (nu+p-mu-pi)^2",200,0,4);
 TH1F * Enudiff1 = new TH1F("Enudiff1","Enu err from total hadron energy",200,-2,2);
 TH2F * Etest1 = new TH2F("Etest1","Enu vs. Total hadronic energy",500,0,10,500,0,10);
 TH1F * Enudiff2 = new TH1F("Enudiff2","Enu err real - from muon",200,-2,2);
 TH2F * Etest2 = new TH2F("Etest2","Enu vs.Enu from muon ",500,0,10,500,0,10);
 TH2F * Enudiff_v_W = new TH2F("Enudiff_v_W","Enudiff vs.W ",500,0,5,500,-2,5);
 TH2F * Enudiff_v_KErecoil = new TH2F("Enudiff_v_KErecoil","Enudiff vs.Erecoil ",200,0,.1,500,-2,5);
 */

 
 TFile *file1 = new TFile(file1name.c_str(),"file1name");
 TTree *gst1 = (TTree*)file1->Get("gst");
 TFile *file2 = new TFile(file2name.c_str(),"file2name");
 TTree *gst2 = (TTree*)file2->Get("gst");
 const long n_entries1 = gst1->GetEntries();
 const long n_entries2 = gst2->GetEntries();
 const double Mp=0.9383, Mmu= 0.1057, MD=1.232, Mpi=0.1396;
 int nf1,nf2,nfpip1,nfpip2,nfp1,nfp2,nfn1,nfn2,nfpi01,nfpi02,nfpim1,nfpim2,pdgf1[100],pdgf2[100];
 double W1,W2,Ws1,Ws2,Q21,Q22,Enu1,Enu2,El1,El2,Emu1,Emu2,pxl1,pxl2,pyl1,pyl2,pzl1,pzl2,pl1,pl2,ptotl1,ptotl2,cthl1=0.,cthl2=0.,pxn1,pxn2,pyn1,pyn2,pzn1,pzn2;
 double Epi1,Epi2,ptotpi1,ptotpi2,Ef1[100],Ef2[100],pxf1[100],pxf2[100],pyf1[100],pyf2[100],pzf1[100],pzf2[100],pxn1,pxn2,pyn1,pyn2,pzn1,pzn2,pf1[100],pf2[100];
 double cthpi1=0.,cthpi2=0.,cthpimu1=0.,cthpimu2=0.,Mmupi21=0.,Mmupi22=0.,Mrecoil1=0.,Mrecoil2=0.,KErecoil1=0.,KErecoil2=0.,ptotf1=0.,ptotf2=0.;
 bool qel1,qel2,res1,res2,dis1,dis2,coh1,coh2;

 //gst->Draw("Ev","Ws<2&&!qel");

 //initialization of variables
   
 gst1->SetBranchAddress("qel",&qel1);
 gst1->SetBranchAddress("res",&res1);
 gst1->SetBranchAddress("dis",&dis1);
 gst1->SetBranchAddress("coh",&coh1);
 gst1->SetBranchAddress("nf",&nf1);
 gst1->SetBranchAddress("nfp",&nfp1);
 gst1->SetBranchAddress("nfn",&nfn1);
 gst1->SetBranchAddress("nfpip",&nfpip1);
 gst1->SetBranchAddress("nfpi0",&nfpi01);
 gst1->SetBranchAddress("nfpim",&nfpim1);
 gst1->SetBranchAddress("Ev",&Enu1);
 gst1->SetBranchAddress("El",&Emu1);
 gst1->SetBranchAddress("pl",&pl1);
 gst1->SetBranchAddress("pxl",&pxl1);
 gst1->SetBranchAddress("pyl",&pyl1);
 gst1->SetBranchAddress("pzl",&pzl1);
 gst1->SetBranchAddress("W",&W1);
 gst1->SetBranchAddress("Ws",&Ws1);
 gst1->SetBranchAddress("Q2",&Q21);
 gst1->SetBranchAddress("Ef",&Ef1[0]);
 gst1->SetBranchAddress("pf",&pf1[0]);
 gst1->SetBranchAddress("pxf",&pxf1[0]);
 gst1->SetBranchAddress("pyf",&pyf1[0]);
 gst1->SetBranchAddress("pzf",&pzf1[0]);
 gst1->SetBranchAddress("pxn",&pxn1);
 gst1->SetBranchAddress("pyn",&pyn1);
 gst1->SetBranchAddress("pzn",&pzn1);
 gst1->SetBranchAddress("pdgf",&pdgf1[0]);

 gst2->SetBranchAddress("qel",&qel2);
 gst2->SetBranchAddress("res",&res2);
 gst2->SetBranchAddress("dis",&dis2);
 gst2->SetBranchAddress("coh",&coh2);
 gst2->SetBranchAddress("nf",&nf2);
 gst2->SetBranchAddress("nfp",&nfp2);
 gst2->SetBranchAddress("nfn",&nfn2);
 gst2->SetBranchAddress("nfpip",&nfpip2);
 gst2->SetBranchAddress("nfpi0",&nfpi02);
 gst2->SetBranchAddress("nfpim",&nfpim2);
 gst2->SetBranchAddress("Ev",&Enu2);
 gst2->SetBranchAddress("El",&Emu2);
 gst2->SetBranchAddress("pl",&pl2);
 gst2->SetBranchAddress("pxl",&pxl2);
 gst2->SetBranchAddress("pyl",&pyl2);
 gst2->SetBranchAddress("pzl",&pzl2);
 gst2->SetBranchAddress("W",&W2); 
 gst2->SetBranchAddress("Ws",&Ws2);
 gst2->SetBranchAddress("Q2",&Q22);
 gst2->SetBranchAddress("Ef",&Ef2[0]);
 gst2->SetBranchAddress("pf",&pf2[0]);
 gst2->SetBranchAddress("pxf",&pxf2[0]);
 gst2->SetBranchAddress("pyf",&pyf2[0]);
 gst2->SetBranchAddress("pzf",&pzf2[0]);
 gst2->SetBranchAddress("pxn",&pxn2);
 gst2->SetBranchAddress("pyn",&pyn2);
 gst2->SetBranchAddress("pzn",&pzn2);
 gst2->SetBranchAddress("pdgf",&pdgf2[0]);


 for(int j=0;j<n_entries1;j++){
   gst1->GetEntry(j);
   //cout << j << " ";
   if (qel1==1){
     Q21qel->Fill(Q21);
     W1qel->Fill(W1);
     Enu1qel->Fill(Enu1);
     FermiMomqel1->Fill(TMath::Sqrt(pxn1*pxn1+pyn1*pyn1+pzn1*pzn1));
     muonKEqel1->Fill(Emu1-.106);
     muonMomqel1->Fill(pl1);
     nfqel1->Fill(nf1);
   }

   for(int jj=0;jj<nf1;jj++){
     if (pdgf1[jj]==2212)
       if (qel1==1){
	 protonKEqel1->Fill(Ef1[jj]-.9383);
	 protonMomqel1->Fill(pf1[jj]);
       }
   }
 }

 for(int k=0;k<n_entries2;k++){
   //cout << k << " ";
   gst2->GetEntry(k);
   if (qel2==1){
     Q22qel->Fill(Q22);
     W2qel->Fill(W2);
     Enu2qel->Fill(Enu2);
     FermiMomqel2->Fill(TMath::Sqrt(pxn2*pxn2+pyn2*pyn2+pzn2*pzn2));
     muonKEqel2->Fill(Emu2-.106);
     muonMomqel2->Fill(pl2);
     nfqel2->Fill(nf2);
   }
   for(int jj=0;jj<nf2;jj++){
     if (pdgf2[jj]==2212){
       if (qel2==1){
	 protonKEqel2->Fill(Ef2[jj]-.9383);
	 protonMomqel2->Fill(pf2[jj]);
       }
     }
   }

 }

 gStyle->SetPalette(1); 
 TCanvas * c1 = new TCanvas("c1","c1",800,600);
 Q21qel->SetStats(0);
 Q22qel->SetStats(0);
 Q21qel->SetTitle("Q2 qel");
 Q21qel->SetLineWidth(2);
 Q22qel->SetLineWidth(2);
 Q21qel->SetLineColor(2);
 Q22qel->SetLineColor(4);
 Q21qel->GetXaxis()->SetTitle("Q^2 (GeV^2)");
 Q21qel->GetYaxis()->SetTitle("Counts");
 Q21qel->Draw();
 Q22qel->Draw("same");
 leg1->Draw();
 c1->SaveAs("Q2qel.png");
 
 TCanvas * c2 = new TCanvas("c2","c2",800,600);
 W1qel->SetStats(0);
 W2qel->SetStats(0);
 W1qel->SetTitle("W qel");
 W1qel->SetLineWidth(2);
 W2qel->SetLineWidth(2);
 W1qel->SetLineColor(2);
 W2qel->SetLineColor(4);
 W1qel->GetXaxis()->SetTitle("W (GeV)");
 W1qel->GetYaxis()->SetTitle("Counts");
 W1qel->Draw();
 W2qel->Draw("same");
 leg1->Draw();
 c2->SaveAs("Wqel.png");

 TCanvas * c3 = new TCanvas("c3","c3",800,600);
 Enu1qel->SetStats(0);
 Enu2qel->SetStats(0);
 Enu1qel->SetTitle("Enu qel");
 Enu1qel->SetLineWidth(2);
 Enu2qel->SetLineWidth(2);
 Enu1qel->SetLineColor(2);
 Enu2qel->SetLineColor(4);
 Enu1qel->GetXaxis()->SetTitle("Enu (GeV)");
 Enu1qel->GetYaxis()->SetTitle("Counts");
 Enu1qel->Draw();
 Enu2qel->Draw("same");
 leg1->Draw();
 c3->SaveAs("Enuqel.png");

 TCanvas * c4 = new TCanvas("c4","c4",800,600);
 muonKEqel1->SetStats(0);
 muonKEqel2->SetStats(0);
 muonKEqel1->SetTitle("muon KE qel");
 muonKEqel1->SetLineWidth(2);
 muonKEqel2->SetLineWidth(2);
 muonKEqel1->SetLineColor(2);
 muonKEqel2->SetLineColor(4);
 muonKEqel1->GetXaxis()->SetTitle("E (GeV)");
 muonKEqel1->GetYaxis()->SetTitle("Counts");
 muonKEqel1->Draw();
 muonKEqel2->Draw("same");
 leg1->Draw();
 c4->SaveAs("muonKEqel.png");

 TCanvas * c5 = new TCanvas("c5","c5",800,600);
 protonKEqel1->SetStats(0);
 protonKEqel2->SetStats(0);
 protonKEqel1->SetTitle("Proton KE qel");
 protonKEqel1->SetLineWidth(2);
 protonKEqel2->SetLineWidth(2);
 protonKEqel1->SetLineColor(2);
 protonKEqel2->SetLineColor(4);
 protonKEqel1->GetXaxis()->SetTitle("E (GeV)");
 protonKEqel1->GetYaxis()->SetTitle("Counts");
 protonKEqel1->Draw();
 protonKEqel2->Draw("same");
 leg1->Draw();
 c5->SaveAs("protonKEqel.png");


 TCanvas * c6 = new TCanvas("c6","c6",800,600);
 muonMomqel1->SetStats(0);
 muonMomqel2->SetStats(0);
 muonMomqel1->SetTitle("muon Momentum qel");
 muonMomqel1->SetLineWidth(2);
 muonMomqel2->SetLineWidth(2);
 muonMomqel1->SetLineColor(2);
 muonMomqel2->SetLineColor(4);
 muonMomqel1->GetXaxis()->SetTitle("E (GeV)");
 muonMomqel1->GetYaxis()->SetTitle("Counts");
 muonMomqel1->Draw();
 muonMomqel2->Draw("same");
 leg1->Draw();
 c6->SaveAs("muonMomqel.png");

 TCanvas * c7 = new TCanvas("c7","c7",800,600);
 protonMomqel1->SetStats(0);
 protonMomqel2->SetStats(0);
 protonMomqel1->SetTitle("Proton Momentum qel");
 protonMomqel1->SetLineWidth(2);
 protonMomqel2->SetLineWidth(2);
 protonMomqel1->SetLineColor(2);
 protonMomqel2->SetLineColor(4);
 protonMomqel1->GetXaxis()->SetTitle("p (GeV)");
 protonMomqel1->GetYaxis()->SetTitle("Counts");
 c7->SetLogy();
 protonMomqel1->Draw();
 protonMomqel2->Draw("same");
 leg1->Draw();
 c7->SaveAs("protonMomqel.png");

 TCanvas * c8 = new TCanvas("c8","c8",800,600);
 nfqel1->SetStats(0);
 nfqel2->SetStats(0);
 nfqel1->SetTitle("nf qel");
 nfqel1->SetLineWidth(2);
 nfqel2->SetLineWidth(2);
 nfqel1->SetLineColor(2);
 nfqel2->SetLineColor(4);
 nfqel1->GetXaxis()->SetTitle("Number of particles");
 nfqel1->GetYaxis()->SetTitle("Counts");
 nfqel1->Draw();
 nfqel2->Draw("same");
 leg1->Draw();
 c8->SaveAs("nfqel.png");

 TCanvas * c9 = new TCanvas("c9","c9",800,600);
 FermiMomqel1->SetStats(0);
 FermiMomqel2->SetStats(0);
 FermiMomqel1->SetTitle("Fermi Mom qel");
 FermiMomqel1->SetLineWidth(2);
 FermiMomqel2->SetLineWidth(2);
 FermiMomqel1->SetLineColor(2);
 FermiMomqel2->SetLineColor(4);
 FermiMomqel1->GetXaxis()->SetTitle("E (GeV)");
 FermiMomqel1->GetYaxis()->SetTitle("Counts");
 FermiMomqel1->Draw();
 FermiMomqel2->Draw("same");
 leg1->Draw();
 c9->SaveAs("FermiMomqel.png");

 delete file1;
 delete file2;
}

int main()
{
  Drawplots_lfg();
  return 0;
}
