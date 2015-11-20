// target is the target symbol (eg C12 or Pb208) and Type is LS or N and LFG
// or GFG, so options are N_GFG, N_LFG, LS_GFG, LS_LFG.
// Assumed filenames are of the form type_numu_target.gst.root
void graphCurve(TString type, TString target, bool same,
	       TString quantity = "Q2", TString cut = "qel"){
  TFile *file = new TFile(type + "_numu_" + target + ".gst.root");
  if(same){
    gst->Draw(quantity,cut,"same");
  }else{
    gst->Draw(quantity,cut);
  }
}

void graphTargetE1GeV(TString target, TString quantity = "Q2", TString cut = "qel"){
  TString types[4] = {"LS_GFG","LS_LFG","N_GFG","N_LFG"};
  Color_t colors[4] = {kBlue,kCyan,kRed,kMagenta};

  TH1F* curves[4];

  gStyle->SetOptStat(0);

  int i = 0;
  graphCurve(types[i],target,false, quantity, cut);
  curves[i] = (TH1F*) htemp->Clone();
  curves[i]->SetLineColor(colors[i]);
  curves[i]->SetTitle(types[i] + ", " + target + " target");

  int i = 1;
  graphCurve(types[i],target,false, quantity, cut);
  curves[i] = (TH1F*) htemp->Clone();
  curves[i]->SetLineColor(colors[i]);
  curves[i]->SetTitle(types[i] + ", " + target + " target");

  curves[0]->Draw();
  curves[1]->Draw("same");
  c1->SetTitle(target + " target, Eneutrino = 1 GeV");
  c1->BuildLegend();
  c1->SaveAs(target + "_E1GeV_N_" +quantity +"_" + cut + ".png");
  c1->Close();

  int i = 2;
  graphCurve(types[i],target,false, quantity, cut);
  curves[i] = (TH1F*) htemp->Clone();
  curves[i]->SetLineColor(colors[i]);
  curves[i]->SetTitle(types[i] + ", " + target + " target");

  int i = 3;
  graphCurve(types[i],target,false, quantity, cut);
  curves[i] = (TH1F*) htemp->Clone();
  curves[i]->SetLineColor(colors[i]);
  curves[i]->SetTitle(types[i] + ", " + target + " target");


  // for some reason the above code causes a segmentation violation the first
  // time it is run if it is put in a for loop
  /*for(int i = 0; i<4;i++){
    graphCurve(types[i],target,false);
    curves[i] = (TH1F*) htemp->Clone();
    curves[i]->SetLineColor(kRed);
    }*/

  curves[2]->Draw();
  curves[3]->Draw("same");

  c1->SetTitle(target + " target, Eneutrino = 1 GeV");
  c1->BuildLegend();
  c1->SaveAs(target + "_E1GeV_LS_" +quantity +"_" + cut + ".png");
  c1->Close();
}

// By default plot all four splines for n, Li6, C12, Cd112, and Pb208
void plotAllTargets(TString quantity = "Q2", TString cut = "qel"){
  TString targets[5] = {"neutron","Li6","C12","Cd112","Pb208"};
  for(int i=0;i<5;i++) graphTargetE1GeV(targets[i],quantity,cut);
}
/*
// Overloaded to work with up to 7 targets. Assumed there is a root file for LS_GFG,LS_LFG,N_GFG, and N_LFG
void plotAllTargets(TString tgt1,TString quantity = "Q2", TString cut = "qel"){
  TString targets[1] = {tgt1};
  graphTargetE1GeV(targets[i]);
}
void plotAllTargets(TString tgt1,TString tgt2,TString quantity = "Q2", TString cut = "qel"){
  TString targets[2] = {tgt1,tgt2};
  for(int i=0;i<2;i++) graphTargetE1GeV(targets[i]);
}
void plotAllTargets(TString tgt1,TString tgt2,TString tgt3,TString quantity = "Q2", TString cut = "qel"){
  TString targets[3] = {tgt1,tgt2,tgt3};
  for(int i=0;i<3;i++) graphTargetE1GeV(targets[i]);
}
void plotAllTarget(TString tgt1,TString tgt2,TString tgt3,TString tgt4,TString quantity = "Q2", 
		   TString cut = "qel"){
  TString targets[4] = {tgt1,tgt2,tgt3,tgt4};
  for(int i=0;i<4;i++) graphTargetE1GeV(targets[i]);
}
void plotAllTarget(TString tgt1,TString tgt2,TString tgt3,TString tgt4,TString tgt5,
		   TString quantity = "Q2", TString cut = "qel"){
  TString targets[5] = {tgt1,tgt2,tgt3,tgt4,tgt5};
  for(int i=0;i<5;i++) graphTargetE1GeV(targets[i]);
}
void plotAllTarget(TString tgt1,TString tgt2,TString tgt3,TString tgt4,TString tgt5,TString tgt6,
		   TString quantity = "Q2", TString cut = "qel"){
  TString targets[6] = {tgt1,tgt2,tgt3,tgt4,tgt5,tgt6};
  for(int i=0;i<6;i++) graphTargetE1GeV(targets[i]);
}
void plotAllTarget(TString tgt1,TString tgt2,TString tgt3,TString tgt4,TString tgt5,TString tgt6,
		   TString tgt7,TString quantity = "Q2", TString cut = "qel"){
  TString targets[7] = {tgt1,tgt2,tgt3,tgt4,tgt5,tgt6,tgt7};
  for(int i=0;i<7;i++) graphTargetE1GeV(targets[i]);
}
*/
