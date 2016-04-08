#include "Algorithm/AlgFactory.h"
#include "Algorithm/AlgConfigPool.h"
#include "Nuclear/NuclearModel.h"
#include "Nuclear/NuclearModelI.h"
#include "GHEP/GHepRecord.h"
#include "Interaction/Interaction.h"

// In whatever methods use the nuclear model, GenerateNucleon must be called
// with a radius if the model could be LFG.
void SomeMethod(GHepRecord * evrec){
  // For example, this code gets the position of the struck nucleon
  double radius = evrec->HitNucleon()->X4()->Vect().Mag();

  // Now generate the nucleon
  Target * tgt = evrec->Summary()->InitStatePtr()->TgtPtr();
  fNuclModel->GenerateNucleon(*tgt,radius);

  // Store the generated nucleon momentum and removal energy
  TVector3 p3 = fNuclModel->Momentum3();
  double w    = fNuclModel->RemovalEnergy();
}

// In LoadConfig, use the key "NuclearModel" to instantiate the default
// nuclear model
void LoadConfig(){
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();
  
  // Create a nuclear model object to check the model type
  RgKey nuclkey = "NuclearModel";
  RgAlg nuclalg = gc->GetAlg(nuclkey);
  AlgFactory * algf = AlgFactory::Instance();
  const NuclearModelI* fNuclModel = 
    dynamic_cast<const NuclearModelI*>(
			     algf->GetAlgorithm(nuclalg.name,nuclalg.config));

  // Stop here if the default nuclear model should be used. Continue to
  // check whether the nuclear model is LFG, and use LocalFGM if it is
  // and FGMBodekRitchie (Relativistic Fermi gas) otherwise

  // Check if the model is a local Fermi gas
  bool fLFG = (nuclModel && nuclModel->ModelType(Target()) == kNucmLocalFermiGas);
  if(fLFG){
    fNuclModel = dynamic_cast<const NuclearModelI *>
      (algf->GetAlgorithm("genie::FGMBodekRitchie","Default"));
  }else{
    fNuclModel = dynamic_cast<const NuclearModelI *>
      (algf->GetAlgorithm("genie::FGMBodekRitchie","Default"));
  }
}

// Code to calculate the fermi momentum using LFG or RFG depending on the
// chosen nuclear model in UserPhysicsOptions.xml
double KFermi(GHepRecord * evrec, double r){
  // Create a nuclear model object to check the model type
  RgKey nuclkey = "NuclearModel";
  RgAlg nuclalg = gc->GetAlg(nuclkey);
  AlgFactory * algf = AlgFactory::Instance();
  const NuclearModelI* nuclModel = 
    dynamic_cast<const NuclearModelI*>(
			     algf->GetAlgorithm(nuclalg.name,nuclalg.config));
  // Check if the model is a local Fermi gas
  bool fLFG = (nuclModel && nuclModel->ModelType(Target()) == kNucmLocalFermiGas);
  
  // Calculate the fermi momentum
  double kf;
  if(fLFG){
    Interaction * interaction = evrec->Summary();
    GHepParticle * hit = evrec->HitNucleon();
    int nucleon_pdgc = hit->Pdg();
    assert(pdg::IsProton(nucleon_pdgc) || pdg::IsNeutron(nucleon_pdgc));
    Target* tgt = interaction->InitStatePtr()->TgtPtr();
    int A = tgt->A();
    bool is_p = pdg::IsProton(nucleon_pdgc);
    double numNuc = (is_p) ? (double)tgt->Z():(double)tgt->N();
    double hbarc = kLightSpeed*kPlankConstant/fermi;
    kf= TMath::Power(3*kPi2*numNuc*
		     genie::utils::nuclear::Density(radius,A),1.0/3.0) *hbarc;
  }else{
    kf = fKFTable->FindClosestKF(tgt_pdgc, nuc_pdgc);
  }
}
