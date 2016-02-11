//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - October 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Mar 03, 2009 - CA
   Moved into the new QEL package from its previous location (EVGModules)

*/
//____________________________________________________________________________

#include <TVector3.h>
#include <TLorentzVector.h>

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Constants.h"
#include "Conventions/KineVar.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepStatus.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"
#include "QEL/QELPrimaryLeptonGenerator.h"
#include "Utils/PrintUtils.h"
#include "Utils/MathUtils.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
QELPrimaryLeptonGenerator::QELPrimaryLeptonGenerator() :
PrimaryLeptonGenerator("genie::QELPrimaryLeptonGenerator")
{

}
//___________________________________________________________________________
QELPrimaryLeptonGenerator::QELPrimaryLeptonGenerator(string config):
PrimaryLeptonGenerator("genie::QELPrimaryLeptonGenerator", config)
{

}
//___________________________________________________________________________
QELPrimaryLeptonGenerator::~QELPrimaryLeptonGenerator()
{

}
//___________________________________________________________________________
void QELPrimaryLeptonGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  //PrimaryLeptonGenerator::ProcessEventRecord(evrec);
  Interaction * interaction = evrec->Summary();

  // Boost vector for [LAB] <-> [Nucleon Rest Frame] transforms
  TVector3 beta = this->NucRestFrame2Lab(evrec);

  // Neutrino 4p
  TLorentzVector * p4v = evrec->Probe()->GetP4(); // v 4p @ LAB
  p4v->Boost(-1.*beta);                           // v 4p @ Nucleon rest frame

  // Get the final state primary lepton energy and momentum components
  // along and perpendicular to the neutrino direction 
  double Q2  = interaction->Kine().Q2(true);
  double y   = interaction->Kine().y(true);
  double Ev  = p4v->E(); 
  double ml  = interaction->FSPrimLepton()->Mass();
  double ml2 = TMath::Power(ml,2);

  LOG("LeptonicVertex", pNOTICE)
    << "Ev = " << Ev << ", Q2 = " << Q2 << ", y = " << y;

  double El  = (1-y)*Ev;
  double plp = El - 0.5*(Q2+ml2)/Ev;                          // p(//)
  double plt = TMath::Sqrt(TMath::Max(0.,El*El-plp*plp-ml2)); // p(-|)


  LOG("LeptonicVertex", pNOTICE)
      << "fsl: E = " << El << ", |p//| = " << plp << ", [pT] = " << plt;
    
  double phi;
  if(interaction->Kine().KVSet(kKVSelphil)){
    // Get stored angle of transverse components
    LOG("LeptonicVertex",pNOTICE) << "Getting phi";
    phi  = interaction->Kine().GetKV(kKVSelphil);
  }else{
    // Randomize transverse components
    LOG("LeptonicVertex",pNOTICE) << "Generating phi";
    RandomGen * rnd = RandomGen::Instance();
    phi  = 2*kPi * rnd->RndLep().Rndm();
  }
  double pltx = plt * TMath::Cos(phi);
  double plty = plt * TMath::Sin(phi);

  // Take a unit vector along the neutrino direction @ the nucleon rest frame
  TVector3 unit_nudir = p4v->Vect().Unit(); 

  // Rotate lepton momentum vector from the reference frame (x'y'z') where 
  // {z':(neutrino direction), z'x':(theta plane)} to the nucleon rest frame
  TVector3 p3l(pltx,plty,plp);
  p3l.RotateUz(unit_nudir);

  // Lepton 4-momentum in the nucleon rest frame
  TLorentzVector p4l(p3l,El);

  LOG("LeptonicVertex", pNOTICE)
       << "fsl @ NRF: " << utils::print::P4AsString(&p4l);

  // Boost final state primary lepton to the lab frame
  p4l.Boost(beta); // active Lorentz transform

  LOG("LeptonicVertex", pNOTICE)
       << "fsl @ LAB: " << utils::print::P4AsString(&p4l);

  // Figure out the Final State Lepton PDG Code
  int pdgc = interaction->FSPrimLepton()->PdgCode();

  // Create a GHepParticle and add it to the event record
  this->AddToEventRecord(evrec, pdgc, p4l);

  // Set final state lepton polarization
  this->SetPolarization(evrec);

  delete p4v;  
}
//___________________________________________________________________________
