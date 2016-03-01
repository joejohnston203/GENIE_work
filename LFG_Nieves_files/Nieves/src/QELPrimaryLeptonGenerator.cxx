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
 @ Feb 18, 2016 - JJ (SD)
   A lepton is now generated and saved by the QELKinematicsGenerator.
   Updated this class to access and store that lepton.

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
  Kinematics * kine = interaction->KinePtr();

  // Check that a lepton was stored. If not, call the parent method.
  if(!(kine->KVSet(kKVSelTl) && kine->KVSet(kKVSelctl) 
       && kine->KVSet(kKVSelphikq))){
    this->PrimaryLeptonGenerator::ProcessEventRecord(evrec);
    return;
  }
  // Get the final state primary lepton energy and momentum components
  // along and perpendicular to the neutrino direction 
  double ml  = interaction->FSPrimLepton()->Mass();
  double ml2 = TMath::Power(ml,2);

  // Get the components stored by the QEL kinematics generator
  double Tl = kine->GetKV(kKVSelTl); // Used to store momentum magnitude pl
  double ctl = kine->GetKV(kKVSelctl);
  double phi = kine->GetKV(kKVSelphikq);

  double El = TMath::Sqrt(Tl*Tl + ml2);
  double stl = TMath::Sqrt(1-ctl*ctl);
  double plp = Tl * ctl;
  double pltx = Tl * stl * TMath::Cos(phi);
  double plty = Tl * stl * TMath::Sin(phi);

  TLorentzVector p4l(pltx,plty,plp,El);

  LOG("LeptonicVertex", pNOTICE)
       << "fsl @ LAB: " << utils::print::P4AsString(&p4l);

  // Figure out the Final State Lepton PDG Code
  int pdgc = interaction->FSPrimLepton()->PdgCode();

  // Create a GHepParticle and add it to the event record
  this->AddToEventRecord(evrec, pdgc, p4l);

  // Set final state lepton polarization
  this->SetPolarization(evrec);
}
//___________________________________________________________________________
