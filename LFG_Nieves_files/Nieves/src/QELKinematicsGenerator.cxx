//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Mar 03, 2009 - CA
   Moved into the new QEL package from its previous location (EVGModules)
 @ Mar 05, 2010 - CA
   Added a temprorary SpectralFuncExperimentalCode() 
 @ Feb 06, 2013 - CA
   When the value of the differential cross-section for the selected kinematics
   is set to the event, set the corresponding KinePhaseSpace_t value too.
 @ Feb 14, 2013 - CA
   Temporarily disable the kinematical transformation that takes out the
   dipole form from the dsigma/dQ2 p.d.f.
 @ Feb 18, 2016 - JJ (SD)
   Generate a lepton in each iteration before calculating the cross
   section. Save the chosen lepton for use by the QELPrimaryLeptonGenerator
   when a Q2 value is selected.
   Catch NievesQELExceptions when calculating the cross section.
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/GBuild.h"
#include "Conventions/Controls.h"
#include "Conventions/Constants.h"
#include "Conventions/KineVar.h"
#include "Conventions/KinePhaseSpace.h"
#include "EVGCore/EVGThreadException.h"
#include "EVGCore/EventGeneratorI.h"
#include "EVGCore/RunningThreadInfo.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepFlags.h"
//#include "LlewellynSmith/NievesQELException.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGLibrary.h"
#include "QEL/QELKinematicsGenerator.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/PrintUtils.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::constants;
using namespace genie::utils;

//___________________________________________________________________________
QELKinematicsGenerator::QELKinematicsGenerator() :
KineGeneratorWithCache("genie::QELKinematicsGenerator")
{

}
//___________________________________________________________________________
QELKinematicsGenerator::QELKinematicsGenerator(string config) :
KineGeneratorWithCache("genie::QELKinematicsGenerator", config)
{

}
//___________________________________________________________________________
QELKinematicsGenerator::~QELKinematicsGenerator()
{

}
//___________________________________________________________________________
void QELKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  if(fGenerateUniformly) {
    LOG("QELKinematics", pNOTICE)
          << "Generating kinematics uniformly over the allowed phase space";
  }

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();

  //-- Get the interaction and set the 'trust' bits
  Interaction * interaction = evrec->Summary();
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  //-- Note: The kinematic generator would be using the free nucleon cross
  //   section (even for nuclear targets) so as not to double-count nuclear
  //   suppression. This assumes that a) the nuclear suppression was turned
  //   on when computing the cross sections for selecting the current event 
  //   and that b) if the event turns out to be unphysical (Pauli-blocked) 
  //   the next attempted event will be forced to QEL again.
  //   (discussion with Hugh - GENIE/NeuGEN integration workshop - 07APR2006
  interaction->SetBit(kIAssumeFreeNucleon);

  //-- Get the limits for the generated Q2
  const KPhaseSpace & kps = interaction->PhaseSpace();
  Range1D_t Q2 = kps.Limits(kKVQ2);

  if(Q2.max <=0 || Q2.min>=Q2.max) {
     LOG("QELKinematics", pWARN) << "No available phase space";
     evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
     genie::exceptions::EVGThreadException exception;
     exception.SetReason("No available phase space");
     exception.SwitchOnFastForward();
     throw exception;
  }

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  //   If the kinematics are generated uniformly over the allowed phase
  //   space the max xsec is irrelevant
  double xsec_max = (fGenerateUniformly) ? -1 : this->MaxXSec(evrec);

  //-- Try to select a valid Q2 using the rejection method

  // kinematical limits
  double Q2min  = Q2.min+kASmallNum;
  double Q2max  = Q2.max-kASmallNum;
//double QD2min = utils::kinematics::Q2toQD2(Q2min);
//double QD2max = utils::kinematics::Q2toQD2(Q2max);
  double xsec   = -1.;
  double gQ2    =  0.;

  // Store the struck nucleon position for use by the xsec method
  double radius = evrec->HitNucleon()->GetX4()->Vect().Mag();
  interaction->KinePtr()->SetKV(kKVSelRad,radius);
  
  unsigned int iter = 0;
  bool accept = false;
  while(1) {
     iter++;
     if(iter > kRjMaxIterations) {
        LOG("QELKinematics", pWARN)
          << "Couldn't select a valid Q^2 after " << iter << " iterations";
        evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
        genie::exceptions::EVGThreadException exception;
        exception.SetReason("Couldn't select kinematics");
        exception.SwitchOnFastForward();
        throw exception;
     }

     //-- Generate a Q2 value within the allowed phase space
/*
     if(fGenerateUniformly) {
         gQ2 = Q2min + (Q2max-Q2min) * rnd->RndKine().Rndm();
     } else {
         // In unweighted mode - use transform that takes out the dipole form
         double gQD2 = QD2min + (QD2max-QD2min) * rnd->RndKine().Rndm();
         gQ2  = utils::kinematics::QD2toQ2(gQD2);
     }
*/
     gQ2 = Q2min + (Q2max-Q2min) * rnd->RndKine().Rndm();
     interaction->KinePtr()->SetQ2(gQ2);
     LOG("QELKinematics", pINFO) << "Trying: Q^2 = " << gQ2;

     // Generate a lepton before calculating the cross section
     SetRunningLepton(evrec);

     //-- Computing cross section for the current kinematics
     //try{
     xsec = fXSecModel->XSec(interaction, kPSQ2fE);
     //}catch(exceptions::NievesQELException exception) {
     //LOG("QELKinematics",pINFO) << exception;
     //LOG("QELKinematics",pINFO) << "rewinding";
     //xsec = -1.0; // Do not accept
     //}
     
     //-- Decide whether to accept the current kinematics
     if(!fGenerateUniformly) {
        this->AssertXSecLimits(interaction, xsec, xsec_max);

        double t = xsec_max * rnd->RndKine().Rndm();
     //double J = kinematics::Jacobian(interaction,kPSQ2fE,kPSQD2fE);
        double J = 1.;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
        LOG("QELKinematics", pDEBUG)
            << "xsec= " << xsec << ", J= " << J << ", Rnd= " << t;
#endif
        accept = (t < J*xsec);
     } else {
        accept = (xsec>0);
     }

     //-- If the generated kinematics are accepted, finish-up module's job
     if(accept) {
        LOG("QELKinematics", pINFO) << "Selected: Q^2 = " << gQ2;

        // reset bits
        interaction->ResetBit(kISkipProcessChk);
        interaction->ResetBit(kISkipKinematicChk);
        interaction->ResetBit(kIAssumeFreeNucleon);

	// compute the rest of the kinematical variables

	// get neutrino energy at struck nucleon rest frame and the
	// struck nucleon mass (can be off the mass shell)
	const InitialState & init_state = interaction->InitState();
	double E  = init_state.ProbeE(kRfHitNucRest);
	double M = init_state.Tgt().HitNucP4().M();

	LOG("QELKinematics", pNOTICE) << "E = " << E << ", M = "<< M;

	// The hadronic inv. mass is equal to the recoil nucleon on-shell mass.
	// For QEL/Charm events it is set to be equal to the on-shell mass of
	// the generated charm baryon (Lamda_c+, Sigma_c+ or Sigma_c++)
	//
	const XclsTag & xcls = interaction->ExclTag();
	int rpdgc = 0;
	if(xcls.IsCharmEvent()) { rpdgc = xcls.CharmHadronPdg();           }
	else                    { rpdgc = interaction->RecoilNucleonPdg(); }
	assert(rpdgc);
	double gW = PDGLibrary::Instance()->Find(rpdgc)->Mass();

	LOG("QELKinematics", pNOTICE) << "Selected: W = "<< gW;

	// (W,Q2) -> (x,y)
	double gx=0, gy=0;
	kinematics::WQ2toXY(E,M,gW,gQ2,gx,gy);

        // set the cross section for the selected kinematics
        evrec->SetDiffXSec(xsec,kPSQ2fE);

        // for uniform kinematics, compute an event weight as
        // wght = (phase space volume)*(differential xsec)/(event total xsec)
        if(fGenerateUniformly) {
          double vol     = kinematics::PhaseSpaceVolume(interaction,kPSQ2fE);
          double totxsec = evrec->XSec();
          double wght    = (vol/totxsec)*xsec;
          LOG("QELKinematics", pNOTICE)  << "Kinematics wght = "<< wght;

          // apply computed weight to the current event weight
          wght *= evrec->Weight();
          LOG("QELKinematics", pNOTICE) << "Current event wght = " << wght;
          evrec->SetWeight(wght);
        }


        // lock selected kinematics & clear running values
	Kinematics * kine = interaction->KinePtr();
        kine->SetQ2(gQ2, true);
        kine->SetW (gW,  true);
        kine->Setx (gx,  true);
        kine->Sety (gy,  true);
        kine->ClearRunningValues();

	// Lock outgoing lepton
	kine->SetKV(kKVSelTl, kine->GetKV(kKVTl));
	kine->SetKV(kKVSelctl, kine->GetKV(kKVctl));
	kine->SetKV(kKVSelphikq, kine->GetKV(kKVphikq));
	
        return;
     }
  }// iterations
}
//___________________________________________________________________________
void QELKinematicsGenerator::SpectralFuncExperimentalCode(
  GHepRecord * evrec) const
{
  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();

  //-- Get the interaction and set the 'trust' bits
  Interaction * interaction = new Interaction(*evrec->Summary());
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  //-- Note: The kinematic generator would be using the free nucleon cross
  //   section (even for nuclear targets) so as not to double-count nuclear
  //   suppression. This assumes that a) the nuclear suppression was turned
  //   on when computing the cross sections for selecting the current event 
  //   and that b) if the event turns out to be unphysical (Pauli-blocked) 
  //   the next attempted event will be forced to QEL again.
  //   (discussion with Hugh - GENIE/NeuGEN integration workshop - 07APR2006
  interaction->SetBit(kIAssumeFreeNucleon);

  //-- Assume scattering off a nucleon on the mass shell (PWIA prescription)
  double Mn  = interaction->InitState().Tgt().HitNucMass(); // PDG mass, take it to be on-shell
  double pxn = interaction->InitState().Tgt().HitNucP4().Px();
  double pyn = interaction->InitState().Tgt().HitNucP4().Py();
  double pzn = interaction->InitState().Tgt().HitNucP4().Pz();
  double En  = interaction->InitState().Tgt().HitNucP4().Energy();
  double En0 = TMath::Sqrt(pxn*pxn + pyn*pyn + pzn*pzn + Mn*Mn);
  double Eb  = En0 - En;
  interaction->InitStatePtr()->TgtPtr()->HitNucP4Ptr()->SetE(En0);

  //-- Get the limits for the generated Q2
  const KPhaseSpace & kps = interaction->PhaseSpace();
  Range1D_t Q2 = kps.Limits(kKVQ2);

  if(Q2.max <=0 || Q2.min>=Q2.max) {
     LOG("QELKinematics", pWARN) << "No available phase space";
     evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
     genie::exceptions::EVGThreadException exception;
     exception.SetReason("No available phase space");
     exception.SwitchOnFastForward();
     throw exception;
  }

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  //   If the kinematics are generated uniformly over the allowed phase
  //   space the max xsec is irrelevant
//  double xsec_max = (fGenerateUniformly) ? -1 : this->MaxXSec(evrec);
  double xsec_max = this->MaxXSec(evrec);

  // get neutrino energy at struck nucleon rest frame and the
  // struck nucleon mass (can be off the mass shell)
  const InitialState & init_state = interaction->InitState();
  double E  = init_state.ProbeE(kRfHitNucRest);

  LOG("QELKinematics", pNOTICE) << "E = " << E << ", M = "<< Mn;

  //-- Try to select a valid Q2 using the rejection method

  // kinematical limits
  double Q2min  = Q2.min+kASmallNum;
  double Q2max  = Q2.max-kASmallNum;
  double xsec   = -1.;
  double gQ2    =  0.;
  double gW     =  0.;
  double gx     =  0.;
  double gy     =  0.;

  // Store the struck nucleon position for use by the xsec method
  double radius = evrec->HitNucleon()->X4()->Vect().Mag();
  interaction->KinePtr()->SetKV(kKVSelRad,radius);

  unsigned int iter = 0;
  bool accept = false;
  while(1) {
     iter++;
     if(iter > kRjMaxIterations) {
        LOG("QELKinematics", pWARN)
          << "Couldn't select a valid Q^2 after " << iter << " iterations";
        evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
        genie::exceptions::EVGThreadException exception;
        exception.SetReason("Couldn't select kinematics");
        exception.SwitchOnFastForward();
        throw exception;
     }

     //-- Generate a Q2 value within the allowed phase space
     gQ2 = Q2min + (Q2max-Q2min) * rnd->RndKine().Rndm();
     LOG("QELKinematics", pNOTICE) << "Trying: Q^2 = " << gQ2;

     // The hadronic inv. mass is equal to the recoil nucleon on-shell mass.
     // For QEL/Charm events it is set to be equal to the on-shell mass of
     // the generated charm baryon (Lamda_c+, Sigma_c+ or Sigma_c++)
     //
     const XclsTag & xcls = interaction->ExclTag();
     int rpdgc = 0;
     if(xcls.IsCharmEvent()) { rpdgc = xcls.CharmHadronPdg();           }
     else                    { rpdgc = interaction->RecoilNucleonPdg(); }
     assert(rpdgc);
     gW = PDGLibrary::Instance()->Find(rpdgc)->Mass();

     // (W,Q2) -> (x,y)
     kinematics::WQ2toXY(E,Mn,gW,gQ2,gx,gy);

     LOG("QELKinematics", pNOTICE) << "W = "<< gW;
     LOG("QELKinematics", pNOTICE) << "x = "<< gx;
     LOG("QELKinematics", pNOTICE) << "y = "<< gy;

     // v
     double gv  = gy * E;
     double gv2 = gv*gv;

     LOG("QELKinematics", pNOTICE) << "v = "<< gv;

     // v -> v~
     double gvtilde  = gv + Mn - Eb - TMath::Sqrt(Mn*Mn+pxn*pxn+pyn*pyn+pzn*pzn);
     double gvtilde2 = gvtilde*gvtilde;

     LOG("QELKinematics", pNOTICE) << "v~ = "<< gvtilde;

     // Q~^2
     double gQ2tilde = gQ2 - gv2 + gvtilde2;

     LOG("QELKinematics", pNOTICE) << "Q~^2 = "<< gQ2tilde;

     // Set updated Q2
     interaction->KinePtr()->SetQ2(gQ2tilde);

     // Generate a lepton before calculating xsec
     SetRunningLepton(evrec);

     //-- Computing cross section for the current kinematics
     //try{
     xsec = fXSecModel->XSec(interaction, kPSQ2fE);
     //}catch(exceptions::NievesQELException exception) {
     //LOG("QELKinematics",pINFO) << exception;
     //LOG("QELKinematics",pINFO) << "Setting xsec = 0.0";
     //xsec = 0.0;
     //}

     //-- Decide whether to accept the current kinematics
//     if(!fGenerateUniformly) {
        this->AssertXSecLimits(interaction, xsec, xsec_max);

        double t = xsec_max * rnd->RndKine().Rndm();
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
        LOG("QELKinematics", pDEBUG)
            << "xsec= " << xsec << ", Rnd= " << t;
#endif
        accept = (t < xsec);
//     } else {
//        accept = (xsec>0);
//     }

     //-- If the generated kinematics are accepted, finish-up module's job
     if(accept) {
        LOG("QELKinematics", pNOTICE) << "Selected: Q^2 = " << gQ2;

        // reset bits
//        interaction->ResetBit(kISkipProcessChk);
//        interaction->ResetBit(kISkipKinematicChk);
//        interaction->ResetBit(kIAssumeFreeNucleon);

        // set the cross section for the selected kinematics
        evrec->SetDiffXSec(xsec,kPSQ2fE);

        // for uniform kinematics, compute an event weight as
        // wght = (phase space volume)*(differential xsec)/(event total xsec)
//        if(fGenerateUniformly) {
//          double vol     = kinematics::PhaseSpaceVolume(interaction,kPSQ2fE);
//          double totxsec = evrec->XSec();
//          double wght    = (vol/totxsec)*xsec;
//          LOG("QELKinematics", pNOTICE)  << "Kinematics wght = "<< wght;

          // apply computed weight to the current event weight
//          wght *= evrec->Weight();
//          LOG("QELKinematics", pNOTICE) << "Current event wght = " << wght;
//          evrec->SetWeight(wght);
//        }

        // lock selected kinematics & clear running values
//        interaction->KinePtr()->SetQ2(gQ2, true);
//        interaction->KinePtr()->SetW (gW,  true);
//        interaction->KinePtr()->Setx (gx,  true);
//        interaction->KinePtr()->Sety (gy,  true);
//        interaction->KinePtr()->ClearRunningValues();

	Kinematics * kine = evrec->Summary()->KinePtr();
        kine->SetQ2(gQ2, true);
        kine->SetW (gW,  true);
        kine->Setx (gx,  true);
        kine->Sety (gy,  true);
        kine->ClearRunningValues();

	// Lock outgoing lepton
	kine->SetKV(kKVSelTl, kine->GetKV(kKVTl));
	kine->SetKV(kKVSelctl, kine->GetKV(kKVctl));
	kine->SetKV(kKVSelphikq, kine->GetKV(kKVphikq));
	delete interaction;

        return;
     }
  }// iterations
}
//___________________________________________________________________________
void QELKinematicsGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void QELKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void QELKinematicsGenerator::LoadConfig(void)
{
// Load sub-algorithms and config data to reduce the number of registry
// lookups

  //-- Safety factor for the maximum differential cross section
  fSafetyFactor = fConfig->GetDoubleDef("MaxXSec-SafetyFactor", 1.25);

  //-- Minimum energy for which max xsec would be cached, forcing explicit
  //   calculation for lower eneries
  fEMin = fConfig->GetDoubleDef("Cache-MinEnergy", 1.00);

  //-- Maximum allowed fractional cross section deviation from maxim cross
  //   section used in rejection method
  fMaxXSecDiffTolerance = 
               fConfig->GetDoubleDef("MaxXSec-DiffTolerance",999999.);
  assert(fMaxXSecDiffTolerance>=0);

  //-- Generate kinematics uniformly over allowed phase space and compute
  //   an event weight?
  fGenerateUniformly = fConfig->GetBoolDef("UniformOverPhaseSpace", false);
}
//____________________________________________________________________________
double QELKinematicsGenerator::ComputeMaxXSec(
                                       const Interaction * interaction) const
{
// Computes the maximum differential cross section in the requested phase
// space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
// method and the value is cached at a circular cache branch for retrieval
// during subsequent event generation.
// The computed max differential cross section does not need to be the exact
// maximum. The number used in the rejection method will be scaled up by a
// safety factor. But it needs to be fast - do not use a very small dQ2 step.

  double max_xsec = 0.0;

  const KPhaseSpace & kps = interaction->PhaseSpace();
  Range1D_t rQ2 = kps.Limits(kKVQ2);
  if(rQ2.min <=0 || rQ2.max <= rQ2.min) return 0.;

  const double logQ2min = TMath::Log(rQ2.min + kASmallNum);
  const double logQ2max = TMath::Log(rQ2.max - kASmallNum);

  const int N  = 15;
  const int Nb = 10;

  double dlogQ2   = (logQ2max - logQ2min) /(N-1);
  double xseclast = -1;
  bool   increasing;

  for(int i=0; i<N; i++) {
     double Q2 = TMath::Exp(logQ2min + i * dlogQ2);
     interaction->KinePtr()->SetQ2(Q2);
     double xsec;
     //-- Computing cross section for the current kinematics
     //try{
     xsec = fXSecModel->XSec(interaction, kPSQ2fE);
     //}catch(exceptions::NievesQELException exception) {
     //LOG("QELKinematics",pINFO) << exception;
     //LOG("QELKinematics",pINFO) << "Setting xsec = 0.0";
     //xsec = 0.0;
     //}
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("QELKinematics", pDEBUG)  << "xsec(Q2= " << Q2 << ") = " << xsec;
#endif
     max_xsec = TMath::Max(xsec, max_xsec);
     increasing = xsec-xseclast>=0;
     xseclast   = xsec;

     // once the cross section stops increasing, I reduce the step size and
     // step backwards a little bit to handle cases that the max cross section
     // is grossly underestimated (very peaky distribution & large step)
     if(!increasing) {
       dlogQ2/=(Nb+1);
       for(int ib=0; ib<Nb; ib++) {
	 Q2 = TMath::Exp(TMath::Log(Q2) - dlogQ2);
         if(Q2 < rQ2.min) continue;
         interaction->KinePtr()->SetQ2(Q2);
	 //try{
	 xsec = fXSecModel->XSec(interaction, kPSQ2fE);
	 //}catch(exceptions::NievesQELException exception) {
	 //LOG("QELKinematics",pINFO) << exception;
	 //LOG("QELKinematics",pINFO) << "Setting xsec = 0.0";
	 //xsec = 0.0;
	 //}
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
         LOG("QELKinematics", pDEBUG)  << "xsec(Q2= " << Q2 << ") = " << xsec;
#endif
         max_xsec = TMath::Max(xsec, max_xsec);
       }
       break;
     }
  }//Q^2

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy
  max_xsec *= fSafetyFactor;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  SLOG("QELKinematics", pDEBUG) << interaction->AsString();
  SLOG("QELKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("QELKinematics", pDEBUG) << "Computed using alg = " << *fXSecModel;
#endif

  return max_xsec;
}
//___________________________________________________________________________
// Generate a lepton and store it by using the KineVar construct with
// the Kinematics class (variables kKVTl, kKVctl, kKVphikq).
// Stored lepton is in the LAB FRAME.
void QELKinematicsGenerator::SetRunningLepton(GHepRecord * evrec) const{
  Interaction * interaction = evrec->Summary();

  // Velocity for an active Lorentz transform taking the final state primary
  // lepton from the [nucleon rest frame] --> [LAB]
  const InitialState & init_state = interaction->InitState();
  const TLorentzVector & pnuc4 = init_state.Tgt().HitNucP4(); //[@LAB]
  TVector3 beta = pnuc4.BoostVector();


  // Neutrino 4p
  TLorentzVector * p4v = evrec->Probe()->GetP4(); // v 4p @ LAB
  p4v->Boost(-1.*beta);                           // v 4p @ Nucleon rest frame

  // Look-up selected kinematics & other needed kinematical params
  double Q2  = interaction->Kine().Q2(false);

  // get neutrino energy at struck nucleon rest frame and the
  // struck nucleon mass (can be off the mass shell)
  double E  = init_state.ProbeE(kRfHitNucRest);
  double M = init_state.Tgt().HitNucP4().M();

  const XclsTag & xcls = interaction->ExclTag();
  int rpdgc = 0;
  if(xcls.IsCharmEvent()) { rpdgc = xcls.CharmHadronPdg();           }
  else                    { rpdgc = interaction->RecoilNucleonPdg(); }
  assert(rpdgc);
  double W = PDGLibrary::Instance()->Find(rpdgc)->Mass();
  // (W,Q2) -> (x,y)
  double x=0, y=0;
  kinematics::WQ2toXY(E,M,W,Q2,x,y);
  double Ev  = p4v->E(); 
  double ml  = interaction->FSPrimLepton()->Mass();
  double ml2 = TMath::Power(ml,2);

  LOG("QELKinematics", pINFO)
             << "Ev = " << Ev << ", Q2 = " << Q2 << ", y = " << y;

  // Compute the final state primary lepton energy and momentum components
  // along and perpendicular the neutrino direction 
  double El = (1-y)*Ev;
  double plp = El - 0.5*(Q2+ml2)/Ev;                          // p(//)
  double plt = TMath::Sqrt(TMath::Max(0.,El*El-plp*plp-ml2)); // p(-|)

  LOG("QELKinematics", pINFO)
        << "trying fsl: E = " << El << ", |p//| = " << plp << ", [pT] = " << plt;

  // Randomize transverse components
  RandomGen * rnd = RandomGen::Instance();
  double phi  = 2*kPi * rnd->RndLep().Rndm();
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

  LOG("QELKinematics", pINFO)
       << "trying fsl @ NRF: " << utils::print::P4AsString(&p4l);

  // Boost final state primary lepton to the lab frame
  p4l.Boost(beta); // active Lorentz transform

  LOG("QELKinematics", pINFO)
       << "trying fsl @ LAB: " << utils::print::P4AsString(&p4l);

  // Store the outgoing lepton components in the lab frame
  // Assume El^2 = Tl^2 + ml^2
  double Tl = p4l.Vect().Mag(); // used to store 3-vector magnitude
  LOG("QELKinematics",pDEBUG) << "ml = " << ml << ", Tl = " << Tl;
  double ctl = p4l.CosTheta(); // cos(theta) in the lab frame
  phi = p4l.Phi(); // get phi in the lab fram
  interaction->KinePtr()->SetKV(kKVTl,Tl);
  interaction->KinePtr()->SetKV(kKVctl,ctl);
  interaction->KinePtr()->SetKV(kKVphikq,phi);

  delete p4v;
}
//___________________________________________________________________________

