//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Joe Johnston, University of Pittsburgh

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <iostream>
#include <fstream>
#include <exception>
#include <complex>

#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <ctype.h>

#include "Algorithm/AlgConfigPool.h"
#include "Base/XSecIntegratorI.h"
#include "Base/QELFormFactors.h"
#include "Base/QELFormFactorsModelI.h"
#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Conventions/KineVar.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "LlewellynSmith/NievesQELCCPXSec.h"
#include "LlewellynSmith/NievesQELException.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGUtils.h"
#include "QEL/QELKinematicsGenerator.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/NuclearUtils.h"
#include "Utils/PrintUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils;
using namespace genie::units;

//____________________________________________________________________________
NievesQELCCPXSec::NievesQELCCPXSec() :
XSecAlgorithmI("genie::NievesQELCCPXSec")
{

}
//____________________________________________________________________________
NievesQELCCPXSec::NievesQELCCPXSec(string config) :
XSecAlgorithmI("genie::NievesQELCCPXSec", config)
{

}
//____________________________________________________________________________
NievesQELCCPXSec::~NievesQELCCPXSec()
{

 }
//____________________________________________________________________________
double NievesQELCCPXSec::XSec(const Interaction * interaction,
		       KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  // Get kinematics and init-state parameters
  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();
  const Target & target = init_state.Tgt();

  double E  = init_state.ProbeE(kRfLab);//Incoming Neutrino Energy, Lab
  double ml = interaction->FSPrimLepton()->Mass();
  double M  = target.HitNucMass();
  double q2 = kinematics.q2();

  if(q2>=0.0){
    LOG("Nieves", pWARN) << "q2>=0";
    exceptions::NievesQELException exception;
    exception.SetReason("Invalid q, q2>=0.0");
    exception.SetUnphysicalQ2(true);
    throw exception;
  }

  LOG("Nieves", pDEBUG) << "q2 = " << q2;

  // Calculate auxiliary parameters
  double ml2     = TMath::Power(ml,    2);
  double M2      = TMath::Power(M,     2);

  double s = (2*E+M)*M;
  double num = TMath::Power(s-M2,2);
  double Gfactor = kGF2 * fCos8c2 / (8*kPi*num);

  double coulombFactor = 1.0;

  // Compute nucleon differential cross section

  // Lepton momenta in coordinate system with the nuetrino in the z dir
  // (lab frame)

 // Initial neutrino
  TVector3 k3Vec(0,0,E);

  // Generate outgoing lepton, and store for use by primary lepton generator
  TLorentzVector * probe = new TLorentzVector(k3Vec,E);
  try{
    SetRunningOutgoingLepton(interaction, probe);
    delete probe;
  }catch(exceptions::NievesQELException exception) {
    // Invalid kinematics
    // delete the probe, then pass the exception to the calling class
    delete probe;
    throw exception;
  }

  double El = ml + kinematics.GetKV(kKVTl);
  double ctl = kinematics.GetKV(kKVctl);
  double stl = TMath::Sqrt(1-ctl*ctl); // Sin of polar angle
  double phi = kinematics.GetKV(kKVphikq);
  double pl = TMath::Sqrt(El*El-ml*ml);

  // Coulomb Effects
  if(fCoulomb){
    if(kinematics.KVSet(kKVSelRad))
      r = kinematics.GetKV(kKVSelRad);
    else
      r = 0.0;
    // Coulomb potential
    double Vc = vcr(& target, double r);

    // Outgoing lepton energy and momentum including coulomb potential
    int sign = (pdg::IsNeutrino(init_state.ProbePdg())) ? 1 : -1;
    double ElLocal = El - sign*Vc;
    if(ElLocal - ml <= 0.0){
      LOG("Nieves",pINFO) << "Event should be rejected. Coulomb effects "
			    << "push kinematics below threshold";
      exceptions::NievesQELException exception;
      exception.SetReason("Outgoing lepton energy below threshold after coulomb effects");
      exception.SetUnphysicalQ2(true);
      throw exception;
    }
    double plLocal = TMath::Sqrt(ElLocal*ElLocal-ml2);

    // Correction factor
    coulombFactor= plLocal*ElLocal/pl/El;
    
    // Could use corrected outgoing momentum to calculate cross section
    //pl = plLocal;
    //El = ElLocal;
  }
  
  TVector3 kPrime3Vec(pl*stl*TMath::Cos(phi),pl*stl*TMath::Sin(phi),pl*ctl);

  TVector3 q3Vec = k3Vec-kPrime3Vec;

  // Rotate k and kPrime so q is in the z direction
  TVector3 rot = (k3Vec.Cross(kPrime3Vec)).Unit(); // Vector to rotate about
  double angle = k3Vec.Angle(q3Vec); // Angle between the z direction and q

  // Rotate if the rotation vector is not 0
  if(! rot.Mag() < 0.0001){
    k3Vec.Rotate(angle,rot);
    kPrime3Vec.Rotate(angle,rot);
    q3Vec = k3Vec-kPrime3Vec;
  }
  //else{
    //    LOG("Nieves",pDEBUG) << "Rotation angle is 0, qx = " << q3Vec.x()
    //		<< ", qy = " << q3Vec.y()
    //		<< ", qz = " << q3Vec.z();
  //}

  TLorentzVector kVec(k3Vec,E);
  TLorentzVector kPrimeVec(kPrime3Vec,El);
  TLorentzVector qVec = kVec - kPrimeVec;

  // Substract the binding energy per nucleon from q0 to account for the
  // energy needed to remove the outgoing particle from the nucleus
  //double Qvalue = nuclear::BindEnergyPerNucleon(target);
  //LOG("Nieves",pDEBUG) << "Binding Energy Per Nucleon = " << Qvalue;
  double Qvalue = 0.016827;
  qVec.SetE(qVec.E() - Qvalue);

  // Check that the energy tranfer q0 is greater than 0, or else the
  // following equations do not apply. (Note also that the event would
  // be Pauli blocked when suppression is on)
  /*if(qVec.E()<=0){
    LOG("Nieves", pINFO) << "q0<=0.0";
    exceptions::NievesQELException exception;
    exception.SetReason("Invalid q, q0<=0.0");

    // Assume this nucleon and q2 are certain to give an unphysical event
    // Maybe they wouldn't, but is is unlikely that if this has happened
    // that a different phi would give a non-negative q0
    exception.SetUnphysicalQ2(true);
    throw exception;
    }*/

  ofstream miscstream;
  if(fPrintData){
    miscstream.open("initVals.txt", std::ios_base::app);
    miscstream << -q2 << "\t" << El << "\t" << ctl << "\t" << -qVec.Mag2() << "\n";
    miscstream.close();
  }

  // pInit is the initial nucleon momentum in the coordinates where 
  // q is in the z direction.
  TLorentzVector pInit = target.HitNucP4();

  double LmunuAnumuResult = LmunuAnumu(interaction,pInit,kVec,kPrimeVec);

  double xsec = Gfactor*coulombFactor*LmunuAnumuResult;

  LOG("Nieves",pDEBUG) << "RPA=" << fRPA 
		       << ", Coulomb=" << fCoulomb 
		       << ", Nieves Suppression=" << fNievesSuppression
		       << ", q2 = " << q2 << ", xsec = " << xsec;

  ofstream uhoh;
  if( isnan(xsec)){
    LOG("Nieves",pWARN) << "xsec is nan";
    uhoh.open("uhoh.txt", std::ios_base::app);
    uhoh << q2 << "\t" << xsec << "\t" << qVec.E() << "\t" 
	 << qVec.Vect().Mag() << "\t" << E << "\n";
    uhoh.close();
    exceptions::NievesQELException exception;
    exception.SetReason("xsec is nan");
    throw exception;
  }
  
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Nieves", pDEBUG)
     << "dXSec[QEL]/dQ2 [FreeN](E = "<< E << ", Q2 = "<< -q2 << ") = "<< xsec;
#endif

  //----- The algorithm computes dxsec/dQ2
  //      Check whether variable tranformation is needed
  if(kps!=kPSQ2fE) {
    double J = utils::kinematics::Jacobian(interaction,kPSQ2fE,kps);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("Nieves", pDEBUG)
     << "Jacobian for transformation to: " 
                  << KinePhaseSpace::AsString(kps) << ", J = " << J;
#endif
    xsec *= J;
  }
  //----- if requested return the free nucleon xsec even for input nuclear tgt
  //      the xsec will still include RPA corrections
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  // Suppression commented out here because if fNievesSuppression is false,
  // no suppression is used
  /*
  //----- compute nuclear suppression factor
  //      (R(Q2) is adapted from NeuGEN - see comments therein)
  double R;
  if(fNievesSuppression)
    R = 1.0; //Suppresssion factor accounted for by theta functions
  else
    R = nuclear::NuclQELXSecSuppression("Default", 0.5, interaction);
  */

  //----- number of scattering centers in the target
  int nucpdgc = target.HitNucPdg();
  int NNucl = (pdg::IsProton(nucpdgc)) ? target.Z() : target.N(); 
  /*
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Nieves", pDEBUG) 
       << "Nuclear suppression factor R(Q2) = " << R << ", NNucl = " << NNucl;
#endif
xsec *= (R*NNucl); // nuclear xsec*/
  xsec *= NNucl;

  return xsec;
}//____________________________________________________________________________
double NievesQELCCPXSec::Integral(const Interaction * interaction) const
{
  LOG("Nieves",pDEBUG) << "genie::NievesQELCCPXSec::Integral() called";
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  LOG("Nieves",pDEBUG) << "Integral xsec = " << xsec << "\n";
  return xsec;
}
//____________________________________________________________________________
bool NievesQELCCPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();

  if(!proc_info.IsQuasiElastic()) return false;

  int  nuc = init_state.Tgt().HitNucPdg();
  int  nu  = init_state.ProbePdg();

  bool isP   = pdg::IsProton(nuc);
  bool isN   = pdg::IsNeutron(nuc);
  bool isnu  = pdg::IsNeutrino(nu);
  bool isnub = pdg::IsAntiNeutrino(nu);

  bool prcok = proc_info.IsWeakCC() && ((isP&&isnub) || (isN&&isnu));
  if(!prcok) return false;

  return true;
}
//____________________________________________________________________________
void NievesQELCCPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NievesQELCCPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NievesQELCCPXSec::LoadConfig(void)
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();
  
  double thc = fConfig->GetDoubleDef(
                              "CabbiboAngle", gc->GetDouble("CabbiboAngle"));
  fCos8c2 = TMath::Power(TMath::Cos(thc), 2);

  // variables 
  double intervalFractions[10] = {.9931285991,.9639719272,.9122344282,
  				  .8391169718,.7463319064,.6360536807,
  		       .5108670019,.3737060887,.2277858511,.0765265211};
  fIntervalFractions.assign(intervalFractions,intervalFractions + 10);

  double w[10] = {.0176140071,.0406014298,.0626720483,.0832767415,.1019301198,
		  .1181945319,.1316886384,.1420961093,.1491729864,.1527533871};
  fW.assign(w,w+10);
  // hbarc for unit conversion, GeV*fm
  fhbarc = kLightSpeed*kPlankConstant/fermi;

   // load QEL form factors model
  fFormFactorsModel = dynamic_cast<const QELFormFactorsModelI *> (
                                             this->SubAlg("FormFactorsAlg"));
  assert(fFormFactorsModel);
  fFormFactors.SetModel(fFormFactorsModel); // <-- attach algorithm

   // load XSec Integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);

  // Load settings for RPA and Coulomb effects

  // RPA corrections will not effect a free nucleon
  fRPA = fConfig->GetBoolDef("RPA",true);
  
  // Coulomb Correction- adds a correction factor, and alters outgoing lepton 
  // 3-momentum magnitude (but not direction)
  // Correction only becomes sizeable near threshold and/or for heavy nuclei
  fCoulomb = fConfig->GetBoolDef("Coulomb",true);

  LOG("Nieves",pNOTICE) << "RPA=" << fRPA << ", useCoulomb=" << fCoulomb;

  fPrintData = fConfig->GetBoolDef("PrintData",false);

  // Create a nuclear model object to check the model type
  RgKey nuclkey = "NuclearModel";
  RgAlg nuclalg = gc->GetAlg(nuclkey);
  AlgFactory * algf = AlgFactory::Instance();
  const NuclearModelI* nuclModel = 
    dynamic_cast<const NuclearModelI*>(
			     algf->GetAlgorithm(nuclalg.name,nuclalg.config));
  // Check if the model is a local Fermi gas
  fLFG = (nuclModel && nuclModel->ModelType(Target()) == kNucmLocalFermiGas);
  
  if(!fLFG){
    // get the Fermi momentum table for relativistic Fermi gas
    fKFTableName = fConfig->GetStringDef ("FermiMomentumTable",
					  gc->GetString("FermiMomentumTable"));
    fKFTable = 0;
    
    FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
    fKFTable = kftp->GetTable(fKFTableName);
    assert(fKFTable);
  }
}
//___________________________________________________________________________
void NievesQELCCPXSec::SetRunningOutgoingLepton(
	  const Interaction * interaction, TLorentzVector * probe) const {
  // -----------------------------------------------------------------------
  // Uses Q2, initial neutrino 4-momentum, and the initial nucleon
  // 4-momentum. Calculates outgoing lepton momentum by boosting the initial 
  // conditions to the nucleon rest frame, calculating the outgoing lepton 
  // energy and polar angle (which are fixed given the initial conditions),
  // and randomizing the azimuthal angle. This 4-momentum is then boosted 
  // back to the lab frame and stored to by used by the xsec model
  // -----------------------------------------------------------------------

  // Boost vector for [LAB] <-> [Nucleon Rest Frame] transforms
  const InitialState & init_state = interaction->InitState();

  const TLorentzVector & pnuc4 = init_state.Tgt().HitNucP4(); //[@LAB]
  
  TVector3 beta = pnuc4.BoostVector();

  // Neutrino 4p
  //TLorentzVector * probe = evrec->Probe()->GetP4(); // v 4p @ LAB
  probe->Boost(-1.*beta);                           // v 4p @ Nucleon rest frame


  // get neutrino energy at struck nucleon rest frame and the
  // struck nucleon mass (can be off the mass shell)
  double Ev  = init_state.ProbeE(kRfHitNucRest);
  double M = init_state.Tgt().HitNucP4().M();
  
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
  
  //LOG("Nieves", pNOTICE) << "Running: W = "<< gW;

  double gQ2  = interaction->Kine().Q2();

  //if(fPrintData){
  //ofstream nucstream;
  //nucstream.open("nucleon_lab.txt", std::ios_base::app);
  //nucstream << gQ2 << "\t" << pnuc4.E() << "\t" << pnuc4.Px() << "\t"
  //	    << pnuc4.Py() << "\t" << pnuc4.Pz() << "\t"
  //    << pnuc4.Vect().Mag() << "\n";
  //nucstream.close();
  //}

  // (W,Q2) -> (x,y)
  double gx=0, gy=0;
  kinematics::WQ2toXY(Ev,M,gW,gQ2,gx,gy);

  // Store kinematics variables to be stored if this lepton is selected
  interaction->KinePtr()->SetKV(kKVW, gW);
  interaction->KinePtr()->SetKV(kKVx, gx);
  interaction->KinePtr()->SetKV(kKVy, gy);

  //lepstream << gQ2 << "\t" <<  Ev << "\t" << M << "\t" << gW <<"\t" << gy;

  // Look-up selected kinematics & other needed kinematical params
  double ml  = interaction->FSPrimLepton()->Mass();
  double ml2 = TMath::Power(ml,2);

  //LOG("Nieves", pINFO)
  //           << "Ev = " << Ev << ", Q2 = " << Q2 << ", y = " << y;

  // Compute the final state primary lepton energy and momentum components
  // along and perpendicular the neutrino direction 
  double El  = (1-gy)*Ev;
  if(El < ml){
    LOG("Nieves", pNOTICE) << "El < ml";
    exceptions::NievesQELException exception;
    exception.SetReason("Q2 and initial nucleon give El < ml");
    exception.SetUnphysicalQ2(true);
    throw exception;
  }
  double plp = El - 0.5*(gQ2+ml2)/Ev;                          // p(//)
  double plt = TMath::Sqrt(TMath::Max(0.,El*El-plp*plp-ml2)); // p(-|)

  //LOG("Nieves", pINFO)
  //      << "fsl: E = " << El << ", |p//| = " << plp << ", [pT] = " << plt;

  // Randomize transverse components
  RandomGen * rnd = RandomGen::Instance();
  double phi  = 2*kPi * rnd->RndLep().Rndm();
  double pltx = plt * TMath::Cos(phi);
  double plty = plt * TMath::Sin(phi);

  // Take a unit vector along the neutrino direction @ the nucleon rest frame
  TVector3 unit_nudir = probe->Vect().Unit(); 

  // Rotate lepton momentum vector from the reference frame (x'y'z') where 
  // {z':(neutrino direction), z'x':(theta plane)} to the nucleon rest frame
  TVector3 p3l(pltx,plty,plp);
  p3l.RotateUz(unit_nudir);

  // Lepton 4-momentum in the nucleon rest frame
  TLorentzVector p4l(p3l,El);

  TLorentzVector qNRF = *probe-p4l;

  // Boost final state primary lepton to the lab frame
  p4l.Boost(beta); // active Lorentz transform

  LOG("Nieves", pINFO)
       << "(running) fsl @ LAB: " << utils::print::P4AsString(&p4l);

  double Elab  = init_state.ProbeE(kRfLab);
  TVector3 k3Vec(0,0,Elab);    
 
  TLorentzVector * probelab = new TLorentzVector(k3Vec,Elab);
  
  TLorentzVector qLab = *probelab - p4l;

  ofstream lepstream;
  if(fPrintData){
    lepstream.open("QELKinGenLep.txt", std::ios_base::app);
    
    lepstream << -qNRF.Mag2() << "\t" << -qLab.Mag2() << "\t" 
	      << qNRF.E() << "\t" << qLab.E() << "\t"
	      << qNRF.Px() << "\t" << qLab.Px() << "\t"
	      << qNRF.Py() << "\t" << qLab.Py() << "\t"
	      << qNRF.Pz() << "\t" << qLab.Pz() << "\t"
	      << qNRF.Vect().Mag() << "\t" << qLab.Vect().Mag() << "\t"
	      << pnuc4.E() << "\t" << pnuc4.Vect().Mag() << "\t"
	      << pnuc4.Mag2() << "\t" << beta.Mag() << "\t";
    qNRF.Boost(beta);
    lepstream << qNRF.E() << "\t" << qNRF.Vect().Mag()/qNRF.E() << "\t";
  }
    
  // Store lab 4-momentum using SetKV
  double Tl = p4l.E() - ml;
  if(Tl < 0.0){
    LOG("Nieves", pNOTICE) << "Tl < 0 in lab frame";
    exceptions::NievesQELException exception;
    exception.SetReason("Q2 and initial nucleon give El < ml");
    exception.SetUnphysicalQ2(true);
    delete probelab;
    throw exception;
  }
  interaction->KinePtr()->SetKV(kKVTl, p4l.E() - ml);
  interaction->KinePtr()->SetKV(kKVctl ,p4l.CosTheta());
  interaction->KinePtr()->SetKV(kKVphikq, p4l.Phi());

  if(fPrintData){
    lepstream << interaction->KinePtr()->GetKV(kKVTl) << "\t";
    lepstream << interaction->KinePtr()->GetKV(kKVctl) << "\t";
    lepstream << interaction->KinePtr()->GetKV(kKVphikq) << "\n";
    lepstream.close();
  }

  delete probelab;

}
//___________________________________________________________________________
void NievesQELCCPXSec::CNCTCLimUcalc(const Target * target,
				     radius r,
				     TLorentzVector qVec,double q2,
				     bool is_neutrino,
				     double & CN, double & CT, double & CL,
				     double & imaginaryU,
				     double & t0,
				     double & r00) const
{
  if(target->IsNucleus()){
    double M = target->HitNucMass();
    double dq = qVec.Vect().Mag();
    double dq2 = TMath::Power(dq,2);
    //Terms for polarization coefficients CN,CT, and CL
    double hbarc2 = TMath::Power(fhbarc,2);
    double c0 = 0.380/fhbarc;//Constant for CN in natural units

    //Density gives the nuclear density, normalized to 1
    //Input radius r must be in fm
    int A = target->A();
    int Z = target->Z();
    int N = target->N();
    double rhop = nuclear::Density(r,A)*Z;
    double rhon = nuclear::Density(r,A)*N;
    double rho = rhop + rhon;
    double rho0 = A*nuclear::Density(0,A);
    
    double fPrime = 0.33*rho/rho0+0.45*(1-rho/rho0);
    
    // Get Fermi momenta
    double kF1, kF2;
    if(fLFG){
      if(pdg::IsProton(target->HitNucPdg())){
	kF1 = TMath::Power(3*kPi2*rhop, 1.0/3.0) *fhbarc;
	kF2 = TMath::Power(3*kPi2*rhon, 1.0/3.0) *fhbarc;
      }else{
	kF1 = TMath::Power(3*kPi2*rhon, 1.0/3.0) *fhbarc;
	kF2 = TMath::Power(3*kPi2*rhop, 1.0/3.0) *fhbarc;
      }
    }else{
      int tgt_pdgc = target.Pdg();
      if(pdg::IsProton(target->HitNucPdg())){
	kF1 = fKFTable->FindClosestKF(tgt_pdgc, kPdgProton);
	kF2 = fKFTable->FindClosestKF(tgt_pdgc, kPdgNeutron);
      }else{
	kF1 = fKFTable->FindClosestKF(tgt_pdgc, kPdgNeutron);
	kF2 = fKFTable->FindClosestKF(tgt_pdgc, kPdgProton);	
      }
    }
    //LOG("Nieves",pDEBUG) << "r=" << r << ",kF1=" << kF1 << ",kF2=" << kF2;
    double kF = kF1 + kF2;

    std::complex<double> imU(relLindhardIm(qVec.E(),dq,kF1,kF2,
					   M,is_neutrino,t0,r00));

    imaginaryU = imag(imU);
    //LOG("Nieves",pDEBUG) << "imU = " << imaginaryU;
    std::complex<double> relLin(0,0),udel(0,0);

    /* LOG("Nieves",pDEBUG) << "q0 = " <<qVec.E() << ", dq = " << dq
			 << ", kF = " << kF << ", M = " << M 
			 << ", is_neu = " << is_neutrino
			 << ", imU = " << imU;*/

    if(abs(imU)>pow(10.0,-10)){
      relLin = relLindhard(qVec.E(),dq,kF,M,is_neutrino,imU);
      udel = deltaLindhard(qVec.E(),dq,rho,kF);
    }
    std::complex<double> relLinTot(relLin + udel);


  /* CRho = 2
     DeltaRho = 2500 MeV, (2.5 GeV)^2 = 6.25 GeV^2
     mRho = 770 MeV, (0.770 GeV)^2 = 0.5929 GeV^2
     g' = 0.63 */
    double Vt = 0.08*4*kPi/kPionMass2 *  
      (2* TMath::Power((6.25-0.5929)/(6.25-q2),2)*dq2/(q2-0.5929) + 0.63);
  /* f^2/4/Pi = 0.08
     DeltaSubPi = 1200 MeV, (1.2 GeV)^2 = 1.44 GeV^2
     g' = 0.63 */
    double Vl = 0.08*4*kPi/kPionMass2 * 
      (TMath::Power((1.44-kPionMass2)/(1.44-q2),2)*dq2/(q2-kPionMass2)+0.63);
    
    CN = 1.0/TMath::Power(abs(1.0-c0*fPrime*relLin/hbarc2),2);

    CT = 1.0/TMath::Power(abs(1.0-relLinTot*Vt),2);
    CL = 1.0/TMath::Power(abs(1.0-relLinTot*Vl),2);
    //LOG("Nieves",pDEBUG) <<"CN = " << CN <<",CT = " << CT << ",CL = " << CL;
    }else{
    //Polarization Coefficients: all equal to 1.0 for free nucleon
    CN = 1.0;
    CT = 1.0;
    CL = 1.0;
    imaginaryU = 0.0;
  }
}
//____________________________________________________________________________
// Gives the imaginary part of the relativistic lindhard function in GeV^2
// and sets the values of t0 and r00
std::complex<double> NievesQELCCPXSec::relLindhardIm(double q0, double dq, 
						     double kFn, double kFp, 
						     double M,
						     bool isNeutrino,
						     double & t0,
						     double & r00) const
{
  double M2 = TMath::Power(M,2);
  double EF1,EF2;
  if(isNeutrino){
    EF1 = TMath::Sqrt(M2+TMath::Power(kFn,2)); //EFn
    EF2 = TMath::Sqrt(M2+TMath::Power(kFp,2)); //EFp
  }else{
    EF1 = TMath::Sqrt(M2+TMath::Power(kFp,2)); //EFp
    EF2 = TMath::Sqrt(M2+TMath::Power(kFn,2)); //EFn
  }

  double q2 = TMath::Power(q0,2) - TMath::Power(dq,2);
  double a = (-q0+dq*TMath::Sqrt(1-4.0*M2/q2))/2.0;
  double epsRP = TMath::Max(TMath::Max(M,EF2-q0),a);

  // Other theta functions for q are handled by nuclear suppression
  // That is, q0>0 and -q2>0 are always handled, and q0>EF2-EF1 is
  // handled if pauli blocking is on, because otherwise the final
  // nucleon would be below the fermi sea
  // This theta function should only be used if suppression is on
  //if(fNievesSuppression && !interaction->TestBit(kIAssumeFreeNucleon )
  //&& !EF1-epsRP<0){
  //LOG("Nieves", pINFO) << "Average value of E(p) above Fermi sea";
  //exceptions::NievesQELException exception;
  //exception.SetReason("Average nucleon energy above Fermi sea");
  //throw exception;
  //}else{
    t0 = 0.5*(EF1+epsRP);
    r00 = (TMath::Power(EF1,2)+TMath::Power(epsRP,2)+EF1*epsRP)/3.0;
    std::complex<double> result(0.0,-M2/2.0/kPi/dq*(EF1-epsRP));
    return result;
    //}
}
//____________________________________________________________________________
//Following obtained from fortran code by J Nieves, which contained the following comment:
/*
 NUCLEON relativistic Lindhard Function
 Same normalization as ULIN
 Real part
 taken from Eur.Phys.J.A25:299-318,2005 (Barbaro et al)
 Eq. 61

 Im. part: Juan. 

 INPUT: Real*8 
  q0:Energy   [fm]
  qm: modulus 3mom [fm]
  kf: Fermi mom [fm]

 OUTPUT: Complex*16 [fm]

 USES: ruLinRelX, relLindhardIm
 */
//Takes inputs in GeV (with imU in GeV^2), and gives output in GeV^2
std::complex<double> NievesQELCCPXSec::relLindhard(double q0gev, 
		        double dqgev, double kFgev, double M, 
			bool isNeutrino, std::complex<double> relLindIm) const
{
  double q0 = q0gev/fhbarc;
  double qm = dqgev/fhbarc;
  double kf = kFgev/fhbarc;
  double m = kFgev/fhbarc;

  if(q0>qm){
    LOG("Nieves", pWARN) << "relLindhard() failed";
    exceptions::NievesQELException exception;
    exception.SetReason("relLindhard not valid for q2 > 0");
    throw exception;
  }

  std::complex<double> RealLinRel(ruLinRelX(q0,qm,kf,m)+ruLinRelX(-q0,qm,kf,m));
  //Units of GeV^2
  return(RealLinRel*TMath::Power(fhbarc,2) + 2.0*relLindIm);
}
//____________________________________________________________________________
//Inputs assumed to be in natural units
std::complex<double> NievesQELCCPXSec::ruLinRelX(double q0, double qm, 
						 double kf, double m) const
{
  double q02 = TMath::Power(q0, 2);
  double qm2 = TMath::Power(qm, 2);
  double kf2 = TMath::Power(kf, 2);
  double m2  = TMath::Power(m,  2);
  double m4  = TMath::Power(m,  4);

  double ef = TMath::Sqrt(m2+kf2);
  double q2 = q02-qm2;
  double q4 = TMath::Power(q2,2);
  double ds = TMath::Sqrt(1.0-4.0*m2/q2);
  double L1 = log((kf+ef)/m);
  std::complex<double> uL2(
       TMath::Log(TMath::Abs(
		    (ef + q0 - TMath::Sqrt(m2+TMath::Power(kf-qm,2)))/
		    (ef + q0 - TMath::Sqrt(m2 + TMath::Power(kf + qm,2))))) + 
       TMath::Log(TMath::Abs(
		    (ef + q0 + TMath::Sqrt(m2 + TMath::Power(kf - qm,2)))/
		    (ef + q0 + TMath::Sqrt(m2 + TMath::Power(kf + qm,2))))));

  std::complex<double> uL3(
       TMath::Log(TMath::Abs((TMath::Power(2*kf + q0*ds,2)-qm2)/
			     (TMath::Power(2*kf - q0*ds,2)-qm2))) + 
       TMath::Log(TMath::Abs((TMath::Power(kf-ef*ds,2) - (4*m4*qm2)/q4)/
			     (TMath::Power(kf+ef*ds,2) - (4*m4*qm2)/q4))));

  std::complex<double> RlinrelX(-L1/(16.0*kPi2)+
				uL2*(2.0*ef+q0)/(32.0*kPi2*qm)-
				uL3*ds/(64.0*kPi2));

  return RlinrelX*16.0*m2;
}
//____________________________________________________________________________
//Following obtained from fortran code by J Nieves, which contained the following comment:
/*
   complex Lindhard function for symmetric nuclear matter:
                    from Appendix of
                    E.Oset et al Phys. Rept. 188:79, 1990
                    formula A.4 

   input variables: 
     q_zero [fm^-1] : Energy
     q_mod  [fm^-1] : Momentum
     rho    [fm^3]  : Nuclear density
     k_fermi[fm^-1] : Fermi momentum

   All variables are real*8

   output variable: 
     delta_lind [fm^-2]

            ATTENTION!!!
 Only works properly for real q_zero,
 if q_zero has an imaginary part calculates the L. function
 assuming Gamma= 0.
 Therefore this subroutine provides two different functions
 depending on whether q_zero is real or not!!!!!!!!!!!
*/
std::complex<double> NievesQELCCPXSec::deltaLindhard(double q0, 
                                double dq, double rho, double kF) const
{
  double q_zero = q0/fhbarc;
  double q_mod = dq/fhbarc;
  double k_fermi = kF/fhbarc;
  //Divide by hbarc in order to use natural units (rho is already in the correct units)
  
  //m = 939/197.3, md = 1232/197.3, mpi = 139/197.3
  double m = 4.7592;
  double md = 6.2433;
  double mpi = 0.7045;
  
  double fdel_f = 2.13;
  double wr = md-m;
  double gamma = 0;
  double gammap = 0;

  double q_zero2 =  TMath::Power(q_zero,  2);
  double q_mod2 =   TMath::Power(q_mod,   2);
  double k_fermi2 = TMath::Power(k_fermi, 2);

  double m2 =       TMath::Power(m,       2);
  double m4 =       TMath::Power(m,       4);
  double mpi2 =     TMath::Power(mpi,     2);
  double mpi4 =     TMath::Power(mpi,     4);

  double fdel_f2 =  TMath::Power(fdel_f,  2);
  
  //For the current code q_zero is always real
  //If q_zero can have an imaginary part then only the real part is used
  //until z and zp are calculated

  double s = m2+q_zero2-q_mod2+
    2.0*q_zero *TMath::Sqrt(m2+3.0/5.0*k_fermi2);

  if(s>TMath::Power(m+mpi,2)){
    double srot = TMath::Sqrt(s);
    double qcm = TMath::Sqrt(TMath::Power(s,2)+mpi4+m4-2.0*(s*mpi2+s*m2+
	      	mpi2*m2)) /(2.0*srot);
    gamma = 1.0/3.0 * 1.0/(4.0*kPi) * fdel_f2*
     TMath::Power(qcm,3)/srot*(m+TMath::Sqrt(m2+TMath::Power(qcm,2)))/mpi2;
  }
  double sp = m2+q_zero2-q_mod2-
    2.0*q_zero *TMath::Sqrt(m2+3.0/5.0*k_fermi2);


  if(sp > TMath::Power(m+mpi,2)){
    double srotp = TMath::Sqrt(sp);
    double qcmp = TMath::Sqrt(TMath::Power(sp,2)+mpi4+m4-2.0*(sp*mpi2+sp*m2+
		 mpi2*m2))/(2.0*srotp);
    gammap = 1.0/3.0 * 1.0/(4.0*kPi) * fdel_f2*
      TMath::Power(qcmp,3)/srotp*(m+TMath::Sqrt(m2+TMath::Power(qcmp,2)))/mpi2;
  }
  //}//End if statement
  const std::complex<double> iNum(0,1.0);

  std::complex<double> z(md/(q_mod*k_fermi)*(q_zero-q_mod2/(2.0*md)
                         -wr +iNum*gamma/2.0));
  std::complex<double> zp(md/(q_mod*k_fermi)*(-q_zero-q_mod2/(2.0*md)
                          -wr +iNum*gammap/2.0));
  
  std::complex<double> pzeta(0.0);
  if(abs(z) > 50.0){
    pzeta = 2.0/(3.0*z)+2.0/(15.0*z*z*z);
  }else if(abs(z) < TMath::Power(10.0,-2)){
    pzeta = 2.0*z-2.0/3.0*z*z*z-iNum*kPi/2.0*(1.0-z*z);
  }else{
    pzeta = z + (1.0-z*z) * log((z+1.0)/(z-1.0))/2.0;
  }

  std::complex<double> pzetap(0);
  if(abs(zp) > 50.0){
    pzetap = 2.0/(3.0*zp)+2.0/(15.0*zp*zp*zp);
  }else if(abs(zp) < TMath::Power(10.0,-2)){
    pzetap = 2.0*zp-2.0/3.0*zp*zp*zp-iNum*kPi/2.0*(1.0-zp*zp);
  }else{
    pzetap = zp+ (1.0-zp*zp) * log((zp+1.0)/(zp-1.0))/2.0;
  }

  //Multiply by hbarc^2 to give answer in units of GeV^2
  return 2.0/3.0 * rho * md/(q_mod*k_fermi) * (pzeta +pzetap) * fdel_f2 * 
    TMath::Power(fhbarc,2);
}

//____________________________________________________________________________
// Gives coulomb potential in units of GeV
double NievesQELCCPXSec::vcr(const Target * target, double Rcurr) const{
  if(target->IsNucleus()){
    int A = target->A();
    int Z = target->Z();
    // RMax calculated using formula from Nieves' fortran code and default
    // charge and neutron matter density paramters from NuclearUtils.cxx
    double Rmax;
    if(A > 20){
      double c = TMath::Power(A,0.35), z = 0.54;
      Rmax = c + 9.25*z;
    }else{
      // c = 1.75 for A <= 20
      Rmax = TMath::Sqrt(20.0)*1.75;
    }
    //LOG("Nieves",pDEBUG) "A = " << A
    //  << ", Rcurr = " << Rcurr << ", Rmax = " << Rmax;

    if(Rcurr >= Rmax){
      LOG("Nieves",pNOTICE) << "Radius greater than maximum radius for coulomb corrections."
			  << " Integrating to max radius.";
      Rcurr = Rmax;
    }
    
    int nInts = 1; // Number of intervals
    int nElts = 20*nInts; // Number of elements in x arrays

    std::vector<double> rvec = integrationSetup(0.0, Rcurr, nInts);
    std::vector<double> Vvec;
    Vvec.resize(nElts);

    for(int i=0; i<nElts; i++){
      double r = rvec[i];
      
      //Density gives the nuclear density, normalized to number p or n
      //Input radius r must be in fm
      double rhop = nuclear::Density(r,A,true);
      //double rhon = nuclear::Density(r,A,false);
      //double rho = rhop + rhon;

      Vvec[i] = rhop*r*r / Rcurr;
    }
    double result1 = integrate(0.0, Rcurr, nInts, Vvec);

    rvec = integrationSetup(Rcurr, Rmax, nInts);
    for(int i=0; i<nElts; i++){
      double r = rvec[i];
      double rhop = nuclear::Density(r,A,true);
      Vvec[i] = rhop*r;
    }
    double result2 = integrate(Rcurr, Rmax, nInts, Vvec);

    // Multiply by Z to normalize densities to number of protons
    // Multiply by hbarc to put result in GeV instead of fm
    return -kAem*4*kPi*(result1 + result2)*Z*fhbarc;
  }else{
    // If target is not a nucleus the potential will be 0
    return 0.0;
  }
}
  /* 
     To integrate, first call integrationSetup with the lower and upper limits
     a and b, and the number of intervals n. The returned vector x will be of
     size 20*n. For each element of x, the value of the function to be 
     integrated can be stored in a new vector y. ( y[i] = f(x[i]) )
     Finally, integrate is called with a,b,n, and y as the parameters,and 
     integrate will return the definite integral from a to b of the function.
  */
//____________________________________________________________________________
std::vector<double> NievesQELCCPXSec::integrationSetup(double a, double b, 
						       int n) const{
  int np = 20*n;
  double interval = (b-a)/double(n);
  double delt = interval*0.5;
  double orig = a-delt;

  int i1 = -21;
  int i2;
  int j1;
  int j2;

  std::vector<double> x;
  x.resize(np);
  for(int i=1;i<=n;i++){
    orig+=interval;
    i1 += 20;
    i2 = i1 + 21;
    for(int j=1;j<=10;j++){
      j1 = i1+j;
      j2 = i2-j;
      x[j1] = orig - delt*fIntervalFractions[j-1];
      x[j2] = orig + delt*fIntervalFractions[j-1];
      // fIntervalFractions is used to choose at which values in the interval
      // the function is evaluated.
    }
  }
  return x;
}
//____________________________________________________________________________
double NievesQELCCPXSec::integrate(double a, double b, int n,
				   std::vector<double> y) const{
  if(a == b)
    return 0.0;

  double weightedAvgSum = 0;
  // For each i value (ie each interval), the weighted avg of the function in 
  // the corresponding interval is added to weightedAvgSum. For example, if the
  // function is a constant c, then after each i, weightedAvgSum will be equal 
  // to c*i. We will need to multiply this by the interval length to get the 
  // integral.

  int i1 = -21;
  int i2;
  int j1;
  int j2;
  for(int i=1;i<=n;i++){
    i1 +=20;
    i2 = i1+21;
    for(int j=1;j<=10;j++){
      j1 = i1+j;
      j2 = i2-j;
      weightedAvgSum += fW[j-1]*0.5*(y[j1]+y[j2]);
      // fW is used to get a weighted sum in each interval
      // fW is normalized so the sum of its elements is 1
    }
  }
  return weightedAvgSum*(b-a)/double(n);
}
/*
//____________________________________________________________________________
std::complex<double> NievesQELCCPXSec::integrate(double a, double b, int n,
			         std::vector<std::complex<double> > y) const{
  std::complex<double> weightedAvgSum(0,0);
  // For each i value (ie each interval), the weighted avg of the function in 
  // the corresponding interval is added to weightedAvgSum. For example, if the
  // function is a constant c, then after each i, weightedAvgSum will be equal 
  // to c*i. We will need to multiply this by the interval length to get the 
  // integral.

  int i1 = -21;
  int i2;
  int j1;
  int j2;
  for(int i=1;i<=n;i++){
    i1 +=20;
    i2 = i1+21;
    for(int j=1;j<=10;j++){
      j1 = i1+j;
      j2 = i2-j;
      weightedAvgSum += fW[j-1]*0.5*(y[j1]+y[j2]);
      // fW is used to get a weighted sum in each interval 
      // fW is normalized so the sum of its elements is 1
    }
  }
  return weightedAvgSum*(b-a)/double(n);
}
*/
//____________________________________________________________________________
int NievesQELCCPXSec::leviCivita(int input[]) const{
  int copy[4] = {input[0],input[1],input[2],input[3]};
  int permutations = 0;
  int temp;

  for(int i=0;i<4;i++){
    for(int j=i+1;j<4;j++){
      //If any two elements are equal return 0
      if(input[i] == input[j])
	return 0;
      //If a larger number is before a smaller one, use permutations
      //(exchanges of two adjacent elements) to move the smaller element
      //so it is before the larger element, eg 2341->2314->2134->1234
      if(copy[i]>copy[j]){
	temp = copy[j];
	for(int k=j;k>i;k--){
	  copy[k] = copy[k-1];
	  permutations++;
	}
	copy[i] = temp;
      }
    }
  }

  if(permutations % 2 == 0){
    return 1;
  }else{
    return -1;
  }  
}
//____________________________________________________________________________
double NievesQELCCPXSec::LmunuAnumu(const Interaction * interaction,
				    const TLorentzVector piVec,
				    const TLorentzVector kVec,
				    const TLorentzVector kPrimeVec) const
{
  const TLorentzVector qVec = kVec-kPrimeVec;

  const Kinematics & kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();
  const Target & target = init_state.Tgt();

  double M = target.HitNucMass();
  double q2 = kinematics.q2();

  const double k[4] = {kVec.E(),kVec.Px(),kVec.Py(),kVec.Pz()};
  const double kPrime[4] = {kPrimeVec.E(),kPrimeVec.Px(),
			    kPrimeVec.Py(),kPrimeVec.Pz()};
	  
  const double q[4] = {qVec.E(),qVec.Px(),qVec.Py(),qVec.Pz()};
  double q0 = q[0];

  double dq = TMath::Sqrt(TMath::Power(q[1],2)+
			  TMath::Power(q[2],2)+TMath::Power(q[3],2));

  bool is_neutrino = pdg::IsNeutrino(init_state.ProbePdg());
  int sign = (is_neutrino) ? 1 : -1;

  // Calculate the QEL form factors
  fFormFactors.Calculate(interaction);    

  double F1V   = 0.5*fFormFactors.F1V();
  double xiF2V = 0.5*fFormFactors.xiF2V();
  double FA    = -fFormFactors.FA();

  // According to Nieves' paper, Fp = 2.0*M*FA/(kPionMass2-q2), but Llewelyn-
  // Smith uses Fp = 2.0*M^2*FA/(kPionMass2-q2), so I divide by M
  // This gives units of GeV^-1
  double Fp    = -1.0/M*fFormFactors.Fp(); 
  
  //double Fp    = 2.0*M*FA/(kPionMass2-q2); // Has units of GeV^-1

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Nieves", pDEBUG) << "\n" << fFormFactors;
#endif

  // Calculate auxiliary parameters
  double M2      = TMath::Power(M,     2);
  double FA2     = TMath::Power(FA,    2);
  double Fp2     = TMath::Power(Fp,    2);
  double F1V2    = TMath::Power(F1V,   2);
  double xiF2V2  = TMath::Power(xiF2V, 2);
  double q2_M2   = q2/M2;
  double q02     = TMath::Power(q[0],  2);
  double dq2     = TMath::Power(dq,    2);
  double q4      = TMath::Power(q2,    2);

  //Terms for xsec terms that don't have RPA corrections

  double a1 = 8.0*q2*((F1V2+2*F1V*xiF2V+xiF2V2)+FA2*(1.0/4.0-M2/q2));
  double a2 = 32.0*F1V2 - 8.0*xiF2V2*q2/M2 + 8.0*FA2;
  double a3 = sign*16.0*FA*(F1V+xiF2V);
  double a4 = -8.0*q2_M2*xiF2V2*(M2/q2+1.0/4.0) - 
    8.0*M2*FA2/(kPionMass2-q2)*(q2/(kPionMass2-q2)+2.0) - 16.0*F1V*xiF2V;

  double r;
  if(kinematics.KVSet(kKVSelRad))
    r = kinematics.GetKV(kKVSelRad);
  else
    r = 0.0;
  
  double t0,r00;
  double CN,CT,CL,imU;
  CNCTCLimUcalc(& target,r,qVec,q2,is_neutrino,CN,CT,CL,imU,t0,r00);

  if(! fRPA){
    CN=1.0;
    CT=1.0;
    CL=1.0;
  }

  //CN = 1.0, CT = 1.0, CL = 1.0;
  //LOG("Nieves",pDEBUG) << "CN=" << CN << ",CT=" << CT << ",CL=" << CL << ",imU=" << imU;

  // Use average values for initial momentum to calculate A, as given
  // in Appendix B of Nieves' paper. T gives average values of components
  // of p, and R gives the average value of two components multiplied
  // together
  double t3 = (0.5*q2 + q0*t0)/dq; // Average pz

  // Vector of p
  double tulin[4] = {t0,0.0,0.0,t3};

  // R is a 4x4 matrix, with R[mu][nu] is the average 
  // value of p[mu]*p[nu]
  double aR = r00-M2;
  double bR = (q4+4.0*r00*q02+4.0*q2*q0*t0)/(4.0*dq2);
  double gamma = (aR-bR)/2.0;
  double delta = (-aR+3.0*bR)/2.0/dq2;

  double r03 = (0.5*q2*t0 + q0*r00)/dq; // Average E(p)*pz

  double rulin[4][4] = { {r00,   0.0, 0.0,   r03}, 
			 {0.0, gamma, 0.0,   0.0},
			 {0.0,   0.0, gamma, 0.0},
			 {r03,   0.0, 0.0,   gamma+delta*dq2} };

  /*
  // Old values for testing purposes
  //double tulin[4];
  tulin[0] = piVec.E();
  tulin[1] = piVec.Px();
  tulin[2] = piVec.Py();
  tulin[3] = piVec.Pz();

  for(int i=0; i<4; i++)
    for(int j=0; j<4; j++)
      rulin[i][j] = tulin[i]*tulin[j];
  */

  /*
  for(int i=0; i<4; i++)
    LOG("Nieves",pDEBUG) <<  tulin[i];

  for(int i=0; i<4; i++)
    for(int j=0; j<4; j++)
      LOG("Nieves",pDEBUG) << rulin[i][j];
    */

  //Additional constants and variables
  const int g[4][4] = {{1,0,0,0},{0,-1,0,0},{0,0,-1,0},{0,0,0,-1}};
  const std::complex<double> iNum(0,1);
  int leviCivitaIndexArray[4];
  int lcNum = 0;
  double imaginaryPart = 0;

  std::complex<double> sum(0.0,0.0);

  double kPrimek = k[0]*kPrime[0]-k[1]*kPrime[1]-k[2]*kPrime[2]-k[3]*kPrime[3];

  std::complex<double> Lmunu(0.0,0.0),Lnumu(0.0,0.0);
  std::complex<double> Amunu(0.0,0.0),Anumu(0.0,0.0);
  
  // Write values to a file for testing purposes
  ofstream ffstream;
  if(fPrintData){
    ffstream.open("formfactors.txt", std::ios_base::app);
    ffstream << -q2 << "\t" << q[0] << "\t" << dq << "\t"
	     << F1V << "\t" << xiF2V << "\t" << FA << "\t" 
	     << Fp;
    if(-qVec.Mag2() > 0)
      ffstream << "\t" << -qVec.Mag2() << "\n";
    else
      ffstream << "\n";
    ffstream.close();
    
    ffstream.open("cnctcl.txt", std::ios_base::app);
    ffstream << -q2 << "\t" << CN << "\t" << CT << "\t" << CL << "\n";
    ffstream.close();
    
    ffstream.open("tr.txt", std::ios_base::app);
    ffstream << -q2 << "\t" << tulin[0] << "\t" << tulin[3] << "\t" << rulin[0][0] << "\t"
	     << rulin[0][3] << "\t" << rulin[1][1] << "\t" << rulin[3][3] << "\n";
    ffstream.close();
    
    ffstream.open("Amunu.txt", std::ios_base::app);
    ffstream << -q2 << "\t";
  }

  // Calculate LmunuAnumu by iterating over mu and nu
  // In each iteration, add LmunuAnumu to sum if mu=nu, and add
  // LmunuAnumu + LnumuAmunu if mu != nu, since we start nu at mu
  for(int mu=0;mu<4;mu++){
    for(int nu=mu;nu<4;nu++){
      imaginaryPart = 0;
      if(mu == nu){
	//if mu==nu then levi-civita = 0, so imaginary part = 0
	Lmunu = g[mu][mu]*kPrime[mu]*g[nu][nu]*k[nu]+g[nu][nu]*kPrime[nu]*g[mu][mu]*k[mu]-g[mu][nu]*kPrimek;
      }else{
	//if mu!=nu, then g[mu][nu] = 0
	//This same leviCivitaIndex array can be used in the else portion when
	//calculating Anumu
	leviCivitaIndexArray[0] = mu;
	leviCivitaIndexArray[1] = nu;
	for(int a=0;a<4;a++){
	  for(int b=0;b<4;b++){
	    leviCivitaIndexArray[2] = a;
	    leviCivitaIndexArray[3] = b;
	    imaginaryPart += - leviCivita(leviCivitaIndexArray)*kPrime[a]*k[b];
	  }
	}
	//real(Lmunu) is symmetric, and imag(Lmunu) is antisymmetric
	//std::complex<double> num(g[mu][mu]*kPrime[mu]*g[nu][nu]*k[nu]+g[nu][nu]*kPrime[nu]*g[mu][mu]*k[mu],imaginaryPart);
	Lmunu = g[mu][mu]*kPrime[mu]*g[nu][nu]*k[nu]+g[nu][nu]*kPrime[nu]*g[mu][mu]*k[mu] + iNum*imaginaryPart;
	Lnumu = g[nu][nu]*kPrime[nu]*g[mu][mu]*k[mu]+g[mu][mu]*kPrime[mu]*g[nu][nu]*k[nu ]- iNum*imaginaryPart;
      } // End Lmunu calculation
      
      if(mu ==0 && nu == 0){
	Amunu = 16.0*F1V2*(2.0*rulin[0][0]*CN+2.0*q[0]*tulin[0]+q2/2.0)+ 
	  2.0*q2*xiF2V2*
	  (4.0-4.0*rulin[0][0]/M2-4.0*q[0]*tulin[0]/M2-q02*(4.0/q2+1.0/M2)) +
	  4.0*FA2*(2.0*rulin[0][0]+2.0*q[0]*tulin[0]+(q2/2.0-2.0*M2))-
	  (2.0*CL*Fp2*q2+8.0*FA*Fp*CL*M)*q02-16.0*F1V*xiF2V*(-q2+q02)*CN;
	if(fPrintData) ffstream << real(Amunu) << "\t";
	sum += Lmunu*Amunu;
      }else if(mu == 0 && nu == 3){
	Amunu = 16.0*F1V2*((2.0*rulin[0][3]+tulin[0]*dq)*CN+tulin[3]*q[0])+
	  2.0*q2*xiF2V2*
	  (-4.0*rulin[0][3]/M2-2.0*(dq*tulin[0]+q[0]*tulin[3])/M2-dq*q[0]*(4.0/q2+1.0/M2))+
	  4.0*FA2*((2.0*rulin[0][3]+dq*tulin[0])*CL+q[0]*tulin[3])-
	  (2.0*CL*Fp2*q2+8.0*FA*Fp*CL*M)*dq*q[0]-
	  16.0*F1V*xiF2V*dq*q[0];
	if(fPrintData) ffstream << real(Amunu) << "\t";
	Anumu = Amunu;
	sum += Lmunu*Anumu + Lnumu*Amunu;
      }else if(mu == 3 && nu == 3){
	Amunu = 16.0*F1V2*(2.0*rulin[3][3]+2.0*dq*tulin[3]-q2/2.0)+ 
	  2.0*q2*xiF2V2*(-4.0-4.0*rulin[3][3]/M2-4.0*dq*tulin[3]/M2-dq2*(4.0/q2+1.0/M2))+
	  4.0*FA2*(2.0*rulin[3][3]+2.0*dq*tulin[3]-(q2/2.0-2.0*CL*M2))-
	  (2.0*CL*Fp2*q2+8.0*FA*Fp*CL*M)*dq2-
	  16.0*F1V*xiF2V*(q2+dq2);
	if(fPrintData) ffstream << real(Amunu) << "\t";
	sum += Lmunu*Amunu;
      }else if(mu ==1 && nu == 1){
	Amunu = 16.0*F1V2*(2.0*rulin[1][1]-q2/2.0)+ 
	  2.0*q2*xiF2V2*(-4.0*CT-4.0*rulin[1][1]/M2) +
	  4.0*FA2*(2.0*rulin[1][1]-(q2/2.0-2.0*CT*M2))-
	  16.0*F1V*xiF2V*CT*q2;
	if(fPrintData) ffstream << real(Amunu) << "\t";
	sum += Lmunu*Amunu;
      }else if(mu == 2 && nu == 2){
	// Ayy not explicitly listed in paper. This is included so rotating the
	// coordinates of k and k' about the z-axis does not change the xsec.
	Amunu = 16.0*F1V2*(2.0*rulin[2][2]-q2/2.0)+ 
	  2.0*q2*xiF2V2*(-4.0*CT-4.0*rulin[2][2]/M2) +
	  4.0*FA2*(2.0*rulin[2][2]-(q2/2.0-2.0*CT*M2))-
	  16.0*F1V*xiF2V*CT*q2;
	sum += Lmunu*Amunu;
      }else if(mu ==1 && nu == 2){
	Amunu = sign*16.0*iNum*FA*(xiF2V+F1V)*(-dq*tulin[0]*CT + q[0]*tulin[3]);
	Anumu = -Amunu; // Im(A) is antisymmetric
	if(fPrintData) ffstream << imag(Amunu) << "\t";
	sum += Lmunu*Anumu+Lnumu*Amunu;
      }else{
	// No RPA corrections to the remaining terms, so A is not r dependent
	if(mu == nu){
	  Amunu = a1*g[mu][nu]+
	    a2*(rulin[mu][nu]+tulin[mu]*q[nu]/2.0+tulin[nu]*q[mu]/2.0)+a4*q[mu]*q[nu];
	  sum += Lmunu*Amunu;
	}else{
	  imaginaryPart = 0;
	  leviCivitaIndexArray[0] = mu;
	  leviCivitaIndexArray[1] = nu;
	  for(int a=0;a<4;a++){
	    for(int b=a+1;b<4;b++){
	      leviCivitaIndexArray[2] = a;
	      leviCivitaIndexArray[3] = b;
	      //Switching a and b reverses the sign of the leviCivita
	      //symbol, and a==b gives leviCivita = 0
	      lcNum = leviCivita(leviCivitaIndexArray);
	      imaginaryPart += lcNum*g[a][a]*tulin[a]*g[b][b]*q[b] 
		-lcNum*g[b][b]*tulin[b]*g[a][a]*q[a];
	    }
	  }
	  Amunu = a1*g[mu][nu]+
	    a2*(rulin[mu][nu]+tulin[mu]*q[nu]/2.0+tulin[nu]*q[mu]/2.0)+
	    iNum*a3*imaginaryPart+a4*q[mu]*q[nu];
	  Anumu = a1*g[mu][nu]+
	    a2*(rulin[mu][nu]+tulin[mu]*q[nu]/2.0+tulin[nu]*q[mu]/2.0)-
	    iNum*a3*imaginaryPart+a4*q[mu]*q[nu];
	  sum += Lmunu*Anumu+Lnumu*Amunu;
	}
      }
    } // End for loop over nu
  } // End for loop over mu

  if(fPrintData){
    ffstream << "\n";
    ffstream.close();
  }

  // Since the real parts of A and L are both symmetric and the imaginary
  // parts are antisymmetric, the contraction should be real
  if(imag(sum) > TMath::Power(10.0,-4))
    LOG("Nieves",pWARN) << "Imaginary part of tensor contraction is nonzero "
			<< "in QEL XSec, real(sum) = " << real(sum) 
			<< "imag(sum) = " << imag(sum);
  return real(sum);
}
//____________________________________________________________________________
