//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Joe Johnston, University of Pittsburgh
         Steven Dytman, University of Pittsburgh

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <iostream>
#include <fstream>
#include <complex>

#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <ctype.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>

#include "Algorithm/AlgFactory.h"
#include "Algorithm/AlgConfigPool.h"
#include "Base/XSecIntegratorI.h"
#include "Base/QELFormFactors.h"
#include "Base/QELFormFactorsModelI.h"
#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Conventions/KineVar.h"
#include "Conventions/Units.h"
#include "EVGModules/VertexGenerator.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepStatus.h"
#include "Messenger/Messenger.h"
#include "Nuclear/NuclearModelI.h"
#include "Nuclear/FermiMomentumTablePool.h"
#include "Nuclear/FermiMomentumTable.h"
#include "LlewellynSmith/NievesQELCCPXSec.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGUtils.h"
#include "QEL/QELPrimaryLeptonGenerator.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/NuclearUtils.h"
#include "Utils/PrintUtils.h"
#include "Utils/GSLUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
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
  // TESTING CODE:
  // The first time this method is called, output tensor elements and other
  // kinmeatics variables for various kinematics. This can the be compared
  // to Nieves' fortran code for validation purposes
  //bool exit_prog = false;
  /*if(fCompareNievesTensors){
    LOG("Nieves",pNOTICE) << "Printing tensor elements for specific " 
			  << "kinematics for testing purposes";
    CompareNievesTensors(interaction);
    fCompareNievesTensors = false;
    exit(0);
    }*/
  // END TESTING CODE


  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;


  // Get kinematics and init-state parameters
  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();
  const Target & target = init_state.Tgt();

  // Get the four kinematic vectors and caluclate GFactor
  // Create copies of all kinematics, so they can be rotated
  // and boosted to the nucleon rest frame (because the tensor
  // constraction below only applies for the initial nucleon
  // at rest and q in the z direction)
  TLorentzVector neutrinoMom,inNucleonMom,leptonMom, outNucleonMom;
  double Gfactor = 0.;
  if(kps == kPSTnctnBnctl){
    // All kinematics will already be stored
    TLorentzVector * tempNeutrino = init_state.GetProbeP4(kRfLab);
    neutrinoMom = *tempNeutrino;
    delete tempNeutrino;
    inNucleonMom = target.HitNucP4();
    leptonMom = kinematics.FSLeptonP4();
    outNucleonMom = kinematics.HadSystP4();

    // Use these lab kinematics to calculate the Gfactor, in order to make
    // the XSec differential in initial nucleon momentum and energy
    // Divide by 4.0 because Nieves' conventions for the leptonic and hadronic
    // tensor contraction differ from LwlynSmith by a factor of 4
    Gfactor = kGF2*fCos8c2 / (8.0*kPi*kPi*inNucleonMom.E()*neutrinoMom.E()*outNucleonMom.E()*leptonMom.E()) / 4.0;
  }else{
    // Initial Neutrino, Lab frame
    TLorentzVector * tempNeutrino = init_state.GetProbeP4(kRfLab);
    neutrinoMom = *tempNeutrino;
    delete tempNeutrino;
    // Initial Nucleon
    inNucleonMom = target.HitNucP4();

    // Generate outgoing lepton in the lab frame
    leptonMom = GenerateOutgoingLepton(interaction,neutrinoMom);

    // Set outgoing nucleon using conservation of energy
    outNucleonMom = neutrinoMom + inNucleonMom - leptonMom;

    // Calculate G factor
    double E = init_state.ProbeE(kRfLab);
    double M  = target.HitNucMass();
    double M2 = TMath::Power(M,     2);
    double s = (2*E+M)*M;
    double num = TMath::Power(s-M2,2);
    Gfactor = kGF2 * fCos8c2 / (8*kPi*num);
  }

  // Boost to nucleon rest frame to calculate the nucleon rest frame cross section
  TVector3 beta = -1.0 * inNucleonMom.BoostVector(); // boost from lab to nucRest
  neutrinoMom.Boost(beta);
  inNucleonMom.Boost(beta);
  leptonMom.Boost(beta);
  outNucleonMom.Boost(beta);
  
  // Rotate vectors so q is in the z direction, to use Nieves'
  // explicit form of the Amunu tensor
  TVector3 neutrinoMom3 = neutrinoMom.Vect();
  TVector3 leptonMom3 = leptonMom.Vect();
  TVector3 q3Vec = neutrinoMom3-leptonMom3; // TESTING: Use q rather than qTilde

  TVector3 inNucleonMom3 = inNucleonMom.Vect();
  TVector3 outNucleonMom3 = outNucleonMom.Vect();
  //TVector3 q3Vec = outNucleonMom3-inNucleonMom3; // qTilde

  TVector3 zvec(0,0,1.0);
  TVector3 rot = (q3Vec.Cross(zvec)).Unit(); // Vector to rotate about
  double angle = zvec.Angle(q3Vec); // Angle between the z direction and q

  // Rotate if the rotation vector is not 0
  if(rot.Mag() >= kASmallNum){
    neutrinoMom3.Rotate(angle,rot);
    neutrinoMom.SetVect(neutrinoMom3);
    leptonMom3.Rotate(angle,rot);
    leptonMom.SetVect(leptonMom3);
    inNucleonMom3.Rotate(angle,rot);
    inNucleonMom.SetVect(inNucleonMom3);
    outNucleonMom3.Rotate(angle,rot);
    outNucleonMom.SetVect(outNucleonMom3);
  }

  // Calculate q and qTilde
  //TLorentzVector qP4(0,0,0,0);
  TLorentzVector qTildeP4(0,0,0,0);
  //qP4 = neutrinoMom - leptonMom;
  //qTildeP4 = outNucleonMom - inNucleonMom;

  //TESTING: use q instead of qTilde
  qTildeP4 = neutrinoMom - leptonMom;

  double Q2tilde = -1 * qTildeP4.Mag2();
  if(kps==kPSTnctnBnctl) // otherwise q2 is already stored
    interaction->KinePtr()->SetQ2(Q2tilde);
  
  /*LOG("Nieves",pINFO) << "qTildep4 = " << utils::print::P4AsString(&qTildeP4);
    LOG("Nieves",pINFO	<< "neutrinop4 = " << utils::print::P4AsString(&neutrinoMom);
    LOG("Nieves",pINFO  << "inNucleon = " << utils::print::P4AsString(&inNucleonMom);
    LOG("Nieves",pINFO  << "lepton = " << utils::print::P4AsString(&leptonMom);
    LOG("Nieves",pINFO  << "outNucleon = " << utils::print::P4AsString(&outNucleonMom);*/

  /*std::cout << "Q2tilde = " <<Q2tilde << ", " <<  utils::print::P4AsString(&qTildeP4)<< std::endl;
  std::cout << "Q2 (not tilde) = " << -1 * qP4.Mag2() << ", " <<  utils::print::P4AsString(&qP4)<< std::endl;
  std::cout << "Q2 difference (tilde - not) = " << Q2tilde + qP4.Mag2() << std::endl;
  std::cout << "Q2 stored = " << interaction->KinePtr()->Q2() << std::endl;*/

  double q2 = -Q2tilde;

  // Check that q2 < 0 (accounting for rounding errors)
  if(q2>=kASmallNum){
    LOG("Nieves", pWARN) << "q2>=0, returning xsec = 0.0";
    return 0.0;
  }
  // Check that the energy tranfer q0 is greater than 0, or else the
  // following equations do not apply. (Note also that the event would
  // be Pauli blocked )
  if(qTildeP4.E()<=-kASmallNum){
    LOG("Nieves", pWARN) << "q0<=0.0, returning xsec = 0.0";
    return 0.0;
  }

  // Calculate tensor product

  fFormFactors.Calculate(interaction);
  bool is_neutrino = pdg::IsNeutrino(init_state.ProbePdg());
  double M = inNucleonMom.Mag(); // Struck nucleon mass
  double r = target.HitNucPosition();
  bool tgtIsNucleus = target.IsNucleus();
  int tgt_pdgc = target.HitNucPdg();
  int A = target.A();
  int Z = target.Z();
  int N = target.N();
  bool hitNucIsProton = pdg::IsProton(target.HitNucPdg());
  double LmunuAnumuResult = LmunuAnumu(neutrinoMom,inNucleonMom,
				       leptonMom,outNucleonMom,
				       M,r,is_neutrino,tgtIsNucleus,
				       tgt_pdgc,A,Z,N,hitNucIsProton);

  // Calculate Coulomb corrections
  double ml = interaction->FSPrimLepton()->Mass();
  double ml2 = TMath::Power(ml,    2);
  double coulombFactor = 1.0;
  if(fCoulomb){
    // Coulomb potential
    double Vc = vcr(& target, r);

    // Outgoing lepton energy and momentum including coulomb potential
    int sign = is_neutrino ? 1 : -1;
    double El = leptonMom.E();
    double ElLocal = El - sign*Vc;
    if(ElLocal - ml <= 0.0){
      LOG("Nieves",pWARN) << "Event should be rejected. Coulomb effects "
			  << "push kinematics below threshold. Returning "
			  << "xsec = 0.0";
      return 0.0;
    }
    double plLocal = TMath::Sqrt(ElLocal*ElLocal-ml2);

    // Correction factor
    coulombFactor= plLocal*ElLocal/leptonMom.Vect().Mag()/El;
    fElocal = ElLocal;
    fPlocal = plLocal;
    fElep = El;
    fPlep = leptonMom.Vect().Mag();
  }

  // Calculate xsec
  double xsec = Gfactor*coulombFactor*LmunuAnumuResult;

  LOG("Nieves",pDEBUG) << "TESTING: RPA=" << fRPA 
		       << ", Coulomb=" << fCoulomb
		       << ", q2 = " << q2 << ", xsec = " << xsec;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Nieves", pDEBUG)
     << "dXSec[QEL]/dQ2 [FreeN](E = "<< E << ", Q2 = "<< -q2 << ") = "<< xsec;
#endif

  //----- The algorithm computes dxsec/dQ2 or kPSTnctnBnctl
  //      Check whether variable tranformation is needed
  if(kps!=kPSQ2fE && kps!=kPSTnctnBnctl) {
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
  if( interaction->TestBit(kIAssumeFreeNucleon) )    return xsec;

  //----- compute nuclear suppression factor
  //      (R(Q2) is adapted from NeuGEN - see comments therein)
  double R;
  R = nuclear::NuclQELXSecSuppression("Default", 0.5, interaction);

  //----- number of scattering centers in the target
  int nucpdgc = target.HitNucPdg();
  int NNucl = (pdg::IsProton(nucpdgc)) ? target.Z() : target.N(); 

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Nieves", pINFO) 
       << "Nuclear suppression factor R(Q2) = " << R << ", NNucl = " << NNucl;
#endif

  xsec *= (R*NNucl); // nuclear xsec*/

  return xsec;
}
//____________________________________________________________________________
double NievesQELCCPXSec::Integral(const Interaction * in) const
{
  bool nuclear_target = in->InitState().Tgt().IsNucleus();
  double E = in->InitState().ProbeE(kRfHitNucRest);
  if(!nuclear_target || !fDoAvgOverNucleonMomentum) {
    return fXSecIntegrator->Integrate(this,in);
  }



  // If the nuclear model is LFG or if RPA effects are on, then the xsec
  // is dependent on the position in the nucleus, and possible radii and
  // initial nucleon momenta should always be averaged over
  if(fLFG || fRPA || E < fEnergyCutOff) {
    // clone the input interaction so as to tweak the
    // hit nucleon 4-momentum in the averaging loop
    Interaction in_curr(*in);

    // hit target
    Target * tgt = in_curr.InitState().TgtPtr();

    // get nuclear masses (init & final state nucleus)
    int nucleon_pdgc = tgt->HitNucPdg();
    bool is_p = pdg::IsProton(nucleon_pdgc);
    int Zi = tgt->Z();
    int Ai = tgt->A();
    int Zf = (is_p) ? Zi-1 : Zi;
    int Af = Ai-1;
    PDGLibrary * pdglib = PDGLibrary::Instance();
    TParticlePDG * nucl_i = pdglib->Find( pdg::IonPdgCode(Ai, Zi) );
    TParticlePDG * nucl_f = pdglib->Find( pdg::IonPdgCode(Af, Zf) );
    if(!nucl_f) {
      LOG("QELXSec", pFATAL)
	<< "Unknwown nuclear target! No target with code: "
	<< pdg::IonPdgCode(Af, Zf) << " in PDGLibrary!";
      exit(1);
    }
    double Mi  = nucl_i -> Mass(); // initial nucleus mass
    double Mf  = nucl_f -> Mass(); // remnant nucleus mass

    // throw nucleons with fermi momenta and binding energies 
    // generated according to the current nuclear model for the
    // input target and average the cross section
    double xsec_sum = 0.;
    const int nnuc = 2000;
    // VertexGenerator for generating a position before generating
    // each nucleon
    VertexGenerator * vg = new VertexGenerator();
    vg->Configure("Default");
    for(int inuc=0;inuc<nnuc;inuc++){
      // Generate a position in the nucleus
      TVector3 nucpos = vg->GenerateVertex(&in_curr,tgt->A());
      tgt->SetHitNucPosition(nucpos.Mag());

      // Generate a nucleon
      fNuclModel->GenerateNucleon(*tgt, nucpos.Mag());
      TVector3 p3N = fNuclModel->Momentum3();
      double   EN  = Mi - TMath::Sqrt(p3N.Mag2() + Mf*Mf);
      TLorentzVector* p4N = tgt->HitNucP4Ptr();
      p4N->SetPx (p3N.Px());
      p4N->SetPy (p3N.Py());
      p4N->SetPz (p3N.Pz());
      p4N->SetE  (EN);

      double xsec;
      xsec = fXSecIntegrator->Integrate(this,&in_curr);

      xsec_sum += xsec;
    }
    double xsec_avg = xsec_sum / nnuc;
    delete vg;
    return xsec_avg;
  }else{
    return fXSecIntegrator->Integrate(this,in);
  }
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

  fPrintData = fConfig->GetBoolDef("PrintDebugData",false);
  //fPrintData = false;
  fCompareNievesTensors = fPrintData;

  // Get nuclear model for use in Integral()
  RgKey nuclkey = "IntegralNuclearModel";
  fNuclModel = dynamic_cast<const NuclearModelI *> (this->SubAlg(nuclkey));
  assert(fNuclModel);

  // Check if the model is a local Fermi gas
  fLFG = fNuclModel->ModelType(Target()) == kNucmLocalFermiGas;
  
  if(!fLFG){
    // get the Fermi momentum table for relativistic Fermi gas
    fKFTableName = fConfig->GetStringDef ("FermiMomentumTable",
					  gc->GetString("FermiMomentumTable"));
    fKFTable = 0;
    
    FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
    fKFTable = kftp->GetTable(fKFTableName);
    assert(fKFTable);
  }

  // Always average over initial nucleons if the nuclear model is LFG
  fDoAvgOverNucleonMomentum =
    fLFG || fConfig->GetBoolDef("IntegralAverageOverNucleonMomentum", false);

  fEnergyCutOff = 0.;

  if(fDoAvgOverNucleonMomentum) {
    // Get averaging cutoff energy
    fEnergyCutOff = 
      fConfig->GetDoubleDef("IntegralNuclearInfluenceCutoffEnergy", 2.0);
  }

}
//___________________________________________________________________________
void NievesQELCCPXSec::CNCTCLimUcalc(TLorentzVector qTildeP4,
				     double M, double r, bool is_neutrino,
				     bool tgtIsNucleus, int tgt_pdgc,
				     int A, int Z, int N, bool hitNucIsProton,
				     double & CN, double & CT, double & CL,
				     double & imaginaryU,
				     double & t0, double & r00) const
{
  if(tgtIsNucleus){
    double dq = qTildeP4.Vect().Mag();
    double dq2 = TMath::Power(dq,2);
    double q2 = 1 * qTildeP4.Mag2();
    //Terms for polarization coefficients CN,CT, and CL
    double hbarc2 = TMath::Power(fhbarc,2);
    double c0 = 0.380/fhbarc;//Constant for CN in natural units

    //Density gives the nuclear density, normalized to 1
    //Input radius r must be in fm
    double rhop = nuclear::Density(r,A)*Z;
    double rhon = nuclear::Density(r,A)*N;
    double rho = rhop + rhon;
    double rho0 = A*nuclear::Density(0,A);

    // TESTING CODE
    rhopStored = rhop;
    rhonStored = rhon;
    rhoStored = rho;
    rho0Stored = rho0;
    
    double fPrime = (0.33*rho/rho0+0.45*(1-rho/rho0))*c0;
    
    // Get Fermi momenta
    double kF1, kF2;
    if(fLFG){
      if(hitNucIsProton){
	kF1 = TMath::Power(3*kPi2*rhop, 1.0/3.0) *fhbarc;
	kF2 = TMath::Power(3*kPi2*rhon, 1.0/3.0) *fhbarc;
      }else{
	kF1 = TMath::Power(3*kPi2*rhon, 1.0/3.0) *fhbarc;
	kF2 = TMath::Power(3*kPi2*rhop, 1.0/3.0) *fhbarc;
      }
    }else{
      if(hitNucIsProton){
	kF1 = fKFTable->FindClosestKF(tgt_pdgc, kPdgProton);
	kF2 = fKFTable->FindClosestKF(tgt_pdgc, kPdgNeutron);
      }else{
	kF1 = fKFTable->FindClosestKF(tgt_pdgc, kPdgNeutron);
	kF2 = fKFTable->FindClosestKF(tgt_pdgc, kPdgProton);	
      }
    }

    /*fc0 = c0;
    fPrimeStored = fPrime;
    fKF1 = kF1;
    fKF2 = kF2;*/

    //LOG("Nieves",pDEBUG) << "r=" << r << ",kF1=" << kF1 << ",kF2=" << kF2;
    //double kF = kF1 + kF2;
    double kF = TMath::Power(1.5*kPi2*rho, 1.0/3.0) *fhbarc;

    std::complex<double> imU(relLindhardIm(qTildeP4.E(),dq,kF1,kF2,
					   M,is_neutrino,t0,r00));

    imaginaryU = imag(imU);
    //LOG("Nieves",pDEBUG) << "imU = " << imaginaryU;
    std::complex<double> relLin(0,0),udel(0,0);

    /* LOG("Nieves",pDEBUG) << "q0 = " <<qTildeP4.E() << ", dq = " << dq
			 << ", kF = " << kF << ", M = " << M 
			 << ", is_neu = " << is_neutrino
			 << ", imU = " << imU;*/

    // By comparison with Nieves' fortran code
    if(imaginaryU < 0.){
      relLin = relLindhard(qTildeP4.E(),dq,kF,M,is_neutrino,imU);
      udel = deltaLindhard(qTildeP4.E(),dq,rho,kF);
    }
    std::complex<double> relLinTot(relLin + udel);

    // UPDATE THESE MASSES TO USE GENIE VALUES

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
    
    CN = 1.0/TMath::Power(abs(1.0-fPrime*relLin/hbarc2),2);

    CT = 1.0/TMath::Power(abs(1.0-relLinTot*Vt),2);
    CL = 1.0/TMath::Power(abs(1.0-relLinTot*Vl),2);
    LOG("Nieves",pDEBUG) <<"CN = " << CN <<",CT = " << CT << ",CL = " << CL;

    /*fVt = Vt;
    fVl = Vl;
    fRelLinReal = real(relLin);
    fRelLinIm = imag(relLin);
    fRelLinTotReal = real(relLinTot);
    fRelLinTotIm = imag(relLinTot);*/
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
  //if(fNievesSuppression && !interaction->TestBit(kIAssumeFreeNucleon )
  //&& !EF1-epsRP<0){
  //LOG("Nieves", pINFO) << "Average value of E(p) above Fermi sea";
  //return 0;
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
  double m = M/fhbarc;

  if(q0>qm){
    LOG("Nieves", pWARN) << "relLindhard() failed";
    return 0.0;
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
  
  /*q2rellin = q2;
  fl1 = L1;
  fl2 = real(uL2);
  fl3 = real(uL3);
  fkf = kf;
  fef = ef;
  fl2im = imag(uL2);
  fl3im = imag(uL3);*/

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
    //Rmax = 7.558; // TESTING PURPOSES
    //LOG("Nieves",pDEBUG) "A = " << A
    //  << ", Rcurr = " << Rcurr << ", Rmax = " << Rmax;

    if(Rcurr >= Rmax){
      LOG("Nieves",pNOTICE) << "Radius greater than maximum radius for coulomb corrections."
			  << " Integrating to max radius.";
      Rcurr = Rmax;
    }

    ROOT::Math::IBaseFunctionOneDim * func = new
      utils::gsl::wrap::NievesQELvcrIntegrand(Rcurr,A,Z);
    ROOT::Math::IntegrationOneDim::Type ig_type =
      utils::gsl::Integration1DimTypeFromString("adaptive");

    double abstol = 1; // We mostly care about relative tolerance;
    double reltol = 1E-4;
    int nmaxeval = 100000;
    ROOT::Math::Integrator ig(*func,ig_type,abstol,reltol,nmaxeval);
    double result = ig.Integral(0,Rmax);
    delete func;
    
    /*    int nInts = 1; // Number of intervals
    int nElts = 20*nInts; // Number of elements in x arrays

    std::vector<double> rvec = integrationSetup(0.0, Rcurr, nInts);
    std::vector<double> Vvec;
    Vvec.resize(nElts);

    //ofstream rhostream;
    //rhostream.open("r_rho.txt", std::ios_base::app);
    for(int i=0; i<nElts; i++){
      double r = rvec[i];
      
      //Density gives the nuclear density, normalized to number p or n
      //Input radius r must be in fm
      double rhop = Z*nuclear::Density(r,A);

      Vvec[i] = rhop*r*r / Rcurr;
      //if(fCompareNievesTensors) rhostream << r << "\t" << rhop << "\t" << Vvec[i] << "\t" << "\n";
    }
    double result1 = integrate(0.0, Rcurr, nInts, Vvec);

    rvec = integrationSetup(Rcurr, Rmax, nInts);
    for(int i=0; i<nElts; i++){
      double r = rvec[i];
      double rhop = Z*nuclear::Density(r,A);
      Vvec[i] = rhop*r;
      //if(fCompareNievesTensors) rhostream << r << "\t" << rhop << "\t" << Vvec[i] << "\t" << "\n";
    }
    double result2 = integrate(Rcurr, Rmax, nInts, Vvec);*/
    //rhostream.close();
    // Multiply by Z to normalize densities to number of protons
    // Multiply by hbarc to put result in GeV instead of fm
    frmax = Rmax;
    frcurr = Rcurr;
    //fresult1 = result1;
    //fresult2 = result2;
    fresult1 = result;
    falpha = kAem;

    //return -kAem*4*kPi*(result1 + result2)*fhbarc;
    return -kAem*4*kPi*result*fhbarc;
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
// neutrinoMom and leptonMom only affect the leptonic tensor
// inNucleonMom and outNucleonMom are only used to calculate q
// Nieves' averages for intial nucleon momenta are used in all other places
double NievesQELCCPXSec::LmunuAnumu(const TLorentzVector neutrinoMom,
				    const TLorentzVector inNucleonMom,
				    const TLorentzVector leptonMom,
				    const TLorentzVector outNucleonMom,
				    double M, double r, bool is_neutrino, 
				    bool tgtIsNucleus,
				    int tgt_pdgc, int A, int Z, int N,
				    bool hitNucIsProton) const
{
  const double k[4] = {neutrinoMom.E(),neutrinoMom.Px(),neutrinoMom.Py(),neutrinoMom.Pz()};
  const double kPrime[4] = {leptonMom.E(),leptonMom.Px(),
			    leptonMom.Py(),leptonMom.Pz()};
	  
  //const TLorentzVector qTildeP4 = outNucleonMom-inNucleonMom;
  const TLorentzVector qTildeP4 = neutrinoMom-leptonMom;
  double q2 = qTildeP4.Mag2();


  const double q[4] = {qTildeP4.E(),qTildeP4.Px(),qTildeP4.Py(),qTildeP4.Pz()};
  double q0 = q[0];
  double dq = TMath::Sqrt(TMath::Power(q[1],2)+
			  TMath::Power(q[2],2)+TMath::Power(q[3],2));

  int sign = (is_neutrino) ? 1 : -1;

  // Get the QEL form factors (were calculated before this method was called)
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

  double t0,r00;
  double CN=1.,CT=1.,CL=1.,imU=0; // NOTE: imU can be removed as an argument here
  CNCTCLimUcalc(qTildeP4,M,r,is_neutrino,tgtIsNucleus,
		tgt_pdgc,A,Z,N,hitNucIsProton,
		CN,CT,CL,imU,t0,r00);
  
  if(! fRPA){
    CN=1.0;
    CT=1.0;
    CL=1.0;
  }

  //CN = 1.0, CT = 1.0, CL = 1.0;
  //LOG("Nieves",pDEBUG) << "CN=" << CN << ",CT=" << CT << ",CL=" << CL << ",imU=" << imU;

  
  double tulin[4] = {0.,0.,0.,0.};
  double rulin[4][4] = { {0.,0.,0.,0.}, 
			 {0.,0.,0.,0.}, 
			 {0.,0.,0.,0.}, 
			 {0.,0.,0.,0.} };
  if(fCompareNievesTensors){
    // Use average values for initial momentum to calculate A, as given
    // in Appendix B of Nieves' paper. T gives average values of components
    // of p, and R gives the average value of two components multiplied
    // together
    double t3 = (0.5*q2 + q0*t0)/dq; // Average pz
    
    // Vector of p
    
    tulin[0] = t0;
    tulin[3] = t3;
    
    // R is a 4x4 matrix, with R[mu][nu] is the average 
    // value of p[mu]*p[nu]
    double aR = r00-M2;
    double bR = (q4+4.0*r00*q02+4.0*q2*q0*t0)/(4.0*dq2);
    double gamma = (aR-bR)/2.0;
    double delta = (-aR+3.0*bR)/2.0/dq2;
    
    double r03 = (0.5*q2*t0 + q0*r00)/dq; // Average E(p)*pz
    
    rulin[0][0] = r00;
    rulin[0][3] = r03;
    rulin[1][1] = gamma;
    rulin[2][2] = gamma;
    rulin[3][0] = r03;
    rulin[3][3] = gamma+delta*dq2;
  }else{
    // For normal code execulation, tulin is the initial nucleon momentum
    tulin[0] = inNucleonMom.E();
    tulin[1] = inNucleonMom.Px();
    tulin[2] = inNucleonMom.Py();
    tulin[3] = inNucleonMom.Pz();
    
    for(int i=0; i<4; i++)
      for(int j=0; j<4; j++)
	rulin[i][j] = tulin[i]*tulin[j];
  }

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

  // Calculate LmunuAnumu by iterating over mu and nu
  // In each iteration, add LmunuAnumu to sum if mu=nu, and add
  // LmunuAnumu + LnumuAmunu if mu != nu, since we start nu at mu
  double axx=0.,azz=0.,a0z=0.,a00=0.,axy=0.; // Store elts if fPrintData == true
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
	//if(fPrintData)
	a00 = real(Amunu);
	sum += Lmunu*Amunu;
      }else if(mu == 0 && nu == 3){
	Amunu = 16.0*F1V2*((2.0*rulin[0][3]+tulin[0]*dq)*CN+tulin[3]*q[0])+
	  2.0*q2*xiF2V2*
	  (-4.0*rulin[0][3]/M2-2.0*(dq*tulin[0]+q[0]*tulin[3])/M2-dq*q[0]*(4.0/q2+1.0/M2))+
	  4.0*FA2*((2.0*rulin[0][3]+dq*tulin[0])*CL+q[0]*tulin[3])-
	  (2.0*CL*Fp2*q2+8.0*FA*Fp*CL*M)*dq*q[0]-
	  16.0*F1V*xiF2V*dq*q[0];
	//if(fPrintData) 
	a0z= real(Amunu);
	Anumu = Amunu;
	sum += Lmunu*Anumu + Lnumu*Amunu;
      }else if(mu == 3 && nu == 3){
	Amunu = 16.0*F1V2*(2.0*rulin[3][3]+2.0*dq*tulin[3]-q2/2.0)+ 
	  2.0*q2*xiF2V2*(-4.0-4.0*rulin[3][3]/M2-4.0*dq*tulin[3]/M2-dq2*(4.0/q2+1.0/M2))+
	  4.0*FA2*(2.0*rulin[3][3]+2.0*dq*tulin[3]-(q2/2.0-2.0*CL*M2))-
	  (2.0*CL*Fp2*q2+8.0*FA*Fp*CL*M)*dq2-
	  16.0*F1V*xiF2V*(q2+dq2);
	//if(fPrintData) 
	azz = real(Amunu);
	sum += Lmunu*Amunu;
      }else if(mu ==1 && nu == 1){
	Amunu = 16.0*F1V2*(2.0*rulin[1][1]-q2/2.0)+ 
	  2.0*q2*xiF2V2*(-4.0*CT-4.0*rulin[1][1]/M2) +
	  4.0*FA2*(2.0*rulin[1][1]-(q2/2.0-2.0*CT*M2))-
	  16.0*F1V*xiF2V*CT*q2;
	//if(fPrintData) 
	axx = real(Amunu);
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
	//if(fPrintData) 
	axy = imag(Amunu);
	sum += Lmunu*Anumu+Lnumu*Amunu;
      }else{
	// No RPA corrections to the remaining terms, so A is not r dependent
	if(mu == nu){
	  Amunu = a1*g[mu][nu]+
	    a2*(rulin[mu][nu]+tulin[mu]*q[nu]/2.0+tulin[nu]*q[mu]/2.0)+a4*q[mu]*q[nu];
	  sum += Lmunu*Amunu;
	  //
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
    } // End loop over nu
  } // End loop over mu

  if(fCompareNievesTensors){
    // Compute the coulomb factor
    
    // get tmu
    double tmugev = leptonMom.E() - leptonMom.Mag();
  // Print Q2, form factors, and tensor elts
    ofstream ffstream;
    ffstream.open(fTensorsOutFile, std::ios_base::app);
    if(imU < 0){
    ffstream << -q2 << "\t" << q[0] << "\t" << dq
	     << "\t" << axx << "\t" << azz << "\t" << a0z
	     << "\t" << a00 << "\t" << axy << "\t"
	     << CT << "\t" << CL << "\t" << CN << "\t"
	     << tmugev << "\t" << imU << "\t"
	     << F1V << "\t" << xiF2V << "\t" 
	     << FA << "\t" << Fp << "\t"
	     << tulin[0] << "\t"<< tulin[1] << "\t"
	     << tulin[2] << "\t"<< tulin[3] << "\t"
	     << rulin[0][0]<< "\t"<< rulin[0][1]<< "\t"
	     << rulin[0][2]<< "\t"<< rulin[0][3]<< "\t"
	     << rulin[1][0]<< "\t"<< rulin[1][1]<< "\t"
	     << rulin[1][2]<< "\t"<< rulin[1][3]<< "\t"
	     << rulin[2][0]<< "\t"<< rulin[2][1]<< "\t"
	     << rulin[2][2]<< "\t"<< rulin[2][3]<< "\t"
	     << rulin[3][0]<< "\t"<< rulin[3][1]<< "\t"
	     << rulin[3][2]<< "\t"<< rulin[3][3]<< "\t"
	     << fVc << "\t" << fCoulombFactor << "\t"
	     << fElocal << "\t" << fPlocal << "\t"
	     << fElep << "\t" << fPlep << "\t";
      //<< rhopStored << "\t" << rhonStored << "\t"
      //     << rhoStored << "\t" << rho0Stored << "\t" ;
	   /*<< -q2Orig << "\t"
	     << fKF1 << "\t" << fKF2 << "\t"
	     << fc0 << "\t" << fPrimeStored << "\t" << M << "\t"
	     << fVt << "\t" << fVl << "\t" 
	     << fRelLinReal << "\t" << fRelLinIm << "\t"
	     << fRelLinTotReal << "\t" << fRelLinTotIm << "\t"
	     << fl1 << "\t" << fl2 << "\t" << fl3 << "\t"
	     << q2rellin << "\t" << fkf << "\t" << fef << "\t"
	     << fl2im << "\t" << fl3im << "\t";*/

    ffstream << "\n";
    }
    ffstream.close();
  }else  if(fPrintData){
  // Print Q2, form factors, and tensor elts
    ofstream ffstream;
    ffstream.open("tensors_gevgen_allevents.txt", std::ios_base::app);
    ffstream << -q2 << "\t" << q[0] << "\t" << dq << "\t"
	     << F1V << "\t" << xiF2V << "\t" << FA << "\t" << Fp
	     << "\t" << axx << "\t" << azz << "\t" << a0z
	     << "\t" << a00 << "\t" << axy << "\t"
	     << CT << "\t" << CL << "\t" << CN;
    ffstream << "\n";
    ffstream.close();
  }

  
  // Since the real parts of A and L are both symmetric and the imaginary
  // parts are antisymmetric, the contraction should be real
  if(imag(sum) > kASmallNum)
    LOG("Nieves",pWARN) << "Imaginary part of tensor contraction is nonzero "
			<< "in QEL XSec, real(sum) = " << real(sum) 
			<< "imag(sum) = " << imag(sum);

  return real(sum);
}
//____________________________________________________________________________
void QELPrimaryLeptonGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  //PrimaryLeptonGenerator::ProcessEventRecord(evrec);
  Interaction * interaction = evrec->Summary();
  Kinematics * kine = interaction->KinePtr();

  // Check that a lepton was stored. If not, call the parent method.
  if(!(kine->KVSet(kKVSelTl) && kine->KVSet(kKVSelctl) 
       && kine->KVSet(kKVSelphikq))){
    LOG("LeptonicVertex",pINFO) << "Calling parent method to generate lepton";
    PrimaryLeptonGenerator::ProcessEventRecord(evrec);
    return;
  }
  LOG("LeptonicVertex",pINFO) << "Loading stored lepton";
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
// Generate a lepton using PrimaryLeptonGenerator::ProcessEventRecord()
// and store it by using the KineVar construct with
// the Kinematics class (variables kKVTl, kKVctl, kKVphikq).
// Stored lepton is in the LAB FRAME.
TLorentzVector NievesQELCCPXSec::GenerateOutgoingLepton(const Interaction* interaction,
							TLorentzVector p4v) const{
  // Velocity for an active Lorentz transform taking the final state primary
  // lepton from the [nucleon rest frame] --> [LAB]
  const InitialState & init_state = interaction->InitState();
  const TLorentzVector & pnuc4 = init_state.Tgt().HitNucP4(); //[@LAB]
  TVector3 beta = pnuc4.BoostVector();


  // Neutrino 4p
  p4v.Boost(-1.*beta);                           // v 4p @ Nucleon rest frame

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
  utils::kinematics::WQ2toXY(E,M,W,Q2,x,y);
  double Ev  = p4v.E(); 
  double ml  = interaction->FSPrimLepton()->Mass();
  double ml2 = TMath::Power(ml,2);

  LOG("LeptonicVertex", pDEBUG)
    << "Ev = " << Ev << ", Q2 = " << Q2 << ", y = " << y;

  // Compute the final state primary lepton energy and momentum components
  // along and perpendicular the neutrino direction 
  double El = (1-y)*Ev;
  double plp = El - 0.5*(Q2+ml2)/Ev;                          // p(//)
  double plt = TMath::Sqrt(TMath::Max(0.,El*El-plp*plp-ml2)); // p(-|)

  LOG("LeptonicVertex", pDEBUG)
    << "trying fsl: E = " << El << ", |p//| = " << plp << ", [pT] = " << plt;

  // Randomize transverse components
  RandomGen * rnd = RandomGen::Instance();
  double phi  = 2*kPi * rnd->RndLep().Rndm();
  double pltx = plt * TMath::Cos(phi);
  double plty = plt * TMath::Sin(phi);


  // Take a unit vector along the neutrino direction @ the nucleon rest frame
  TVector3 unit_nudir = p4v.Vect().Unit(); 

  // Rotate lepton momentum vector from the reference frame (x'y'z') where 
  // {z':(neutrino direction), z'x':(theta plane)} to the nucleon rest frame
  TVector3 p3l(pltx,plty,plp);
  p3l.RotateUz(unit_nudir);

  // Lepton 4-momentum in the nucleon rest frame
  TLorentzVector p4l(p3l,El);

  LOG("LeptonicVertex", pINFO)
    << "trying fsl @ NRF: " << utils::print::P4AsString(&p4l);

  // Boost final state primary lepton to the lab frame
  p4l.Boost(beta); // active Lorentz transform

  return p4l;
}
//___________________________________________________________________________
// Auxiliary scalar function for internal integration
//____________________________________________________________________________
utils::gsl::wrap::NievesQELvcrIntegrand::NievesQELvcrIntegrand(
					       double Rcurr, int A, int Z):
ROOT::Math::IBaseFunctionOneDim()
{
  fRcurr = Rcurr;
  fA = A;
  fZ = Z;
}
//____________________________________________________________________________
utils::gsl::wrap::NievesQELvcrIntegrand::~NievesQELvcrIntegrand()
{

}
//____________________________________________________________________________
unsigned int utils::gsl::wrap::NievesQELvcrIntegrand::NDim(void) const
{
  return 1;
}
//____________________________________________________________________________
double utils::gsl::wrap::NievesQELvcrIntegrand::DoEval(double rin) const
{
  double rhop = fZ*nuclear::Density(rin,fA);
  if(rin<fRcurr){
    return rhop*rin*rin/fRcurr;
  }else{
    return rhop*rin;
  }
}
//____________________________________________________________________________
ROOT::Math::IBaseFunctionOneDim * 
  utils::gsl::wrap::NievesQELvcrIntegrand::Clone(void) const
{
  return new utils::gsl::wrap::NievesQELvcrIntegrand(fRcurr,fA,fZ);
}
//____________________________________________________________________________

//____________________________________________________________________________
//
// NOTE: THE REMAINING THREE METHODS ARE FOR TESTING PURPOSES ONLY
//
void NievesQELCCPXSec::CompareNievesTensors(const Interaction* in) 
  const {
  //
  // For this comparison to agree with Nieves' fortran code, the form factors
  // should be set as the Dipole form factors class.
  //
  std::cout << "Interating Kinematics and printing tensor elts" << std::endl;
  Interaction * interaction = new Interaction(*in); // copy in
  
  // These should be set by the user- eg read these from file
  double ein = 1.0;
  double ctl = 0.0;
  double rmaxfrac = 0.25;

  if(fRPA)
    fTensorsOutFile = "gen.RPA_E1_ctl0_r0.5";
  else
    fTensorsOutFile = "gen.noRPA_E1_ctl0_r0.5";

  // Calculate radius
  // CARBON
  bool klave = true;
  double rp = 1.692;
  double ap = 1.082;
  double rn = 1.692;
  double an = 1.082;
  double rmax;
  if(!klave)
    rmax = TMath::Max(rp,rn) + 9.25*TMath::Max(ap,an);
  else
    rmax = TMath::Sqrt(20.0)*TMath::Max(rp,rn);
  double r = rmax *  rmaxfrac;

  // Relevant objects and parameters
  //const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();
  const Target & target = init_state.Tgt();

  // Paramters required for LmunuAnumu
  double M  = target.HitNucMass();
  double ml = interaction->FSPrimLepton()->Mass();
  bool is_neutrino = pdg::IsNeutrino(init_state.ProbePdg());
  bool tgtIsNucleus = target.IsNucleus();
  int tgt_pdgc = target.HitNucPdg();
  int A = target.A();
  int Z = target.Z();
  int N = target.N();
  bool hitNucIsProton = pdg::IsProton(target.HitNucPdg());

  // Iterate over lepton energy (which then affects q, which is passed to
  // LmunuAnumu using in and out NucleonMom
  double delta = (ein-0.025)/100.0;
  for(int it=0;it<100;it++){
    double tmu = it*delta;
    double eout = ml + tmu;
    double pout = TMath::Sqrt(eout*eout-ml*ml);

    double pin = TMath::Sqrt(ein*ein); // Assume massless neutrinos

    double qvalue = .016827; // Units of GeV
    double q0 = ein-eout-qvalue;
    double dq = TMath::Sqrt(pin*pin+pout*pout-2.0*ctl*pin*pout);
    double q2 = q0*q0-dq*dq;
    interaction->KinePtr()->SetQ2(-q2);
    //q2Orig = interaction->KinePtr()->q2();

    // When this method is called, inNucleonMom and outNucleonMom are 
    // only used to calulate 
    // q = outNucleonMom - inNucleonMom. I can thus provide the calculated
    // values using outNucleonMom and inNucleonMom and putting q in the
    // z direction, as Nieves does in his paper
    TLorentzVector inNucleonMom(0,0,0,0);
    TLorentzVector outNucleonMom(0,0,dq,q0);
    std::cout << "E = " << outNucleonMom.E() << ", p = " << outNucleonMom.Vect().Mag()
	      << ", q2 = " << outNucleonMom.Mag2() << ", q2 = " << q2 << std::endl;

    // neutrinoMom and leptonMom only directly affect the leptonic tensor, which
    // we are not calculating now. Set to large values so we see if they
    // do make a difference (which would be bad).
    TLorentzVector * neutrinoMom = new TLorentzVector(999999,999999,999999,999999);
    TLorentzVector leptonMom(0,0,pout,eout);

  if(fCoulomb){
    // Coulomb potential
    double Vc = vcr(& target, r);
    fVc = Vc;

    // Outgoing lepton energy and momentum including coulomb potential
    int sign = is_neutrino ? 1 : -1;
    double El = leptonMom.E();
    double ElLocal = El - sign*Vc;
    if(ElLocal - ml <= 0.0){
      LOG("Nieves",pINFO) << "Event should be rejected. Coulomb effects "
			    << "push kinematics below threshold";
      return;
    }
    double plLocal = TMath::Sqrt(ElLocal*ElLocal-ml*ml);

    // Correction factor
    double coulombFactor= plLocal*ElLocal/leptonMom.Vect().Mag()/El;
    fCoulombFactor = coulombFactor;
    fElocal = ElLocal;
    fPlocal = plLocal;
    fElep = El;
    fPlep = leptonMom.Vect().Mag();
  }
    
    
    fFormFactors.Calculate(interaction);
    LmunuAnumu(*neutrinoMom,inNucleonMom,leptonMom,outNucleonMom,
	       M,r,is_neutrino,tgtIsNucleus,tgt_pdgc,A,Z,N,hitNucIsProton);
  }
  std::cout << "Done Interating Kinematics and printing tensor elts" << std::endl;
  return;
}
//____________________________________________________________________________
