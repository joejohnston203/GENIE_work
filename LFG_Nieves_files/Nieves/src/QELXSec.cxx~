//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - May 04, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Jan 18, 2008 - CA
   Add protection against unphysical Q2 limits
 @ Sep 07, 2009 - CA
   Integrated with GNU Numerical Library (GSL) via ROOT's MathMore library.
 @ Nov 30, 2009 - CA
   Added option to integrate over the hit nucleon momentum distribution.
 @ Jan 26, 2016 - JJ
   Choose a radius for the struck nucleon before choosing the momentum, and
   save that radius for use by the XSec method
*/
//____________________________________________________________________________

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "Conventions/Units.h"
#include "Conventions/KineVar.h"
#include "Conventions/RefFrame.h"
#include "CrossSections/QELXSec.h"
#include "CrossSections/GSLXSecFunc.h"
#include "Messenger/Messenger.h"
#include "Nuclear/NuclearModelI.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "Utils/KineUtils.h"
#include "Utils/NuclearUtils.h"
#include "Utils/Range1.h"
#include "Utils/GSLUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
QELXSec::QELXSec() :
XSecIntegratorI("genie::QELXSec")
{

}
//____________________________________________________________________________
QELXSec::QELXSec(string config) :
XSecIntegratorI("genie::QELXSec", config)
{

}
//____________________________________________________________________________
QELXSec::~QELXSec()
{

}
//____________________________________________________________________________
double QELXSec::Integrate(
                  const XSecAlgorithmI * model, const Interaction * in) const
{
  if(! model->ValidProcess(in)) return 0.;

  bool nuclear_target = in->InitState().Tgt().IsNucleus();

  if(!nuclear_target || !fDoAvgOverNucleonMomentum) {
    return this->IntegrateOnce(model,in);
  }

  double E = in->InitState().ProbeE(kRfHitNucRest);
  // Also average over positions in the nucleus if model == lfg
  bool lfg = (fNuclModel && fNuclModel->ModelType(Target()) == kNucmLocalFermiGas);
  if(E < fEnergyCutOff || lfg) {
     // clone the input interaction so as to tweak the
     // hit nucleon 4-momentum in the averaging loop
     Interaction in_curr(*in);

     // hit target
     const Target & tgt = in_curr.InitState().Tgt();

     // get nuclear masses (init & final state nucleus)
     int nucleon_pdgc = tgt.HitNucPdg();
     bool is_p = pdg::IsProton(nucleon_pdgc);
     int Zi = tgt.Z();
     int Ai = tgt.A();
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
     for(int inuc=0; inuc<nnuc; inuc++) {
       // Generate a radius, as in VertexGenerator.cxx
       double r = GenerateRadius(in_curr);
       // Save the radius in the interaction for use by the XSec method
       in_curr.KinePtr()->SetKV(kKVSelRad,r);

       fNuclModel->GenerateNucleon(tgt, r);
       TVector3 p3N = fNuclModel->Momentum3();
       double   EN  = Mi - TMath::Sqrt(p3N.Mag2() + Mf*Mf);

       TLorentzVector* p4N = in_curr.InitStatePtr()->TgtPtr()->HitNucP4Ptr();
       p4N->SetPx (p3N.Px());
       p4N->SetPy (p3N.Py());
       p4N->SetPz (p3N.Pz());
       p4N->SetE  (EN);

       double xsec = this->IntegrateOnce(model,&in_curr);
       xsec_sum += xsec;
     }
     double xsec_avg = xsec_sum / nnuc;
     return xsec_avg;

  } else {
    return this->IntegrateOnce(model,in);
  }

  return 0;
}
//____________________________________________________________________________
double QELXSec::IntegrateOnce(
			      const XSecAlgorithmI * model, const Interaction * in) const
{
  const KPhaseSpace & kps = in->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
     LOG("QELXSec", pDEBUG)  << "*** Below energy threshold";
     return 0;
  }
  Range1D_t rQ2 = kps.Limits(kKVQ2);
  if(rQ2.min<0 || rQ2.max<0) return 0;
  LOG("QELXSec", pDEBUG) 
          << "Q2 integration range = (" << rQ2.min << ", " << rQ2.max << ")";

  Interaction * interaction = new Interaction(*in);
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  ROOT::Math::IBaseFunctionOneDim * func = new 
      utils::gsl::dXSec_dQ2_E(model, interaction);
  ROOT::Math::IntegrationOneDim::Type ig_type = 
      utils::gsl::Integration1DimTypeFromString(fGSLIntgType);
  
  double abstol = 1; //We mostly care about relative tolerance
  ROOT::Math::Integrator ig(*func,ig_type,abstol,fGSLRelTol,fGSLMaxEval);
  double xsec = ig.Integral(rQ2.min, rQ2.max) * (1E-38 * units::cm2);
     
  //LOG("QELXSec", pDEBUG) << "XSec[QEL] (E = " << E << ") = " << xsec;

  delete func;
  delete interaction;

  return xsec;
}
//____________________________________________________________________________
double QELXSec::GenerateRadius(Interaction interaction) const{
  RandomGen * rnd = RandomGen::Instance();
  TVector3 vtx(999999.,999999.,999999.);

  bool uniform = fVtxGenMethod==0;
  bool realistic = !uniform;
  
  double A = interaction.InitStatePtr()->TgtPtr()->A();
  double R = fR0 * TMath::Power(A, 1./3.);

  // ** Set the interaction vertex
  // ** using the requested method: either using a realistic nuclear 
  // ** density or by setting it uniformly within the nucleus
  if(realistic) {
    // Generate the vertex using a realistic nuclear density
    //
    // get inputs to the rejection method
    double ymax = -1;
    double rmax = 3*R;
    double dr   = R/40.;
    for(double r = 0; r < rmax; r+=dr) { 
      ymax = TMath::Max(ymax, r*r * utils::nuclear::Density(r,(int)A)); 
    }
    ymax *= 1.2;
    
    // select a vertex using the rejection method
    unsigned int iter = 0;
    while(1) {
      iter++;
      
      // return 0. if it hasn't find a solution after many attempts
      if(iter > controls::kRjMaxIterations) {
	LOG("QELXSec", pWARN)
	  << "Couldn't generate a vertex position after " << iter << " iterations. Setting radius = 0.0";
	return 0.;
      }
      
      double r = rmax * rnd->RndFsi().Rndm();
      double t = ymax * rnd->RndFsi().Rndm();
      double y = r*r * utils::nuclear::Density(r,(int)A);
      if(y > ymax) {
	LOG("Vtx", pERROR)
	  << "y = " << y << " > ymax = " << ymax 
	  << " for r = " << r << ", A = " << A;
      }
      bool accept = (t < y);
      if(accept) {
	return r;
      }
    }//w(1)
  } //use density?
  if(uniform) {
    // Generate the vertex uniformly within the nucleus
    //
    LOG("Vtx", pINFO) 
      << "Generating intranuclear vertex uniformly in volume";
    while(vtx.Mag() > R) {
      vtx.SetX(-R + 2*R * rnd->RndFsi().Rndm());
      vtx.SetY(-R + 2*R * rnd->RndFsi().Rndm());
      vtx.SetZ(-R + 2*R * rnd->RndFsi().Rndm());
    }
    return vtx.Mag();
  }// uniform? 
  return 0.0;
}
//____________________________________________________________________________
void QELXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void QELXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void QELXSec::LoadConfig(void)
{
  // Get GSL integration type & relative tolerance
  fGSLIntgType = fConfig->GetStringDef("gsl-integration-type", "adaptive");
  fGSLRelTol   = fConfig->GetDoubleDef("gsl-relative-tolerance", 0.01);
  fGSLMaxEval  = (unsigned int) fConfig->GetIntDef ("gsl-max-eval", 100000);

  fDoAvgOverNucleonMomentum =
     fConfig->GetBoolDef("AverageOverNucleonMomentum", false);

  fNuclModel    = 0;
  fEnergyCutOff = 0.;

  if(fDoAvgOverNucleonMomentum) {
    // Get nuclear model
    RgKey nuclkey = "NuclearModel";
    fNuclModel = dynamic_cast<const NuclearModelI *> (this->SubAlg(nuclkey));
    assert(fNuclModel);
   
    // Get averaging cutoff energy
    fEnergyCutOff = 
      fConfig->GetDoubleDef("NuclearInfluenceCutoffEnergy", 2.0);
  }

  // Get parameters for generating a position inside the nucleus
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  fVtxGenMethod = fConfig->GetIntDef (
           "VtxGenerationMethod", gc->GetInt("NUCL-VtxGenerationMethod")); 

  fR0 = fConfig->GetDoubleDef (
           "R0", gc->GetDouble("NUCL-R0")); // fm
}
//____________________________________________________________________________

