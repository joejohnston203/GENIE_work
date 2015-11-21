//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Joe Johnston, University of Pittsburgh

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <sstream>
#include <cstdlib>
#include <TMath.h>

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Messenger/Messenger.h"
#include "Nuclear/LFGMBodekRitchie.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Numerical/RandomGen.h"
#include "Utils/NuclearUtils.h"

using std::ostringstream;
using namespace genie;
using namespace genie::constants;
using namespace genie::utils;

//____________________________________________________________________________
LFGMBodekRitchie::LFGMBodekRitchie() :
NuclearModelI("genie::LFGMBodekRitchie")
{

}
//____________________________________________________________________________
LFGMBodekRitchie::LFGMBodekRitchie(string config) :
NuclearModelI("genie::LFGMBodekRitchie", config)
{

}
//____________________________________________________________________________
LFGMBodekRitchie::~LFGMBodekRitchie()
{
  map<std::pair<string,double> , TH1D*>::iterator iter =fProbDistroMap.begin();
  for( ; iter != fProbDistroMap.begin(); ++iter) {
    TH1D * hst = iter->second;
    if(hst) {
      delete hst;
      hst=0;
    }
  }
  fProbDistroMap.clear();
}
//____________________________________________________________________________
bool LFGMBodekRitchie::GenerateNucleon(const Target & target) const
{
  assert(target.HitNucIsSet());

  fCurrRemovalEnergy = 0;
  fCurrMomentum.SetXYZ(0,0,0);

  //-- set fermi momentum vector
  //
  TH1D * prob = this->ProbDistro(target);
  if(!prob) {
    LOG("LFGMBodekRitchie", pNOTICE)
              << "Null nucleon momentum probability distribution";
    exit(1);
  }
  double p = prob->GetRandom();
  LOG("LFGMBodekRitchie", pINFO) << "|p,nucleon| = " << p;

  RandomGen * rnd = RandomGen::Instance();

  double costheta = -1. + 2. * rnd->RndGen().Rndm();
  double sintheta = TMath::Sqrt(1.-costheta*costheta);
  double fi       = 2 * kPi * rnd->RndGen().Rndm();
  double cosfi    = TMath::Cos(fi);
  double sinfi    = TMath::Sin(fi);

  double px = p*sintheta*cosfi;
  double py = p*sintheta*sinfi;
  double pz = p*costheta;  

  fCurrMomentum.SetXYZ(px,py,pz);

  //-- set removal energy 
  //
  int Z = target.Z();
  map<int,double>::const_iterator it = fNucRmvE.find(Z);
  if(it != fNucRmvE.end()) fCurrRemovalEnergy = it->second;
  else fCurrRemovalEnergy = nuclear::BindEnergyPerNucleon(target);

  return true;
}
//____________________________________________________________________________
double LFGMBodekRitchie::Prob(double p, double w, const Target & target) const
{
  if(w<0) {
    TH1D * prob = this->ProbDistro(target);
    int bin = prob->FindBin(p);
    double y  = prob->GetBinContent(bin);
    double dx = prob->GetBinWidth(bin);
    double p  = y*dx;
    return p;
  }
  return 1;
}
//____________________________________________________________________________
TH1D * LFGMBodekRitchie::ProbDistro(const Target & target) const
{
  // Do not store computed values, because the radius will be different
  // for nearly every event
  //-- return stored /if already computed/
  //std::pair <string,double> targetStr (target.AsString(),
  //				       target.HitNucRadius());
  //map<std::pair<string,double>, TH1D*>::iterator it = 
  //  fProbDistroMap.find(targetStr);
  //if(it != fProbDistroMap.end()) return it->second;

  LOG("LFGMBodekRitchie", pNOTICE)
             << "Computing P = f(p_nucleon) for: " << target.AsString();
  LOG("LFGMBodekRitchie", pNOTICE)
               << "P(cut-off) = " << fPCutOff << ", P(max) = " << fPMax;

  //-- get information for the nuclear target
  int nucleon_pdgc = target.HitNucPdg();
  assert(pdg::IsProton(nucleon_pdgc) || pdg::IsNeutron(nucleon_pdgc));
  int A = target.A();

  assert(target.HitNucIsSet());
  bool is_p = pdg::IsProton(nucleon_pdgc);
  double numNuc = (is_p) ? (double)target.Z():(double)target.N();

  // Calculate Fermi Momentum using LFG model
  double r = target.GetHitNucRadius(), hbarc = .1973269602;
  double KF= TMath::Power(3*kPi2*numNuc*genie::utils::nuclear::Density(r,A),
			    1.0/3.0) *hbarc;

  LOG("LFGMBodekRitchie",pDEBUG) << "r = " << r << ",KF LFG = " << KF;

  double a  = 2.0;
  double C  = 4. * kPi * TMath::Power(KF,3) / 3.;
  // Do not include nucleon correlation tail
  //double R  = 1. / (1.- KF/4.);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("LFGMBodekRitchie", pDEBUG) << "a  = " << a;
  LOG("LFGMBodekRitchie", pDEBUG) << "C  = " << C;
  //LOG("LFGMBodekRitchie", pDEBUG) << "R  = " << R;
#endif

  //-- create the probability distribution

  int npbins = (int) (1000*fPMax);
  TH1D * prob = new TH1D("", "", npbins, 0, fPMax);
  prob->SetDirectory(0);

  double dp = fPMax / (npbins-1);
  double iC = (C>0) ? 1./C : 0.;
  double kfa_pi_2 = TMath::Power(KF*a/kPi,2);

  for(int i = 0; i < npbins; i++) {
     double p  = i * dp;
     double p2 = TMath::Power(p,2);

     // calculate |phi(p)|^2
     double phi2 = 0;
     if (p <= KF)
        phi2 = iC * (1. - 6.*kfa_pi_2);
     // Do not include nucleon correlation tail
     //else if ( p > KF && p < fPCutOff)
     //   phi2 = iC * (2*R*kfa_pi_2*TMath::Power(KF/p,4.));

     // calculate probability density : dProbability/dp
     double dP_dp = 4*kPi * p2 * phi2;
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("LFGMBodekRitchie", pDEBUG) << "p = " << p << ", dP/dp = " << dP_dp;
#endif
     prob->Fill(p, dP_dp);
  }

  //-- normalize the probability distribution
  prob->Scale( 1.0 / prob->Integral("width") );

  //-- store
  //LOG("LFGMBodekRitchie",pDEBUG) << "Target name as string = " << 
  //  target.AsString();
  //std::pair <string,double> targetPair (target.AsString(),
  //					  target.GetHitNucRadius());
  //fProbDistroMap.insert(
  //	map<pair<string,double>, TH1D*>::value_type(targetPair, prob));

  return prob; 
}
//____________________________________________________________________________
void LFGMBodekRitchie::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void LFGMBodekRitchie::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void LFGMBodekRitchie::LoadConfig(void)
{
  fPMax    = fConfig->GetDoubleDef ("MomentumMax", 1.0);
  assert(fPMax > 0);

  // Load removal energy for specific nuclei from either the algorithm's
  // configuration file or the UserPhysicsOptions file.
  // If none is used use Wapstra's semi-empirical formula.
  //
  for(int Z=1; Z<140; Z++) {
    for(int A=Z; A<3*Z; A++) {
      ostringstream key, gckey;
      int pdgc = pdg::IonPdgCode(A,Z);
      gckey << "RFG-NucRemovalE@Pdg=" << pdgc;
      key   << "NucRemovalE@Pdg="     << pdgc;
      RgKey gcrgkey = gckey.str();
      RgKey rgkey   = key.str();
      if (this->GetConfig().Exists(rgkey) || gc->Exists(gcrgkey)) {
        double eb = fConfig->GetDoubleDef(rgkey, gc->GetDouble(gcrgkey));
        eb = TMath::Max(eb, 0.);
        LOG("LFGMBodekRitchie", pINFO)
          << "Nucleus: " << pdgc << " -> using Eb =  " << eb << " GeV";
        fNucRmvE.insert(map<int,double>::value_type(Z,eb));
      }
    }
  }

  // LFG does not use a Fermi Momentum Table
  //AlgConfigPool * confp = AlgConfigPool::Instance();
  //const Registry * gc = confp->GlobalParameterList();
  //fKFTable = fConfig->GetStringDef ("FermiMomentumTable", 
  //gc->GetString("FermiMomentumTable"));
  //
  // Do not include nucleon correlation tail for LFG
  //fPCutOff = fConfig->GetDoubleDef ("MomentumCutOff",  
  //gc->GetDouble("RFG-MomentumCutOff"));
  //assert(fPMax > 0 && fPCutOff > 0 && fPCutOff < fPMax);
}
//____________________________________________________________________________
