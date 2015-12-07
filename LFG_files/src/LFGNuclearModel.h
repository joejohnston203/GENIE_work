//____________________________________________________________________________
/*!

\class    genie::LFGNuclearModel

\brief    local Fermi gas model. Implements the NuclearModelI 
          interface.

\ref      

\author   Joe Johnston, Steven Dytman

\created  December 2015

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _LFG_NUCLEAR_MODEL_H_
#define _LFG_NUCLEAR_MODEL_H_

#include <map>

#include <TH1D.h>
#include "Nuclear/NuclearModelI.h"

using std::map;

namespace genie {

class LFGNuclearModel : public NuclearModelI {

public:
  LFGNuclearModel();
  LFGNuclearModel(string config);
  virtual ~LFGNuclearModel();

  //-- add methods to be called with a nucleon radius;
  bool   GenerateNucleon (const Target & t, double r) const;
  double Prob            (double p, double w, const Target & t, double r) const;

  //-- implement the NuclearModelI interface
  bool           GenerateNucleon (const Target & t) const {
    return GenerateNucleon(t, 0.0);
  }
  double         Prob (double p, double w, const Target & t) const{
    return Prob(p,w,t,0.0);
  }
  NuclearModel_t ModelType       (const Target &) const 
  { 
    return kNucmLocalFermiGas; 
  }

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:
  void   LoadConfig (void);
  TH1D * ProbDistro (const Target & t, double r) const;

  map<int, double> fNucRmvE;

  double fPMax;
};

}         // genie namespace
#endif    // _LFG_NUCLEAR_MODEL_H_

