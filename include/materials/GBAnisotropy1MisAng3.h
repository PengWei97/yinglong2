#pragma once
#include "GBAnisotropy1MisAngBase.h"

// GBAnisotropy1MisAng3 -- This version of GB energy and GB mobility anisotropy is introduced by the Read-Shockley law and the sigmoidal law from J. Gao, et al, Acta Materialia. 223 (2022) 117491.

class GBAnisotropy1MisAng3 : public GBAnisotropy1MisAngBase
{
public:
  static InputParameters validParams();  
  GBAnisotropy1MisAng3(const InputParameters & parameters);

protected:
  // calculated GB energy based on the the Read-Shockley
  virtual Real calculatedGBEnergy(const misoriAngle_isTwining & misori_gbType); 

  // calculated GB mobility based on the sigmoidal law
  virtual Real calculatedGBMobility(const misoriAngle_isTwining & misori_gbType);  

  const Real _delta_theta_HAB;
  const Real _GBsigma_HAB;
  const Real _GBmob_HAB;
  const Real _TT1_sigma;
  const Real _CT1_sigma;
};