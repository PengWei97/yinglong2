#pragma once
#include "GBAnisotropy1MisAngBase.h"
#include "LinearInterpolation.h"

// GBAnisotropy1MisAng5 -- This version of the GB energy can be introduced based on an interpolation function, and additionally assumes that the GB mobility is proportional to the GB energy. Meanwhile, add the function to identify and set the energy and mobility of twin boundaries.

class GBAnisotropy1MisAng5 : public GBAnisotropy1MisAngBase
{
public:
  static InputParameters validParams();  
  GBAnisotropy1MisAng5(const InputParameters & parameters);

protected:
  // calculated GB energy based on the the Read-Shockley
  virtual Real calculatedGBEnergy(const misoriAngle_isTwining & misori_gbType); 

  // calculated GB mobility based on the sigmoidal law
  virtual Real calculatedGBMobility(const misoriAngle_isTwining & misori_gbType);  

  const Real _delta_theta_HAB;
  const Real _GBsigma_LAGB;
  const Real _GBsigma_HAGB;
  const Real _GBmob_HAGB;
  const Real _TT1_sigma;
  const Real _CT1_sigma;

  LinearInterpolation _piecewise_func_sigma;
};