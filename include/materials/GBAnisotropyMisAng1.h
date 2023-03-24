#pragma once
#include "GBAnisotropyMisAngBase.h"

/**
 * This version of GB energy and GB mobility anisotropy is introduced by the Read-Shockley  
 * law and the sigmoidal law from J. Gao, et al, Acta Materialia. 223 (2022) 117491.
 * And set the key parameters (_delta_theta_HAGB, _GBsigma_HAGB, _GBmob_HAGB) 
 * directly through the input file
**/

class GBAnisotropyMisAng1 : public GBAnisotropyMisAngBase
{
public:
  static InputParameters validParams();  
  GBAnisotropyMisAng1(const InputParameters & parameters);

protected:
  // calculated GB energy based on the the Read-Shockley
  virtual Real calculatedGBEnergy(const MisorientationAngleData & _misori_s) override; 

  // calculated GB mobility based on the sigmoidal law
  virtual Real calculatedGBMobility(const MisorientationAngleData & _misori_s) override;  

  const Real _delta_theta_HAGB;
  const Real _GBsigma_HAGB;
  const Real _GBmob_HAGB;
};