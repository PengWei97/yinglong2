#pragma once
#include "GBAnisotropyMisAngBase.h"

/**
 * This version of GB anisotropy is introduced by the Read-Shockley  
 * law and the sigmoidal law from J. Gao, et al, Acta Materialia. 223 (2022) 117491.
 * And set the key parameters (_delta_theta_HAGB, _GBsigma_HAGB, _GBmob_HAGB) 
 * directly through the input file
 * Meanwhile, add the function to identify and set GB energy of twin boundaries.
**/

class GBAnisotropyMisAng2 : public GBAnisotropyMisAngBase
{
public:
  static InputParameters validParams();  
  GBAnisotropyMisAng2(const InputParameters & parameters);

protected:
  // calculated GB energy based on the the Read-Shockley
  virtual Real calculatedGBEnergy(const MisorientationAngleData & _misori_s) override; 

  // TODO 如何有效设定孪晶界的mobility
  // calculated GB mobility based on the sigmoidal law
  virtual Real calculatedGBMobility(const MisorientationAngleData & _misori_s) override;  

  const Real _delta_theta_HAGB;
  const Real _GBsigma_HAGB;
  const Real _GBmob_HAGB;

  const Real _TT1_sigma;
  const Real _CT1_sigma;
};