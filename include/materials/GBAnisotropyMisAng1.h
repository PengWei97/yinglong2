#pragma once
#include "GBAnisotropy1MisAngBase.h"

/**
 * This version of GB energy and GB mobility anisotropy is introduced by the Read-Shockley  
 * law and the sigmoidal law from J. Gao, et al, Acta Materialia. 223 (2022) 117491.
 * And set the key parameters (_delta_theta_HAGB, _GBsigma_HAGB, _GBmob_HAGB) 
 * directly through the input file
**/

class GBAnisotropy1MisAng : public GBAnisotropy1MisAngBase
{
public:
  static InputParameters validParams();  
  GBAnisotropy1MisAng(const InputParameters & parameters);

protected:
  // calculated GB energy based on the the Read-Shockley
  virtual Real calculatedGBEnergy(const MisorientationAngleData & s_misorientation_angle) override; 

  // calculated GB mobility based on the sigmoidal law
  virtual Real calculatedGBMobility(const MisorientationAngleData & s_misorientation_angle) override;  

  const Real _delta_theta_HAGB;
  const Real _GBsigma_HAGB;
  const Real _GBmob_HAGB;
};