//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GBAnisotropyMisoriBase.h"

#include "EulerAngleProvider.h"
#include "GrainTracker.h"
#include "MisorientationAngleCalculator.h"

// Forward Declarations

/**
 * Function[kappa, gamma, m, L] = parameters (sigma, mob, w_GB, sigma0)
 * Parameter determination method is elaborated in Phys. Rev. B, 78(2), 024113, 2008, by N. Moelans
 * Thanks to Prof. Moelans for the explanation of her paper.
 */
class GBAnisotropyMisori : public GBAnisotropyMisoriBase
{
public:
  static InputParameters validParams();

  GBAnisotropyMisori(const InputParameters & parameters);

protected:
  // Set sigma_ij, mob_ij, Q_ij for specific grain boundaries
  virtual void computeGBProperties() override;

  // calculated GB energy based on the the Read-Shockley
  virtual Real calculatedGBEnergy(const MisorientationAngleData & misori_s);

  // calculated GB mobility based on the sigmoidal law
  virtual Real calculatedGBMobility(const MisorientationAngleData & misori_s);

  // used to store orientation structure, including misorientation angle, istwinnig, twinning type;
  MisorientationAngleData _misori_s;
  
  const Real _GBsigma_HAGB;
  const Real _GBmob_HAGB;

  const GrainTracker & _grain_tracker;
  const EulerAngleProvider & _euler; 

  const bool _gb_energy_anisotropy;
  const bool _gb_mobility_anisotropy;    

  MaterialProperty<Real> & _misori_angle;
  MaterialProperty<Real> & _twinning_type;
};