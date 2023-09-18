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
 * GBAnisotropyMisoriBase is created based on GrainTracker for real-time acquisition of
 * sigma_ij, mob_ij based on misorientation angle.
 * function 1: GB anisotropy based classic theroies
 * function 2: Low-energy and low-mobility properties of twin interfaces
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
  
  // for HCP_Ti
  const Real _TT1_sigma;
  const Real _CT1_sigma;
  const Real _TT1_mob;
  const Real _CT1_mob;

  // for FCC_Ni
  const Real _Sigma3_sigma;
  const Real _Sigma9_sigma;
  const Real _Sigma3_mob;
  const Real _Sigma9_mob;

  const GrainTracker & _grain_tracker;
  const EulerAngleProvider & _euler; 

  const bool _gb_energy_anisotropy;
  const bool _gb_mobility_anisotropy;    

  MaterialProperty<Real> & _misori_angle;
  MaterialProperty<Real> & _twinning_type;
};