//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "EBSDReader.h"

// Forward Declarations
class GrainTrackerInterface;

/**
 * Calculates the Deformation Energy associated with a specific geometrically necessary 
 * dislocation density (GNDs) from EBSD data file.
 * MaterialProperty [_beta, _rho_eff] are used to calculated the stored energy.
 */
class DeformedGrainEBSDMaterial : public Material
{
public:
  static InputParameters validParams();

  DeformedGrainEBSDMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  virtual Real getGNDsFromEBSD(const unsigned int & grain_id);

  const unsigned int _op_num; // total number of grains
  const std::vector<const VariableValue *> _vals; // order parameter values

  const Real _length_scale;
  const Real _time_scale;
  const Real _Elas_Mod; // the elastic modulus
  const Real _Burg_vec; // the Length of Burger's Vector
  const Real _stored_factor;
  const Real _JtoeV; // Joule to eV conversion
  const bool _concurrent_recovery;
 
  MaterialProperty<Real> & _beta;
  MaterialProperty<Real> & _rho_eff; // the average effective dislocation density
 
  const GrainTrackerInterface & _grain_tracker; // Grain tracker object
  const EBSDReader & _GNDs_provider;
};
