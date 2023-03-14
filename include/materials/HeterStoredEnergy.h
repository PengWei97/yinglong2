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
 * Calculates the Deformation Energy associated with a specific dislocation density from EBSD data file.
 * The rest of parameters are the same as in the grain growth model
 */
class HeterStoredEnergy : public Material
{
public:
  static InputParameters validParams();

  HeterStoredEnergy(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  const unsigned int _op_num; // total number of grains
  const std::vector<const VariableValue *> _vals; // order parameter values

  const Real _length_scale;
  const Real _time_scale;
  const Real _JtoeV;
  const Real _Elas_Mod; // the elastic modulus
  const Real _Burg_vec; // the Length of Burger's Vector
 
  const GrainTrackerInterface & _grain_tracker; // Grain tracker object
  const EBSDReader & _GNDs_provider;
  
  MaterialProperty<Real> & _beta;
  MaterialProperty<Real> & _rho_eff; // the average/effective dislocation density
  MaterialProperty<Real> & _Def_Eng; // the deformation energy
};
