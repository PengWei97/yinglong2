//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DeformedGrainGG.h"
#include "GrainTrackerInterface.h"

registerMooseObject("yinglongApp", DeformedGrainGG);

InputParameters
DeformedGrainGG::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  params.addParam<Real>("length_scale", 1.0e-9, "Length scale in m, where default is nm");
  params.addParam<Real>("time_scale", 1.0e-6, "Time scale in sec, where default is micro sec");
  params.addParam<Real>("Elas_Mod", 2.50e10, "Elastic Modulus in J/m^3");
  params.addParam<Real>("Burg_vec", 3.0e-10, "Length of Burger Vector in m");
  params.addRequiredParam<UserObjectName>("grain_tracker",
                                          "The GrainTracker UserObject to get values from.");
  params.addRequiredParam<UserObjectName>("GNDs_provider", "Geometric necessary dislocation densities (GNDs) provider for EBSD reader");
  return params;
}

DeformedGrainGG::DeformedGrainGG(const InputParameters & parameters)
  : Material(parameters),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _length_scale(getParam<Real>("length_scale")),
    _time_scale(getParam<Real>("time_scale")),
    _Elas_Mod(getParam<Real>("Elas_Mod")),
    _Burg_vec(getParam<Real>("Burg_vec")),
    _beta(declareProperty<Real>("beta")),
    _rho_eff(declareProperty<Real>("rho_eff")),
    _Def_Eng(declareProperty<Real>("Def_Eng")),
    _grain_tracker(getUserObject<GrainTrackerInterface>("grain_tracker")),
    _GNDs_provider(getUserObject<EBSDReader>("GNDs_provider")),    
    _JtoeV(6.24150974e18) // Joule to eV conversion
{
  if (_op_num == 0)
    paramError("op_num", "Model requires op_num > 0");
}

void
DeformedGrainGG::computeQpProperties()
{
  Real SumEtai2 = 0.0;
  for (unsigned int i = 0; i < _op_num; ++i)
    SumEtai2 += (*_vals[i])[_qp] * (*_vals[i])[_qp];

  // calculate effective dislocation density
  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());

  Real rho_i;
  Real rho0 = 0.0;
  _rho_eff[_qp] = 0.0;
  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    if (op_to_grains[op_index] == FeatureFloodCount::invalid_id)
      continue;

    auto grain_id = op_to_grains[op_index];

    rho_i = _GNDs_provider.getAvgData(grain_id)._custom[0]; // GNDs for each grain, 1/m^2 
    rho0 += rho_i * (*_vals[op_index])[_qp] * (*_vals[op_index])[_qp];
  }

  _rho_eff[_qp] = rho0 / SumEtai2;

  _beta[_qp] = 0.5 * _Elas_Mod * _Burg_vec * _Burg_vec * _JtoeV * _length_scale;

  // Compute the deformation energy
  _Def_Eng[_qp] = _beta[_qp] * _rho_eff[_qp];
}