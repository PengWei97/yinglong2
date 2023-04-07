//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ACSDDFMPoly.h"
#include "Material.h"
#include "GrainTrackerInterface.h"

registerMooseObject("yinglongApp", ACSDDFMPoly);

InputParameters
ACSDDFMPoly::validParams()
{
  InputParameters params = ACBulk<Real>::validParams();
  params.addClassDescription("Stored Energy contribution to grain growth");
  params.addRequiredCoupledVar("v", "Array of coupled variable names");
  params.addRequiredParam<Real>("length_scale", "Length scale in m");
  params.addRequiredParam<UserObjectName>("grain_tracker",
                                          "The GrainTracker UserObject to get values from.");
  params.addRequiredParam<UserObjectName>("GNDs_provider", "Geometric necessary dislocation densities (GNDs) provider for EBSD reader");                                          
  params.addRequiredParam<unsigned int>("op_index", "The index for the current order parameter");
  return params;
}

ACSDDFMPoly::ACSDDFMPoly(const InputParameters & parameters)
  : ACBulk<Real>(parameters),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _vals_var(coupledIndices("v")),
    _length_scale(getParam<Real>("length_scale")),
    _beta(getMaterialProperty<Real>("beta")),
    _rho_eff(getMaterialProperty<Real>("rho_eff")),
    _grain_tracker(getUserObject<GrainTrackerInterface>("grain_tracker")),
    _GNDs_provider(getUserObject<EBSDReader>("GNDs_provider")),  
    _op_index(getParam<unsigned int>("op_index"))
{
}

Real
ACSDDFMPoly::computeDFDOP(PFFunctionType type)
{
  Real SumEtaj = 0.0;
  for (unsigned int i = 0; i < _op_num; ++i)
    SumEtaj += (*_vals[i])[_qp] * (*_vals[i])[_qp];

  Real SumEtai2 = SumEtaj + _u[_qp] * _u[_qp];

  const auto & op_to_grain = _grain_tracker.getVarToFeatureVector(_current_elem->id());
  const auto grn_index = op_to_grain[_op_index];
  Real rho_i = _rho_eff[_qp];

  if (grn_index != FeatureFloodCount::invalid_id && grn_index < _GNDs_provider.getAvgDataSize())
    rho_i = _GNDs_provider.getAvgData(grn_index)._custom[0] * (_length_scale * _length_scale);

  // Calculate the contributions of the deformation energy to the residual and Jacobian
  Real drho_eff_detai = 2.0 * _u[_qp] * (rho_i - _rho_eff[_qp]) / SumEtai2;

  Real scaling_factor = 1.0;

  // Calculate the Stored Energy contribution to either the residual or Jacobian of the grain growth free energy
  switch (type)
  {
    case Residual:
      return scaling_factor * _beta[_qp] * drho_eff_detai;

    case Jacobian:
      return scaling_factor * _beta[_qp] * _phi[_j][_qp] *
             (2.0 * SumEtai2 * ((rho_i - _rho_eff[_qp]) - _u[_qp] * drho_eff_detai) -
              4.0 * _u[_qp] * _u[_qp] * (rho_i - _rho_eff[_qp])) /
             (SumEtai2 * SumEtai2);
  }
  mooseError("Invalid type passed in");
}
