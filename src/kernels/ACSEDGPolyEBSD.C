//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ACSEDGPolyEBSD.h"

registerMooseObject("yinglongApp", ACSEDGPolyEBSD);

InputParameters
ACSEDGPolyEBSD::validParams()
{
  InputParameters params = ACBulk<Real>::validParams();
  params.addClassDescription("Stored Energy contribution to grain growth");
  params.addRequiredCoupledVar("v", "Array of coupled variable names");
  params.addParam<Real>("length_scale", 1.0e-9, "Length scale in m, where default is nm");
  params.addParam<bool>("concurrent_recovery", false, "The concurrent recovery would be considered if true");
  params.addRequiredParam<UserObjectName>("grain_tracker",
                                          "The GrainTracker UserObject to get values from.");
  params.addRequiredParam<UserObjectName>("GNDs_provider", "GNDs provider for EBSD reader");                                    
  params.addRequiredParam<unsigned int>("op_index", "The index for the current order parameter");
  return params;
}

ACSEDGPolyEBSD::ACSEDGPolyEBSD(const InputParameters & parameters)
  : ACBulk<Real>(parameters),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _vals_var(coupledIndices("v")),
    _length_scale(getParam<Real>("length_scale")),
    _concurrent_recovery(getParam<bool>("concurrent_recovery")),  
    _beta(getMaterialProperty<Real>("beta")),
    _rho_eff(getMaterialProperty<Real>("rho_eff")),
    _grain_tracker(getUserObject<GrainTrackerInterface>("grain_tracker")),
    _GNDs_provider(getUserObject<EBSDReader>("GNDs_provider")),
    _op_index(getParam<unsigned int>("op_index"))
{
}

Real
ACSEDGPolyEBSD::computeDFDOP(PFFunctionType type)
{
  Real SumEtaj = 0.0;
  for (unsigned int i = 0; i < _op_num; ++i)
    SumEtaj += (*_vals[i])[_qp] * (*_vals[i])[_qp];

  // Add the current OP to the sum
  Real SumEtai2 = SumEtaj + _u[_qp] * _u[_qp];

  const auto & op_to_grain = _grain_tracker.getVarToFeatureVector(_current_elem->id());
  const auto grn_index = op_to_grain[_op_index]; // Grain ID
  Real rho_i = getGNDsFromEBSD(grn_index);

  // Calculate the contributions of the deformation energy to the residual and Jacobian
  Real drho_eff_detai = 2.0 * _u[_qp] * (rho_i - _rho_eff[_qp]) / SumEtai2;

  // Calculate the Stored Energy contribution to either the residual or Jacobian of the grain growth
  // free energy
  switch (type)
  {
    case Residual:
      return _beta[_qp] * drho_eff_detai;

    case Jacobian:
      return _beta[_qp] * _phi[_j][_qp] *
             (2.0 * SumEtai2 * ((rho_i - _rho_eff[_qp]) - _u[_qp] * drho_eff_detai) -
              4.0 * _u[_qp] * _u[_qp] * (rho_i - _rho_eff[_qp])) /
             (SumEtai2 * SumEtai2);
  }
  mooseError("Invalid type passed in");
}

Real
ACSEDGPolyEBSD::getGNDsFromEBSD(const unsigned int & grain_id)
{
  Real rho_i = _rho_eff[_qp];

  if (grain_id != FeatureFloodCount::invalid_id && grain_id < _GNDs_provider.getGrainNum())
    rho_i = _GNDs_provider.getAvgData(grain_id)._custom[0] * (_length_scale * _length_scale);

  // concurrent recovery
  const Real rho_end = 3.0e7 * (_length_scale * _length_scale); // the minimum GND due to recovery 1/m^2
  auto & time_current = _fe_problem.time(); // current simulation time s

  if (_concurrent_recovery && rho_i > rho_end)
    rho_i = (rho_i - rho_end) * std::exp(-8.535e-4 * time_current) + rho_end; // fitting function considering concurrent recovery

  return rho_i;
}
