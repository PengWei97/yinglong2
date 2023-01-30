//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ACGrGrPolyMOP2.h"

registerMooseObject("yinglongApp", ACGrGrPolyMOP2);

InputParameters
ACGrGrPolyMOP2::validParams()
{
  InputParameters params = ACGrGrBase::validParams();
  params.addClassDescription("Grain-Boundary model poly-crystalline interface Allen-Cahn Kernel for the new MOP-PF model in Nele Moelans-MD-2022");
  params.addRequiredParam<unsigned int>("cur_op","the current id of the order parameter");
  return params;
}

ACGrGrPolyMOP2::ACGrGrPolyMOP2(const InputParameters & parameters)
  : ACGrGrBase(parameters), 
  _gamma_ij(getMaterialProperty<std::vector<std::vector<Real>>>("gamma_asymm_ij")),
  _cur_op(getParam<unsigned int>("cur_op"))
{
}

Real
ACGrGrPolyMOP2::computeDFDOP(PFFunctionType type)
{

  // Sum all other order parameters
  Real SumEtaj_gamma = 0.0;

  for (unsigned int i = 0; i < _op_num; ++i)
    SumEtaj_gamma += (*_vals[i])[_qp] * (*_vals[i])[_qp] * _gamma_ij[_qp][_cur_op][_vals_var[i]];

  // Calculate either the residual or Jacobian of the grain growth free energy
  switch (type)
  {
    case Residual:
    {
      return _mu[_qp] *
             (_u[_qp] * _u[_qp] * _u[_qp] - _u[_qp] + 2.0 * _u[_qp] * SumEtaj_gamma);
    }

    case Jacobian:
    {
      return _mu[_qp] *
             (_phi[_j][_qp] * (3.0 * _u[_qp] * _u[_qp] - 1.0 + 2.0 * SumEtaj_gamma));
    }

    default:
      mooseError("Invalid type passed in");
  }
}

Real
ACGrGrPolyMOP2::computeQpOffDiagJacobian(unsigned int jvar)
{
  for (unsigned int i = 0; i < _op_num; ++i)
    if (jvar == _vals_var[i])
    {
      // Derivative of SumEtaj
      const Real dSumEtaj_gamma = 2.0 * (*_vals[i])[_qp] * _phi[_j][_qp] * _gamma_ij[_qp][_cur_op][_vals_var[i]];

      const Real dDFDOP = _mu[_qp] * 2.0 *  _u[_qp] * dSumEtaj_gamma;

      return _L[_qp] * _test[_i][_qp] * dDFDOP;
    }

  return 0.0;
}