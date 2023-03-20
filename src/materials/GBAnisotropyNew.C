#pragma once

#include "GBAnisotropyNew.h"
#include "MooseMesh.h"

#include <fstream>

registerMooseObject("yinglongApp", GBAnisotropyNew);

InputParameters
GBAnisotropyNew::validParams()
{
  InputParameters params = Material::validParams();
  params.addCoupledVar("T", 300.0, "Temperature in Kelvin");
  params.addParam<Real>("length_scale", 1.0e-9, "Length scale in m, where default is nm");
  params.addParam<Real>("time_scale", 1.0e-9, "Time scale in s, where default is ns");
  params.addRequiredParam<Real>("wGB", "the init diffuse GB width in nm");
  params.addParam<Real>(
      "delta_sigma", 0.1, "factor determining inclination dependence of GB energy");
  params.addParam<Real>(
      "delta_mob", 0.1, "factor determining inclination dependence of GB mobility");
  params.addRequiredParam<FileName>("anisotropic_GB_file_name",
                                    "Name of the file containing: 1)GB mobility prefactor; 2) GB "
                                    "migration activation energy; 3)GB energy");
  params.addRequiredParam<bool>("inclination_anisotropy",
                                "The GB anisotropy inclination would be considered if true");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  return params;
}

GBAnisotropyNew::GBAnisotropyNew(const InputParameters & parameters)
  : Material(parameters),
    _mesh_dimension(_mesh.dimension()),
    _length_scale(getParam<Real>("length_scale")),
    _time_scale(getParam<Real>("time_scale")),
    _delta_sigma(getParam<Real>("delta_sigma")),
    _delta_mob(getParam<Real>("delta_mob")),
    _wGB_init(getParam<Real>("wGB")),
    _anisotropic_GB_file_name(getParam<FileName>("anisotropic_GB_file_name")),
    _inclination_anisotropy(getParam<bool>("inclination_anisotropy")),
    _T(coupledValue("T")),
    _L(declareProperty<Real>("L")),
    _kappa(declareProperty<Real>("kappa_op")),
    _gamma_ij(declareProperty<std::vector<std::vector<Real>>>("gamma_asymm_ij")),
    _mu(declareProperty<Real>("mu")),
    _kb(8.617343e-5),      // Boltzmann constant in eV/K
    _JtoeV(6.24150974e18), // Joule to eV conversion
    _mu_qp(0.0),
    _kappa_qp(0.0),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _grad_vals(coupledGradients("v")),
    _act_wGB(declareProperty<Real>("act_wGB"))
{
  // reshape vectors
  _sigma.resize(_op_num);
  _mob.resize(_op_num);
  _Q.resize(_op_num);

  _L_kappa.resize(_op_num);
  _wGB_g2.resize(_op_num);

  for (unsigned int op = 0; op < _op_num; ++op)
  {
    _sigma[op].resize(_op_num);
    _mob[op].resize(_op_num);
    _Q[op].resize(_op_num);

    _L_kappa[op].resize(_op_num);
    _wGB_g2[op].resize(_op_num);
  }
}

void
GBAnisotropyNew::computeQpProperties()
{
  _gamma_ij[_qp].resize(_op_num);
  for (unsigned int op = 0; op < _op_num; ++op)
    _gamma_ij[_qp][op].assign(_op_num, 0.0);

  getGBAnisotropyFromFile();

  computerModelParameter();

  Real sum_L = 0.0; // for GB mobility anisotropy
  Real sum_wGB = 0.0; // for GB mobility anisotropy

  Real Val = 0.0;
  Real sum_val = 0.0;

  for (unsigned int m = 0; m < _op_num - 1; ++m)
  {
    for (unsigned int n = m + 1; n < _op_num; ++n) // m<n
    {
      Val = (100000.0 * ((*_vals[m])[_qp]) * ((*_vals[m])[_qp]) + 0.01) *
            (100000.0 * ((*_vals[n])[_qp]) * ((*_vals[n])[_qp]) + 0.01); // eta_i^2 * eta_j^2
      sum_val += Val;

      sum_L += _L_kappa[m][n] * Val;
      sum_wGB += _wGB_g2[m][n] * Val;
    }
  }

  _kappa[_qp] = _kappa_qp;
  _mu[_qp] = _mu_qp;
  _L[_qp] = sum_L / sum_val;
  _act_wGB[_qp] = sum_wGB/sum_val; 
}


void 
GBAnisotropyNew::getGBAnisotropyFromFile()
{
  // Read in data from "anisotropic_GB_file_name"
  std::ifstream inFile(_anisotropic_GB_file_name.c_str());

  if (!inFile)
    paramError("anisotropic_GB_file_name", "Can't open GB anisotropy input file");

  for (unsigned int i = 0; i < 2; ++i)
    inFile.ignore(255, '\n'); // ignore line

  Real data;
  for (unsigned int i = 0; i < 3 * _op_num; ++i)
  {
    std::vector<Real> row; // create an empty row of double values
    for (unsigned int j = 0; j < _op_num; ++j)
    {
      inFile >> data;
      row.push_back(data);
    }

    if (i < _op_num)
      _sigma[i] = row; // unit: J/m^2

    else if (i < 2 * _op_num)
      _mob[i - _op_num] = row; // unit: m^4/(J*s)

    else
      _Q[i - 2 * _op_num] = row; // unit: eV
  }

  inFile.close();
}

void
GBAnisotropyNew::computerModelParameter()
{
  Real sigma_init;
  Real g2 = 0.0;
  Real f_interf = 0.0;
  Real y = 0.0;
  Real yyy = 0.0;

  Real sigma_big = 0.0;
  Real sigma_small = 0.0;

  for (unsigned int m = 0; m < _op_num - 1; ++m)
    for (unsigned int n = m + 1; n < _op_num; ++n)
    {
      // Convert units of mobility and energy
      _sigma[m][n] *= _JtoeV * (_length_scale * _length_scale); // eV/nm^2
      _mob[m][n] *= _time_scale / (_JtoeV * (_length_scale * _length_scale * _length_scale * _length_scale)) * std::exp(-_Q[m][n] / (_kb * _T[_qp])); // Convert to nm^4/(eV*ns);

      if (m == 0 && n == 1)
      {
        sigma_big = _sigma[m][n];
        sigma_small = sigma_big;
      }

      else if (_sigma[m][n] > sigma_big)
        sigma_big = _sigma[m][n];

      else if (_sigma[m][n] < sigma_small)
        sigma_small = _sigma[m][n];
    }


  sigma_init = (sigma_big + sigma_small) / 2.0; // how to get the initial sigma and GB width ??

  _mu_qp = 6.0 * sigma_init / _wGB_init;

  _kappa_qp = 3.0/4.0 * sigma_init * _wGB_init;
  for (unsigned int m = 0; m < _op_num - 1; ++m)
  {
    for (unsigned int n = m + 1; n < _op_num; ++n)
    {
      g2 = _sigma[m][n] * _sigma[m][n] / (_kappa_qp * _mu_qp);

      y = -3.0944 * g2 * g2 * g2 * g2 - 1.81694 * g2 * g2 * g2 + 10.323 * g2 * g2 - 8.1819 * g2 + 2.0033; // MD2022

      _gamma_ij[_qp][m][n] = 1 / y; 
      _gamma_ij[_qp][n][m] = 1 / y;

      _L_kappa[m][n] = _sigma[m][n] * _mob[m][n] / _kappa_qp; 
      yyy = y * y * y;

      f_interf = 0.078815 * yyy * yyy - 0.49546 * yyy * y * y + 1.2244 * yyy * y - 1.5281 * yyy + 1.0686 * y * y - 0.55631 * y + 0.29067; // MD2022
      
      _wGB_g2[m][n] = std::sqrt(_kappa_qp/(_mu_qp * f_interf));
      _wGB_g2[n][m] = g2;
    }  
  }
}
