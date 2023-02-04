#pragma once
#include "GBAnisotropy1MisAng3.h"

registerMooseObject("yinglongApp", GBAnisotropy1MisAng3);

InputParameters
GBAnisotropy1MisAng3::validParams()
{
  InputParameters params = GBAnisotropy1MisAngBase::validParams();
  params.addClassDescription(
      "Computes necessary material properties for the isotropic grain growth model");
  params.addParam<Real>("delta_theta_HAB", 15.0, "Benchmark value for GB anisotropy");
  params.addParam<Real>("GBsigma_HAB", 0.9, "GB energy of a high angle GB");
  params.addParam<Real>("GBmob_HAB", 2.5e-6, "GB mobility of a high angle GB");
  params.addParam<Real>("TT1_sigma", 0.1019, "initial value of twin boundary energy for {10-12} TT1 based on MD, J/m^2"); 
  params.addParam<Real>("CT1_sigma", 0.0616, "initial value of twin boundary energy for {11-22} CT1 based on MD, J/m^2");
  return params;  
}

GBAnisotropy1MisAng3::GBAnisotropy1MisAng3(const InputParameters & parameters)
  : GBAnisotropy1MisAngBase(parameters),
    _delta_theta_HAB(getParam<Real>("delta_theta_HAB")),
    _GBsigma_HAB(getParam<Real>("GBsigma_HAB")),
    _GBmob_HAB(getParam<Real>("GBmob_HAB")),
    _TT1_sigma(getParam<Real>("TT1_sigma")),
    _CT1_sigma(getParam<Real>("CT1_sigma"))
{
}

Real
GBAnisotropy1MisAng3::calculatedGBEnergy(const misoriAngle_isTwining & misori_gbType)
{
  const Real B = 5;
  const Real n = 4;
  Real gbSigma = _matrix_sigma; 

  if ((misori_gbType.misor < _delta_theta_HAB) && !_tb_anisotropy || !misori_gbType.isTwinning)
    gbSigma = _GBsigma_HAB * ((misori_gbType.misor / _delta_theta_HAB * (1 - std::log(misori_gbType.misor / _delta_theta_HAB)))); // Eq.7
  else if (_tb_anisotropy && misori_gbType.twinType == "twin_type0")
    gbSigma = _TT1_sigma;
  else if (_tb_anisotropy && misori_gbType.twinType == "twin_type1")
    gbSigma = _CT1_sigma;
  else
    gbSigma = _GBsigma_HAB;

  return gbSigma;
}

Real
GBAnisotropy1MisAng3::calculatedGBMobility(const misoriAngle_isTwining & misori_gbType)
{
  const Real B = 5;
  const Real n = 4;
  Real gbMob = _matrix_mob; 

  if (!_tb_anisotropy || !misori_gbType.isTwinning)
    gbMob = _GBmob_HAB * ((1- std::exp(-B * std::pow( misori_gbType.misor / _delta_theta_HAB, n)))); // Eq.8
  else if (_tb_anisotropy && misori_gbType.twinType == "twin_type0")
    gbMob = _GBmob_HAB * _TT1_sigma / _GBsigma_HAB;
  else if (_tb_anisotropy && misori_gbType.twinType == "twin_type1")
    gbMob = _GBmob_HAB * _CT1_sigma / _GBsigma_HAB;
  else 
    gbMob = _GBmob_HAB;

  return gbMob;
}