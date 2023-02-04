#pragma once
#include "GBAnisotropy1MisAng1.h"

registerMooseObject("yinglongApp", GBAnisotropy1MisAng1);

InputParameters
GBAnisotropy1MisAng1::validParams()
{
  InputParameters params = GBAnisotropy1MisAngBase::validParams();
  params.addClassDescription(
      "Computes necessary material properties for the isotropic grain growth model");
  params.addParam<Real>("delta_theta_HAB", 15.0, "Benchmark value for GB anisotropy");
  params.addParam<Real>("GBsigma_HAB", 0.9, "GB energy of a high angle GB");
  params.addParam<Real>("GBmob_HAB", 2.5e-6, "GB mobility of a high angle GB");
  return params;  
}

GBAnisotropy1MisAng1::GBAnisotropy1MisAng1(const InputParameters & parameters)
  : GBAnisotropy1MisAngBase(parameters),
    _delta_theta_HAB(getParam<Real>("delta_theta_HAB")),
    _GBsigma_HAB(getParam<Real>("GBsigma_HAB")),
    _GBmob_HAB(getParam<Real>("GBmob_HAB"))
{
}

Real
GBAnisotropy1MisAng1::calculatedGBEnergy(const misoriAngle_isTwining & misori_gbType)
{
  const Real B = 5;
  const Real n = 4;
  Real gbSigma = _matrix_sigma; 

  if ((misori_gbType.misor < _delta_theta_HAB) && !_tb_anisotropy || !misori_gbType.isTwinning)
    gbSigma = _GBsigma_HAB * (misori_gbType.misor / _delta_theta_HAB * (1 - std::log(misori_gbType.misor / _delta_theta_HAB) )); // Eq.7
  else
    gbSigma = _GBsigma_HAB;

  return gbSigma;
}

Real
GBAnisotropy1MisAng1::calculatedGBMobility(const misoriAngle_isTwining & misori_gbType)
{
   const Real B = 5;
   const Real n = 4;
   Real gbMob = _matrix_mob; 

   if ((misori_gbType.misor < _delta_theta_HAB) && !_tb_anisotropy || !misori_gbType.isTwinning)
      gbMob = _GBmob_HAB * ((1- std::exp(-B * std::pow( misori_gbType.misor / _delta_theta_HAB, n)))); // Eq.8
    else
      gbMob = _GBmob_HAB;

  return gbMob;
}