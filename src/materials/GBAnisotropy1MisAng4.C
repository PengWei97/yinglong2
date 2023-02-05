#pragma once
#include "GBAnisotropy1MisAng4.h"

registerMooseObject("yinglongApp", GBAnisotropy1MisAng4);

InputParameters
GBAnisotropy1MisAng4::validParams()
{
  InputParameters params = GBAnisotropy1MisAngBase::validParams();
  params.addClassDescription(
      "Computes necessary material properties for the isotropic grain growth model");
  params.addParam<Real>("delta_theta_HAB", 15.0, "Benchmark value for GB anisotropy");
  params.addRequiredParam<std::vector<Real>>("misori_angle", 
                                             "the vector of the misorientation between grains for a piecewise function's independent varivable"); //  a piecewise function ~ x
  params.addRequiredParam<std::vector<Real>>("sigma_energy",
                                             "the vector of grain energy for a piecewise function's dependent varivable"); //  a piecewise function ~ y
  params.addParam<Real>("GBsigma_LAGB", 0.1, "GB energy of a low angle GB");
  params.addParam<Real>("GBsigma_HAGB", 1.0, "GB energy of a high angle GB");
  params.addParam<Real>("GBmob_LAGB", 2.5e-7, "GB mobility of a low angle GB");
  params.addParam<Real>("GBmob_HAGB", 2.5e-6, "GB mobility of a high angle GB");
  params.addParam<Real>("TT1_sigma", 0.1019, "initial value of twin boundary energy for {10-12} TT1 based on MD, J/m^2"); 
  params.addParam<Real>("CT1_sigma", 0.0616, "initial value of twin boundary energy for {11-22} CT1 based on MD, J/m^2");
  return params;  
}

GBAnisotropy1MisAng4::GBAnisotropy1MisAng4(const InputParameters & parameters)
  : GBAnisotropy1MisAngBase(parameters),
    _delta_theta_HAB(getParam<Real>("delta_theta_HAB")),
    _piecewise_func_sigma(getParam<std::vector<Real>>("misori_angle"),
                          getParam<std::vector<Real>>("sigma_energy")),
    _GBsigma_LAGB(getParam<Real>("GBsigma_LAGB")),
    _GBsigma_HAGB(getParam<Real>("GBsigma_HAGB")),
    _GBmob_LAGB(getParam<Real>("GBmob_LAGB")),
    _GBmob_HAGB(getParam<Real>("GBmob_HAGB")),
    _TT1_sigma(getParam<Real>("TT1_sigma")),
    _CT1_sigma(getParam<Real>("CT1_sigma"))
{
}

Real
GBAnisotropy1MisAng4::calculatedGBEnergy(const misoriAngle_isTwining & misori_gbType)
{
  Real gbSigma = _matrix_sigma; 

  if ((misori_gbType.misor < _delta_theta_HAB) && !_tb_anisotropy || !misori_gbType.isTwinning)
    gbSigma = _piecewise_func_sigma.sample(misori_gbType.misor);
  else if (_tb_anisotropy && misori_gbType.twinType == "twin_type0")
    gbSigma = _TT1_sigma;
  else if (_tb_anisotropy && misori_gbType.twinType == "twin_type1")
    gbSigma = _CT1_sigma;
  else if (misori_gbType.misor > _delta_theta_HAB)
    gbSigma = _GBsigma_LAGB;

  return gbSigma;
}

Real
GBAnisotropy1MisAng4::calculatedGBMobility(const misoriAngle_isTwining & misori_gbType)
{
  Real gbMob = _matrix_mob; 

  if ((misori_gbType.misor < _delta_theta_HAB) && !_tb_anisotropy || !misori_gbType.isTwinning)
    gbMob = _GBmob_HAGB * calculatedGBEnergy(misori_gbType) / _GBsigma_HAGB; // Eq.8
  else if (_tb_anisotropy && misori_gbType.twinType == "twin_type0")
    gbMob = _GBmob_HAGB * _TT1_sigma / _GBsigma_HAGB;
  else if (_tb_anisotropy && misori_gbType.twinType == "twin_type1")
    gbMob = _GBmob_HAGB * _CT1_sigma / _GBsigma_HAGB;
  else if (misori_gbType.misor > _delta_theta_HAB)
    gbMob = _GBmob_LAGB;

  return gbMob;
}