#pragma once
#include "GBAnisotropyMisAng2.h"

registerMooseObject("yinglongApp", GBAnisotropyMisAng2);

InputParameters
GBAnisotropyMisAng2::validParams()
{
  InputParameters params = GBAnisotropyMisAngBase::validParams();
  params.addClassDescription(
      "Computes necessary material properties for the anisotropic grain growth model");
  params.addParam<Real>("delta_theta_HAGB", 15.0, "Benchmark value for GB anisotropy");
  params.addParam<Real>("GBsigma_HAGB", 0.9, "GB energy of a high angle GB");
  params.addParam<Real>("GBmob_HAGB", 2.5e-6, "GB mobility of a high angle GB");
  params.addParam<Real>("TT1_sigma", 2.5e-6, "Twin boundary energy for {10-12} tensile twin (type 1) based on MD, J/m^2");
  params.addParam<Real>("CT1_sigma", 2.5e-6, "Twin boundary energy for {11-22} tensile twin (type 1) based on MD, J/m^2");
  return params;  
}

GBAnisotropyMisAng2::GBAnisotropyMisAng2(const InputParameters & parameters)
  : GBAnisotropyMisAngBase(parameters),
    _delta_theta_HAGB(getParam<Real>("delta_theta_HAGB")),
    _GBsigma_HAGB(getParam<Real>("GBsigma_HAGB")),
    _GBmob_HAGB(getParam<Real>("GBmob_HAGB")),
    _TT1_sigma(getParam<Real>("TT1_sigma")),
    _CT1_sigma(getParam<Real>("CT1_sigma"))
{
}

Real
GBAnisotropyMisAng2::calculatedGBEnergy(const MisorientationAngleData & s_misorientation_angle)
{
  Real const & delta_theta = s_misorientation_angle._misor;
  bool const & is_twin = s_misorientation_angle._is_twin;

  // 初始化为高角度晶界的晶界能
  Real gbSigma = _GBsigma_HAGB;

  if (delta_theta <= _delta_theta_HAGB && !is_twin) // 不是孪晶且取向差为小角度晶界时
    gbSigma = _GBsigma_HAGB * ((delta_theta / _delta_theta_HAGB * (1 - std::log(delta_theta / _delta_theta_HAGB)))); // Eq.8
  else if (is_twin)
  {
    if (s_misorientation_angle._twin_type == TwinType::TT1) // 当判定为孪晶界，且为TT1时
      gbSigma = _TT1_sigma;
    else if (s_misorientation_angle._twin_type == TwinType::TT1) // 当判定为孪晶界，且为CT1时
      gbSigma = _CT1_sigma;
  }

  return gbSigma;
}

Real
GBAnisotropyMisAng2::calculatedGBMobility(const MisorientationAngleData & s_misorientation_angle)
{
  // Equation constant
  Real B = 5;
  Real n = 4;

  // 初始化为高角度晶界的晶界能
  Real gbMob = _GBmob_HAGB;

  Real const & delta_theta = s_misorientation_angle._misor;
  bool const & is_twin = s_misorientation_angle._is_twin;

  // 暂时没有考虑孪晶的低能特性
  if (delta_theta <=  _delta_theta_HAGB)
    gbMob = _GBmob_HAGB * ((1- std::exp(-B * std::pow( delta_theta / _delta_theta_HAGB, n)))); // Eq.8

  return gbMob;
}