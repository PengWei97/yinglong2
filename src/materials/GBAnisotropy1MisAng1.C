#pragma once
#include "GBAnisotropy1MisAng1.h"

registerMooseObject("yinglongApp", GBAnisotropy1MisAng1);

InputParameters
GBAnisotropy1MisAng1::validParams()
{
  InputParameters params = GBAnisotropy1MisAngBase::validParams();
  params.addClassDescription(
      "Computes necessary material properties for the anisotropic grain growth model");
  params.addParam<Real>("delta_theta_HAGB", 15.0, "Benchmark value for GB anisotropy");
  params.addParam<Real>("GBsigma_HAGB", 0.9, "GB energy of a high angle GB");
  params.addParam<Real>("GBmob_HAGB", 2.5e-6, "GB mobility of a high angle GB");
  return params;  
}

GBAnisotropy1MisAng1::GBAnisotropy1MisAng1(const InputParameters & parameters)
  : GBAnisotropy1MisAngBase(parameters),
    _delta_theta_HAGB(getParam<Real>("delta_theta_HAGB")),
    _GBsigma_HAGB(getParam<Real>("GBsigma_HAGB")),
    _GBmob_HAGB(getParam<Real>("GBmob_HAGB"))
{
}

Real
GBAnisotropy1MisAng1::calculatedGBEnergy(const MisorientationAngleData & s_misorientation_angle)
{
  Real const & delta_theta = s_misorientation_angle._misor;

  Real gbSigma = _GBsigma_HAGB * ((delta_theta / _delta_theta_HAGB * (1 - std::log(delta_theta / _delta_theta_HAGB)))); // Eq.8

  return gbSigma;
}

Real
GBAnisotropy1MisAng1::calculatedGBMobility(const MisorientationAngleData & s_misorientation_angle)
{
  // Equation constant
  Real B = 5;
  Real n = 4;

  // GB mobility based on misorientation angle
  Real gbMob;

  Real const & delta_theta = s_misorientation_angle._misor;
  if (delta_theta >  _delta_theta_HAGB)
    gbMob = _GBmob_HAGB;
  else
    gbMob = _GBmob_HAGB * ((1- std::exp(-B * std::pow( delta_theta / _delta_theta_HAGB, n)))); // Eq.8

  return gbMob;
}