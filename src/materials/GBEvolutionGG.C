//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GBEvolutionGG.h"

registerMooseObject("yinglongApp", GBEvolutionGG);

InputParameters
GBEvolutionGG::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Computes necessary material properties for the isotropic grain growth model");
  params.addRequiredParam<UserObjectName>(
      "grain_tracker", "Name of GrainTrackerGG user object that provides Grain ID according to element ID");
  params.addRequiredParam<UserObjectName>("euler_angle_provider",
                                          "Name of Euler angle provider user object");
  params.addRequiredCoupledVar("T", "Temperature in Kelvin");
  params.addParam<Real>("f0s", 0.125, "The GB energy constant ");
  params.addRequiredParam<Real>("wGB", "Diffuse GB width in the length scale of the model");
  params.addRequiredParam<Real>("GBEnergy", "Grain boundary energy in J/m^2");
  params.addParam<Real>("length_scale", 1.0e-9, "Length scale in m, where default is nm");
  params.addParam<Real>("time_scale", 1.0e-9, "Time scale in s, where default is ns");
  params.addParam<Real>(
      "GBMobility",
      -1,
      "GB mobility input in m^4/(J*s), that overrides the temperature dependent calculation");
  params.addParam<Real>("GBmob0", 0, "Grain boundary mobility prefactor in m^4/(J*s)");
  params.addParam<Real>("Q", 0, "Grain boundary migration activation energy in eV");
  params.addParam<Real>("molar_volume",
                        24.62e-6,
                        "Molar volume in m^3/mol, needed for temperature gradient driving force");
  return params;
}

GBEvolutionGG::GBEvolutionGG(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _grain_tracker(getUserObject<GrainTrackerGG>("grain_tracker")),
    _euler(getUserObject<EulerAngleProvider>("euler_angle_provider")),
    _f0s(getParam<Real>("f0s")),
    _wGB(getParam<Real>("wGB")),
    _GBEnergy(getParam<Real>("GBEnergy")),
    _length_scale(getParam<Real>("length_scale")),
    _time_scale(getParam<Real>("time_scale")),
    _GBmob0(getParam<Real>("GBmob0")),
    _Q(getParam<Real>("Q")),
    _GBMobility(getParam<Real>("GBMobility")),
    _molar_vol(getParam<Real>("molar_volume")),
    _T(coupledValue("T")),
    _sigma(declareProperty<Real>("sigma")),
    _M_GB(declareProperty<Real>("M_GB")),
    _kappa(declareProperty<Real>("kappa_op")),
    _gamma(declareProperty<Real>("gamma_asymm")),
    _L(declareProperty<Real>("L")),
    _dLdT(parameters.hasDefaultCoupledValue("T")
              ? nullptr
              : &declarePropertyDerivative<Real>("L", getVar("T", 0)->name())),
    _l_GB(declareProperty<Real>("l_GB")),
    _mu(declareProperty<Real>("mu")),
    _entropy_diff(declareProperty<Real>("entropy_diff")),
    _molar_volume(declareProperty<Real>("molar_volume")),
    _act_wGB(declareProperty<Real>("act_wGB")),
    _misAngle(declareProperty<Real>("misAngle")),
    _twType(declareProperty<unsigned int>("twType")),
    _kb(8.617343e-5),     // Boltzmann constant in eV/K
    _JtoeV(6.24150974e18) // Joule to eV conversion
{
  if (_GBMobility == -1 && _GBmob0 == 0)
    mooseError("Either a value for GBMobility or for GBmob0 and Q must be provided");

  // Initialize _s_misoriTwin
  _s_misoriTwin.misor = -1.0;
  _s_misoriTwin.isTwinning = false;
  _s_misoriTwin.twinType = "none";
}

void
GBEvolutionGG::computeQpProperties()
{
  // eV/nm^2
  _sigma[_qp] = _GBEnergy * _JtoeV * (_length_scale * _length_scale);

  const Real length_scale4 = _length_scale * _length_scale * _length_scale * _length_scale;

  // GB mobility Derivative
  Real dM_GBdT;

  if (_GBMobility < 0)
  {
    // Convert to lengthscale^4/(eV*timescale);
    const Real M0 = _GBmob0 * _time_scale / (_JtoeV * length_scale4);

    _M_GB[_qp] = M0 * std::exp(-_Q / (_kb * _T[_qp]));
    dM_GBdT = MetaPhysicL::raw_value(_M_GB[_qp] * _Q / (_kb * _T[_qp] * _T[_qp]));
  }
  else
  {
    // Convert to lengthscale^4/(eV*timescale)
    _M_GB[_qp] = _GBMobility * _time_scale / (_JtoeV * length_scale4);
    dM_GBdT = 0.0;
  }

  // in the length scale of the system
  _l_GB[_qp] = _wGB;

  _L[_qp] = 4.0 / 3.0 * _M_GB[_qp] / _l_GB[_qp];
  if (_dLdT)
    (*_dLdT)[_qp] = MetaPhysicL::raw_value(4.0 / 3.0 * dM_GBdT / _l_GB[_qp]);
  _kappa[_qp] = 3.0 / 4.0 * _sigma[_qp] * _l_GB[_qp];
  _gamma[_qp] = 1.5;
  _mu[_qp] = 3.0 / 4.0 * 1.0 / _f0s * _sigma[_qp] / _l_GB[_qp];

  // J/(K mol) converted to eV(K mol)
  _entropy_diff[_qp] = 8.0e3 * _JtoeV;

  // m^3/mol converted to ls^3/mol
  _molar_volume[_qp] = _molar_vol / (_length_scale * _length_scale * _length_scale);
  _act_wGB[_qp] = 0.5e-9 / _length_scale;

  getMisorientationAngle();
}

void
GBEvolutionGG::getMisorientationAngle()
{
   // get the GB location based on the grainTrackerGG in the quadrature point
  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());

  std::vector<unsigned int> orderParameterIndex; // Create a vector of order parameter indices
  std::vector<unsigned int> grainIDIndex; // Create a vector of grain IDs  
  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index]; // grain id

    if (grain_id == FeatureFloodCount::invalid_id)
      continue;

    orderParameterIndex.push_back(op_index);
    grainIDIndex.push_back(grain_id);
  }

  // calculate the misorientation angle and GB anisotropy (GB energy and mobility)
  for (unsigned int i = 0; i < grainIDIndex.size() - 1; ++i)
    for (unsigned int j = i+1; j < grainIDIndex.size(); ++j)
    {
      auto angles_i = _euler.getEulerAngles(grainIDIndex[i]);
      auto angles_j = _euler.getEulerAngles(grainIDIndex[j]);    
      _s_misoriTwin = CalculateMisorientationAngle::calculateMisorientaion(angles_i, angles_j, _s_misoriTwin, "hcp");

      if (grainIDIndex.size() == 2)
      {
        _misAngle[_qp] = _s_misoriTwin.misor;
        if (_s_misoriTwin.isTwinning)
        {
          if (_s_misoriTwin.twinType == "twin_type0")
            _twType[_qp] = 0;
          else if (_s_misoriTwin.twinType == "twin_type1")
            _twType[_qp] = 1;
        }
      }
    }
}