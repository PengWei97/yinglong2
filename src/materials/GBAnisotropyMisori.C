//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GBAnisotropyMisori.h"
#include "MooseMesh.h"
#include <cmath>

registerMooseObject("yinglongApp", GBAnisotropyMisori);

InputParameters
GBAnisotropyMisori::validParams()
{
  InputParameters params = GBAnisotropyMisoriBase::validParams();
  params.addClassDescription(
      "Computes necessary material properties for the anisotropic grain growth model");
  params.addParam<Real>("TT1_sigma",  0.9, "Twin boundary energy for {10-12} tensile twin (type 1) based on MD, J/m^2");
  params.addParam<Real>("CT1_sigma",  0.9, "Twin boundary energy for {11-22} compresssion twin (type 1) based on MD, J/m^2");  
  params.addParam<Real>("TT1_mob", 2.5e-6, "Twin boundary mobility for {10-12} tensile twin (type 1) based on experiment, m^4/(J*s)");
  params.addParam<Real>("CT1_mob", 2.5e-6, "Twin boundary mobility for {11-22} compresssion twin (type 1) based on experiment, m^4/(J*s)");
  params.addRequiredParam<UserObjectName>(
      "grain_tracker", "Name of GrainTracker user object that provides Grain ID according to element ID");
  params.addRequiredParam<UserObjectName>("euler_angle_provider",
                                          "Name of Euler angle provider user object");
  params.addParam<bool>("gb_energy_anisotropy", false,
      "The GB energy anisotropy based on misorientation would be considered if true");
  params.addParam<bool>("gb_mobility_anisotropy", true,
      "The GB mobility anisotropy would be considered if true");                                          
  return params;
}

GBAnisotropyMisori::GBAnisotropyMisori(const InputParameters & parameters)
  : GBAnisotropyMisoriBase(parameters),
    _TT1_sigma(getParam<Real>("TT1_sigma")),
    _CT1_sigma(getParam<Real>("CT1_sigma")),
    _TT1_mob(getParam<Real>("TT1_mob")),
    _CT1_mob(getParam<Real>("CT1_mob")),
    _grain_tracker(getUserObject<GrainTracker>("grain_tracker")),
    _euler(getUserObject<EulerAngleProvider>("euler_angle_provider")),
    _gb_energy_anisotropy(getParam<bool>("gb_energy_anisotropy")),
    _gb_mobility_anisotropy(getParam<bool>("gb_mobility_anisotropy")),
    _misori_angle(declareProperty<Real>("misori_angle")),
    _twinning_type(declareProperty<Real>("twinning_type"))  
{
}

void
GBAnisotropyMisori::computeGBProperties()
{
  _misori_angle[_qp] = 0.0;
  _twinning_type[_qp] = -1.0;

  // get the GB location based on the GrainTracker in the quadrature point
  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id()); 

  std::vector<unsigned int> var_index;
  std::vector<unsigned int> grain_id_index;
  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index]; // grain id

    if (grain_id == FeatureFloodCount::invalid_id)
      continue;

    var_index.push_back(op_index);
    grain_id_index.push_back(grain_id);
  }

  // Set other unset mob_ij or sigma_ij to this value
  Real mob_min =  _GBmob_HAGB * 0.1;
  Real sigma_min =  _GBsigma_HAGB * 0.1;
   
  // When at grain boundaries or junction
  if (grain_id_index.size() > 1)
  {
    // Traverse to obtain the sigma_ij, mob_ij of the activation order parameters at the orthogonal point
    for (unsigned int i = 0; i < grain_id_index.size() - 1; ++i)
      for (unsigned int j = i+1; j < grain_id_index.size(); ++j)
      {
        // calculate misorientation angle
        auto angles_i = _euler.getEulerAngles(grain_id_index[i]);
        auto angles_j = _euler.getEulerAngles(grain_id_index[j]);
        _misori_s = MisorientationAngleCalculator::calculateMisorientaion(angles_i, angles_j, _misori_s);

        if (_gb_energy_anisotropy)
          _sigma[var_index[i]][var_index[j]] = calculatedGBEnergy(_misori_s);
        else
          _sigma[var_index[i]][var_index[j]] = _GBsigma_HAGB;

        if (_gb_mobility_anisotropy)
          _mob[var_index[i]][var_index[j]] = calculatedGBMobility(_misori_s);
        else
          _mob[var_index[i]][var_index[j]] = _GBmob_HAGB;

        _sigma[var_index[j]][var_index[i]] =  _sigma[var_index[i]][var_index[j]];
        _mob[var_index[j]][var_index[i]] =  _mob[var_index[i]][var_index[j]];

        // Set the twinning type
        if (i == 0 && j == 1)
        {
          _misori_angle[_qp] =  _misori_s._misor;

          if (_misori_s._is_twin && (_misori_s._twin_type == TwinType::TT1))
            _twinning_type[_qp] = 1.0;
          else if (_misori_s._is_twin && (_misori_s._twin_type == TwinType::CT1))
            _twinning_type[_qp] = 2.0;
          else
            _twinning_type[_qp] = 0.0; // GB

          mob_min = _mob[var_index[i]][var_index[j]];
          sigma_min = _sigma[var_index[i]][var_index[j]];
        }

        if (mob_min > _mob[var_index[i]][var_index[j]])
          mob_min = _mob[var_index[i]][var_index[j]];

        if (sigma_min > _sigma[var_index[i]][var_index[j]])
          sigma_min = _sigma[var_index[i]][var_index[j]];
      }

    // Traverse the entire _mob_ij and _sigam_ij and set other values to mob_min or sigma_min
    for (unsigned int i = 0; i < _mob.size(); ++i)
      for (unsigned int j = 0; j < _mob.size(); ++j)
      {
        if (_mob[i][j] < mob_min)
        {
          _mob[i][j] = mob_min;
          _mob[j][i] = mob_min;
        }

        if (_sigma[i][j] < sigma_min)
        {
          _sigma[i][j] = sigma_min;
          _sigma[j][i] = sigma_min;
        }
      }
    
  }
}

Real
GBAnisotropyMisori::calculatedGBEnergy(const MisorientationAngleData & misori_s)
{
  Real gbSigma = _GBsigma_HAGB;

  // transition misorientation angle between low and high-angle grain boundary
  Real trans_misori_angle_HAGB = 15.0;

  if ((_misori_s._misor <= 2.0) && !misori_s._is_twin)
    gbSigma = _GBsigma_HAGB * ((2.0 / trans_misori_angle_HAGB * (1 - std::log(2.0 / trans_misori_angle_HAGB))));    
  else if ((_misori_s._misor <= trans_misori_angle_HAGB) && !misori_s._is_twin)
    gbSigma = _GBsigma_HAGB * ((misori_s._misor / trans_misori_angle_HAGB * (1 - std::log(misori_s._misor / trans_misori_angle_HAGB))));
  else if (misori_s._is_twin)
  {
    if (misori_s._twin_type == TwinType::TT1)
      gbSigma = _TT1_sigma;
    else if (misori_s._twin_type == TwinType::CT1)
      gbSigma = _CT1_sigma;
  }

  return gbSigma; 
}

Real
GBAnisotropyMisori::calculatedGBMobility(const MisorientationAngleData & misori_s)
{
  // Initialize GB mobility
  Real gbMob = _GBmob_HAGB;

  // transition misorientation angle between low and high-angle grain boundary
  Real trans_misori_angle_HAGB = 15;

  // Equation constant
  Real B = 5;
  Real n = 4;
  
  if ((_misori_s._misor <= 2.0) && !misori_s._is_twin)
    gbMob = _GBmob_HAGB * ((1- std::exp(-B * std::pow( 2.0 / trans_misori_angle_HAGB, n)))); // Eq.8    
  else if ((_misori_s._misor <=  trans_misori_angle_HAGB) && !misori_s._is_twin )
    gbMob = _GBmob_HAGB * ((1- std::exp(-B * std::pow( misori_s._misor / trans_misori_angle_HAGB, n)))); // Eq.8
  else if (misori_s._is_twin)
  {
    if (misori_s._twin_type == TwinType::TT1)
      gbMob = _TT1_mob;
    else if (misori_s._twin_type == TwinType::CT1)
      gbMob = _CT1_mob;
  }

  return gbMob;
}