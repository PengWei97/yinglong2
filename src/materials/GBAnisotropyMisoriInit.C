//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GBAnisotropyMisoriInit.h"
#include "MooseMesh.h"
#include <cmath>

registerMooseObject("yinglongApp", GBAnisotropyMisoriInit);

InputParameters
GBAnisotropyMisoriInit::validParams()
{
  InputParameters params = GBAnisotropyMisoriBase::validParams();
  params.addClassDescription(
      "Computes necessary material properties for the anisotropic grain growth model");
  params.addParam<Real>("rate1_HABvsLAB_sigma", 0.5, "Pre-coefficient in the Read-Shockley law");
  params.addParam<Real>("rate2_HABvsLAB_sigma", 0.5, "Post-coefficient in the Read-Shockley law");
  params.addParam<Real>("rate1_HABvsLAB_mob", 0.5, "Pre-coefficient in the sigmoidal law");
  params.addParam<Real>("rate2_HABvsLAB_mob", 0.5, "Post-coefficient in the sigmoidal law");
  params.addRequiredParam<UserObjectName>(
      "grain_tracker", "Name of GrainTracker user object that provides Grain ID according to element ID");
  params.addRequiredParam<UserObjectName>("euler_angle_provider",
                                          "Name of Euler angle provider user object");
  params.addParam<bool>("gb_energy_anisotropy", false,
      "The GB energy anisotropy based on misorientation would be considered if true");
  params.addParam<bool>("gb_mobility_anisotropy", false,
      "The GB mobility anisotropy would be considered if true");                                          
  return params;
}

GBAnisotropyMisoriInit::GBAnisotropyMisoriInit(const InputParameters & parameters)
  : GBAnisotropyMisoriBase(parameters),
    _rate1_HABvsLAB_sigma(getParam<Real>("rate1_HABvsLAB_sigma")),
    _rate2_HABvsLAB_sigma(getParam<Real>("rate2_HABvsLAB_sigma")),
    _rate1_HABvsLAB_mob(getParam<Real>("rate1_HABvsLAB_mob")),
    _rate2_HABvsLAB_mob(getParam<Real>("rate2_HABvsLAB_mob")),
    _grain_tracker(getUserObject<GrainTracker>("grain_tracker")),
    _euler(getUserObject<EulerAngleProvider>("euler_angle_provider")),
    _gb_energy_anisotropy(getParam<bool>("gb_energy_anisotropy")),
    _gb_mobility_anisotropy(getParam<bool>("gb_mobility_anisotropy")),
    _misori_angle(declareProperty<Real>("misori_angle"))
{
}

void
GBAnisotropyMisoriInit::computeGBProperties()
{
  _misori_angle[_qp] = 0.0;

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
   
  Real sigma_min = _GBsigma_HAGB;
  Real sigma_max = _GBsigma_HAGB;
  Real mob_min = _GBmob_HAGB;
  Real mob_max = _GBmob_HAGB;

  // When at grain boundaries or junction
  if (grain_id_index.size() > 1)
  {
    std::fill(_sigma.begin(), _sigma.end(), std::vector<Real>(_op_num, 0.0));
    std::fill(_mob.begin(), _mob.end(), std::vector<Real>(_op_num, 0.0));

    // Traverse to obtain the sigma_ij, mob_ij of the activation order parameters at the orthogonal point
    for (unsigned int i = 0; i < grain_id_index.size() - 1; ++i)
    {
      auto angles_i = _euler.getEulerAngles(grain_id_index[i]);
      for (unsigned int j = i+1; j < grain_id_index.size(); ++j)
      {
        auto angles_j = _euler.getEulerAngles(grain_id_index[j]);

        auto & _sigma_ij = _sigma[var_index[i]][var_index[j]];
        auto & _mob_ij = _mob[var_index[i]][var_index[j]];

        // calculate misorientation angle
        _misori_s = MisorientationAngleCalculator::calculateMisorientaion(angles_i, angles_j, _misori_s);

        if (_gb_energy_anisotropy)
          _sigma_ij = calculatedGBEnergy(_misori_s);
        else
          _sigma_ij = _GBsigma_HAGB;

        if (_gb_mobility_anisotropy)
          _mob_ij = calculatedGBMobility(_misori_s);
        else
          _mob_ij = _GBmob_HAGB;

        _sigma[var_index[j]][var_index[i]] =  _sigma_ij;
        _mob[var_index[j]][var_index[i]] =  _mob_ij;

        // fine the maximum and minimum
        if (sigma_min > _sigma_ij)
          sigma_min = _sigma_ij;
        else if (sigma_max < _sigma_ij)
          sigma_max = _sigma_ij;

        if (mob_min > _mob_ij)
          mob_min = _mob_ij;
        else if (mob_max < _mob_ij)
          mob_max = _mob_ij;
      }
    }

    for (unsigned int i = 0; i < _op_num; ++i)
      for (unsigned int j = 0; j < _op_num; ++j)
      {
        if ((_sigma[i][j] == 0.0))
        {
          _sigma[i][j] = (sigma_max + sigma_min) / 2.0;
          _sigma[j][i] = (sigma_max + sigma_min) / 2.0;
        }

        if ((_mob[i][j] == 0.0))
        {
          _mob[i][j] = (mob_max + mob_min) / 2.0;
          _mob[j][i] = (mob_max + mob_min) / 2.0;
        }
      }
  }
}

Real
GBAnisotropyMisoriInit::calculatedGBEnergy(const MisorientationAngleData & misori_s)
{
  Real gbSigma = _GBsigma_HAGB;

  // transition misorientation angle between low and high-angle grain boundary
  Real trans_misori_HAGB = 15.0;

  if (_misori_s._misor <= trans_misori_HAGB)
    gbSigma = _GBsigma_HAGB * (_rate1_HABvsLAB_sigma * (misori_s._misor / trans_misori_HAGB * (1 - std::log(misori_s._misor / trans_misori_HAGB))) + _rate2_HABvsLAB_sigma);
  else
    gbSigma = _GBsigma_HAGB * (_rate1_HABvsLAB_sigma * (trans_misori_HAGB / trans_misori_HAGB * (1 - std::log(trans_misori_HAGB / trans_misori_HAGB))) + _rate2_HABvsLAB_sigma);

  return gbSigma; 
}

Real
GBAnisotropyMisoriInit::calculatedGBMobility(const MisorientationAngleData & misori_s)
{
  // Initialize GB mobility
  Real gbMob = _GBmob_HAGB;

  // transition misorientation angle between low and high-angle grain boundary
  Real trans_misori_HAGB = 15;

  // Equation constant
  Real B = 5;
  Real n = 4;
  
  if (_misori_s._misor <= trans_misori_HAGB)
    gbMob = _GBmob_HAGB * (_rate1_HABvsLAB_mob * (1- std::exp(-B * std::pow( misori_s._misor / trans_misori_HAGB, n))) + _rate2_HABvsLAB_mob); // Eq.8
  else
    gbMob =  _GBmob_HAGB * (_rate1_HABvsLAB_mob * (1- std::exp(-B * std::pow(trans_misori_HAGB / trans_misori_HAGB, n))) + _rate2_HABvsLAB_mob); // Eq.8

  return gbMob;
}