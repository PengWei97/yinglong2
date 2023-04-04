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
  params.addParam<Real>("GBsigma_HAGB", 0.9, "GB energy of a high angle GB");
  params.addParam<Real>("GBmob_HAGB", 2.5e-6, "GB mobility of a high angle GB");
  params.addRequiredParam<UserObjectName>(
      "grain_tracker", "Name of GrainTracker user object that provides Grain ID according to element ID");
  params.addRequiredParam<UserObjectName>("euler_angle_provider",
                                          "Name of Euler angle provider user object");
  params.addParam<bool>("gb_energy_anisotropy", true,
      "The GB energy anisotropy based on misorientation would be considered if true");
  params.addParam<bool>("gb_mobility_anisotropy", false,
      "The GB mobility anisotropy would be considered if true");                                          
  return params;
}

GBAnisotropyMisori::GBAnisotropyMisori(const InputParameters & parameters)
  : GBAnisotropyMisoriBase(parameters),
    _GBsigma_HAGB(getParam<Real>("GBsigma_HAGB")),
    _GBmob_HAGB(getParam<Real>("GBmob_HAGB")),
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

  std::vector<unsigned int> var_index;  // Create a vector of order parameter indices
  std::vector<unsigned int> grain_id_index; // Create a vector of grain IDs   

  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index]; // grain id

    if (grain_id == FeatureFloodCount::invalid_id)
      continue;

    var_index.push_back(op_index);
    grain_id_index.push_back(grain_id);
  }

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

        // Calculate sigma_ij
        if (_gb_energy_anisotropy)
        {
          _sigma[var_index[i]][var_index[j]] = calculatedGBEnergy(_misori_s);
          _sigma[var_index[j]][var_index[i]] =  _sigma[var_index[i]][var_index[j]];
        }

        // calculate mob_ij
        if (_gb_mobility_anisotropy)
        {
          // _mob[var_index[i]][var_index[j]] = calculatedGBMobility(_misori_s);

          // TODO: The following code is only used to test the applicability of GBAnisotropyMisori and will be deleted later.
          // if (grain_id_index[i] == 1 || grain_id_index[j] == 1)
          //   _mob[var_index[i]][var_index[j]] = 2.4e-14;
          // else
            _mob[var_index[i]][var_index[j]] = _GBmob_HAGB;

            _mob[var_index[j]][var_index[i]] =  _mob[var_index[i]][var_index[j]];
        }

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
        }  
      } 
  } // When at grain boundaries or junction
  else if (grain_id_index.size() == 1)
    std::fill(_mob.begin(), _mob.end(), std::vector<Real>(_op_num, 0.0));

  // TODO: The following code is only used to identify the activation order parameters‘ number identified by GrainTracker, and will be deleted later
  _twinning_type[_qp] = grain_id_index.size();
}

Real
GBAnisotropyMisori::calculatedGBEnergy(const MisorientationAngleData & misori_s)
{
  Real gbSigma(_GBsigma_HAGB);

  // transition misorientation angle between low and high-angle grain boundary
  Real trans_misori_angle_HAGB = 15.0;

  if (_misori_s._misor <= 1.0)
    gbSigma = _GBsigma_HAGB * ((1.0 / trans_misori_angle_HAGB * (1 - std::log(1.0 / trans_misori_angle_HAGB))));    
  else if (_misori_s._misor <= trans_misori_angle_HAGB)
    gbSigma = _GBsigma_HAGB * ((_misori_s._misor / trans_misori_angle_HAGB * (1 - std::log(_misori_s._misor / trans_misori_angle_HAGB))));

  return gbSigma; 
}

Real
GBAnisotropyMisori::calculatedGBMobility(const MisorientationAngleData & misori_s)
{
  // Initialize GB mobility
  Real gbMob(_GBmob_HAGB);

  // transition misorientation angle between low and high-angle grain boundary
  Real trans_misori_angle_HAGB = 15;

  // Equation constant
  Real B = 5;
  Real n = 4;

  if (_misori_s._misor <= 1.0)
    gbMob = _GBmob_HAGB * ((1- std::exp(-B * std::pow( 1.0 / trans_misori_angle_HAGB, n)))); // Eq.8    
  else if (_misori_s._misor <=  trans_misori_angle_HAGB)
    gbMob = _GBmob_HAGB * ((1- std::exp(-B * std::pow( _misori_s._misor / trans_misori_angle_HAGB, n)))); // Eq.8

  return gbMob;
}