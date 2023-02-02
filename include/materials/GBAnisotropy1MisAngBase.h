#pragma once

#include "Material.h"
#include "EulerAngleProvider.h"
#include "GrainTrackerGG.h"
#include "CalculateMisorientationAngle.h"

// GBAnisotropy1MisAngBase -- GB energy and GB mobility anisotropy is introduced by calculating the misorientation angle using Quaternion.

/**
 * Function[kappa, gamma, m, L] = parameters (sigma, mob, w_GB, sigma0)
 * Parameter determination method is elaborated in Phys. Rev. B, 78(2), 024113, 2008, by N. Moelans for old MOP-PF model.
**/

class GBAnisotropy1MisAngBase : public Material
{
public:
  static InputParameters validParams();

  GBAnisotropy1MisAngBase(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  // computer GB energy and mobility matrix based on the misorientation
  virtual void computerGBParameter();

  // calculated GB energy based on the the Read-Shockley
  virtual Real calculatedGBEnergy(const misoriAngle_isTwining & misori_gbType);

  // calculated GB mobility based on the sigmoidal law
  virtual Real calculatedGBMobility(const misoriAngle_isTwining & misori_gbType);

  // Calculation of Phase Field Model Parameters Based on GB Energy (sigma_ij) and GB Mobility (mob_ij)
  virtual void computerModelParameter();

  // used to store orientation structure, including misorientation angle, istwinnig, twinning type;
  misoriAngle_isTwining _s_misoriTwin;  

  const GrainTrackerGG & _grain_tracker;
  const EulerAngleProvider & _euler; 
  bool _is_primary;

  const unsigned int _mesh_dimension;

  const Real _length_scale;
  const Real _time_scale;
  const Real _delta_sigma;
  const Real _delta_mob;
  const Real _wGB;

  const Real _matrix_sigma;
  const Real _matrix_mob;
  const Real _matrix_Q;

  const bool _inclination_anisotropy;
  const bool _misorientation_anisotropy;
  const bool _gbMobility_anisotropy;
  const bool _tb_anisotropy; // Whether to consider twin boundary anisotropy

  const VariableValue & _T;

  std::vector<std::vector<Real>> _sigma; // Obtained by the function computerGBParameter()
  std::vector<std::vector<Real>> _mob;
  std::vector<std::vector<Real>> _Q;
  std::vector<std::vector<Real>> _kappa_gamma;
  std::vector<std::vector<Real>> _a_g2;

  MaterialProperty<Real> & _kappa;
  MaterialProperty<Real> & _gamma;
  MaterialProperty<Real> & _L;
  MaterialProperty<Real> & _mu;
  MaterialProperty<Real> & _misAngle;

  MaterialProperty<Real> & _act_wGB; // needless 

  const Real _kb;
  const Real _JtoeV;
  Real _mu_qp;

  const unsigned int _op_num;

  const std::vector<const VariableValue *> _vals;
  const std::vector<const VariableGradient *> _grad_vals;
};
