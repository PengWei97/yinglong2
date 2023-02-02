#pragma once

#include "Material.h"

// GBAnisotropy1MisAng -- This version of GB energy and GB mobility anisotropy is introduced by the txt file.

/**
 * Function[kappa, gamma, m, L] = parameters (sigma, mob, w_GB, sigma0)
 * Parameter determination method is elaborated in Phys. Rev. B, 78(2), 024113, 2008, by N. Moelans for old MOP-PF model
 * Thanks to Prof. Moelans for the explanation of her paper.
 */

class GBAnisotropy1MisAng : public Material
{
public:
  static InputParameters validParams();

  GBAnisotropy1MisAng(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  // Get GB anisotropy based on file
  virtual void getGBAnisotropyFromFile();

  // Calculation of Phase Field Model Parameters Based on Grain Boundary Energy (sigma_ij) and Grain Boundary Mobility (mob_ij)
  virtual void computerModelParameter();

  const unsigned int _mesh_dimension;

  const Real _length_scale;
  const Real _time_scale;
  const Real _delta_sigma;
  const Real _delta_mob;
  const Real _wGB;

  const FileName _Anisotropic_GB_file_name;

  const bool _inclination_anisotropy;

  const VariableValue & _T;

  std::vector<std::vector<Real>> _sigma;
  std::vector<std::vector<Real>> _mob;
  std::vector<std::vector<Real>> _Q;
  std::vector<std::vector<Real>> _kappa_gamma;
  std::vector<std::vector<Real>> _a_g2;

  MaterialProperty<Real> & _kappa;
  MaterialProperty<Real> & _gamma;
  MaterialProperty<Real> & _L;
  MaterialProperty<Real> & _mu;

  MaterialProperty<Real> & _act_wGB; // needless 

  const Real _kb;
  const Real _JtoeV;
  Real _mu_qp;

  const unsigned int _op_num;

  const std::vector<const VariableValue *> _vals;
  const std::vector<const VariableGradient *> _grad_vals;
};
