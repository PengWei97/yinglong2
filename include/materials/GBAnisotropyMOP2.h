#pragma once

#include "Material.h"

/**
 * Function[kappa, gamma_ij, m, L_ij] = parameters (sigma_ij, mob_ij, w_GB_init, sigma0)
 * Parameter determination method is elaborated in Mater. Des. 217 (2022) 110592, by N. Moelans for new MOP-PF model (MD2022-new)
 * Thanks to Prof. Moelans for the explanation of her paper.
 */

class GBAnisotropyMOP2 : public Material
{
public:
  static InputParameters validParams();

  GBAnisotropyMOP2(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  // Get GB anisotropy based on file
  virtual void getGBAnisotropyFromFile();

  // Calculation of Phase Field Model Parameters Based on GB Energy (sigma_ij) and GB Mobility (mob_ij)
  virtual void computerModelParameter();

  const unsigned int _mesh_dimension;

  const Real _length_scale;
  const Real _time_scale;
  const Real _delta_sigma;
  const Real _delta_mob;
  const Real _wGB_init;

  const FileName _Anisotropic_GB_file_name;

  const bool _inclination_anisotropy;

  const VariableValue & _T;

  std::vector<std::vector<Real>> _sigma;
  std::vector<std::vector<Real>> _mob;
  std::vector<std::vector<Real>> _Q;

  std::vector<std::vector<Real>> _L_kappa; 
  std::vector<std::vector<Real>> _wGB_g2;

  // four model parameters of the MOP-PF model
  MaterialProperty<Real> & _L; 
  MaterialProperty<Real> & _kappa; 
  MaterialProperty<std::vector<std::vector<Real>>> & _gamma_ij;
  MaterialProperty<Real> & _mu; 

  MaterialProperty<Real> & _act_wGB; // may need to be deleted

  const Real _kb;
  const Real _JtoeV;
  Real _mu_qp;
  Real _kappa_qp;

  const unsigned int _op_num;

  const std::vector<const VariableValue *> _vals;
  const std::vector<const VariableGradient *> _grad_vals;
};
