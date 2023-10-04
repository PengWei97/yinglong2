/*
   Author: pengwei
   Date: 2023-09-20
   Constitutive model form:
   Reference: Kalidindi, S. R., Bronkhorst, C. A., & Anand, L. (1992). Crystallographic texture evolution in bulk deformation processing of FCC metals. Journal of the Mechanics and Physics of Solids, 40(3), 537-569.

   Plus an Armstrong-Frederick backstress term:
   Reference: Armstrong, P. J., & Frederick, C. O. (1966). A mathematical representation of the multiaxial Bauschinger effect (Vol. 731). Berkeley, CA: Berkeley Nuclear Laboratories.

   The flow rule of the shearing slip rate:
   \dot{\gamma}^\alpha = \dot{\gamma}_o * |(\tau^\alpha - \chi^\alpha) / g^\alpha|^M * sgn(\tau^\alpha - \chi^\alpha)

   Armstrong-Frederick backstress term:
   \dot{\chi}^\alpha = h * \dot{\gamma}^\alpha - r * |\dot{\gamma}^\alpha| * \chi^\alpha

   Reference: https://github.com/ngrilli/c_pfor_am/blob/main/src/materials/CrystalPlasticityDislocationUpdate.C
*/

#pragma once

#include "CrystalPlasticityKalidindiUpdate.h"

class CrystalPlasticityKalidindiBackstress;

class CrystalPlasticityKalidindiBackstress : public CrystalPlasticityKalidindiUpdate
{
public:
  static InputParameters validParams();

  CrystalPlasticityKalidindiBackstress(const InputParameters & parameters);

protected:
  /**
   * initializes the stateful properties such as
   * stress, plastic deformation gradient, slip system resistances, backstress etc.
   */
  virtual void initQpStatefulProperties() override;

  /**
   * Sets the value of the current and previous substep iteration slip system
   * resistance to the old value at the start of the PK2 stress convergence
   * while loop.
   */
  virtual void setInitialConstitutiveVariableValues() override;

  /**
   * Sets the current slip system resistance value to the previous substep value.
   * In cases where only one substep is taken (or when the first) substep is taken,
   * this method just sets the current value to the old slip system resistance
   * value again.
   */
  virtual void setSubstepConstitutiveVariableValues() override;

  virtual bool calculateSlipRate() override; 

  virtual void calculateConstitutiveSlipDerivative(std::vector<Real> & dslip_dtau) override;

  virtual bool areConstitutiveStateVariablesConverged() override;

  virtual void updateSubstepConstitutiveVariableValues() override;

  // Cache the slip system value before the update for the diff in the convergence check
  virtual void cacheStateVariablesBeforeUpdate() override;

  virtual void calculateStateVariableEvolutionRateComponent() override;

  /*
   * Finalizes the values of the state variables and slip system resistance
   * for the current timestep after convergence has been reached.
   */
  virtual bool updateStateVariables() override;

  // This function updates the backstress
  virtual void ArmstrongFrederickBackstressUpdate();

  const Real _h;
  const Real _h_D;

  // Backstress for each slip system
  MaterialProperty<std::vector<Real>> & _backstress;
  const MaterialProperty<std::vector<Real>> & _backstress_old;
  std::vector<Real> _backstress_increment;

  std::vector<Real> _previous_substep_backstress;
  std::vector<Real> _backstress_before_update;
};
