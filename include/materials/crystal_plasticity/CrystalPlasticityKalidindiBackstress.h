// pengwei ECUST 2023.09.20

// Constitutive model form:
// Kalidindi, S. R., Bronkhorst, C. A., & Anand, L. (1992). Crystallographic texture evolution in bulk deformation processing of FCC metals. Journal of the Mechanics and Physics of Solids, 40(3), 537-569.

// Plus an Armstrong-Frederick backstress term:
// Armstrong, P. J., & Frederick, C. O. (1966). A mathematical representation of the multiaxial Bauschinger effect (Vol. 731). Berkeley, CA: Berkeley Nuclear Laboratories.

// the flow rule of the shearing slip rate\:
// \dot{\gamma}^\alpha=\dot{\gamma}_o\left|\frac{\tau^\alpha-\chi^\alpha}{g^\alpha}\right|^M \operatorname{sgn}\left(\tau^\alpha-\chi^\alpha\right)

// Armstrong-Frederick backstress term:
// \dot{\chi}^\alpha=h \dot{\gamma}^\alpha-r \chi^\alpha\left|\dot{\gamma}^\alpha\right|

// https://github.com/ngrilli/c_pfor_am/blame/main/src/materials/FiniteStrainCrystalPlasticityBackstress.C

#pragma once

#include "CrystalPlasticityStressUpdateBase.h"

class CrystalPlasticityKalidindiBackstress;

/**
 * CrystalPlasticityKalidindiBackstress uses the multiplicative decomposition of the
 * deformation gradient and solves the PK2 stress residual equation at the
 * intermediate configuration to evolve the material state. The internal
 * variables are updated using an interative predictor-corrector algorithm.
 * Backward Euler integration rule is used for the rate equations.
 */

class CrystalPlasticityKalidindiBackstress : public CrystalPlasticityStressUpdateBase
{
public:
  static InputParameters validParams();

  CrystalPlasticityKalidindiBackstress(const InputParameters & parameters);

protected:
  /**
   * initializes the stateful properties such as
   * stress, plastic deformation gradient, slip system resistances, etc.
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

  /**
   * Stores the current value of the slip system resistance into a separate
   * material property in case substepping is needed.
   */
  virtual void updateSubstepConstitutiveVariableValues() override;

  virtual bool calculateSlipRate() override;

  virtual void
  calculateEquivalentSlipIncrement(RankTwoTensor & /*equivalent_slip_increment*/) override;

  virtual void calculateConstitutiveSlipDerivative(std::vector<Real> & dslip_dtau) override;

  // Cache the slip system value before the update for the diff in the convergence check
  virtual void cacheStateVariablesBeforeUpdate() override;

  /**
   * Following the Constitutive model for slip system resistance as given in
   * Kalidindi, S.R., C.A. Bronkhorst, and L. Anand. Crystallographic texture
   * evolution in bulk deformation processing of FCC metals. Journal of the
   * Mechanics and Physics of Solids 40, no. 3 (1992): 537-569. Eqns 40 - 43.
   * The slip system resistant increment is calculated as
   * $\Delta g = \left| \Delta \gamma \cdot q^{\alpha \beta} \cdot h^{\beta} \right|$
   * and a convergence check is performed on the slip system resistance increment
   */
  virtual void calculateStateVariableEvolutionRateComponent() override;

  /*
   * Finalizes the values of the state variables and slip system resistance
   * for the current timestep after convergence has been reached.
   */
  virtual bool updateStateVariables() override;

  /*
   * Determines if the state variables, e.g. defect densities, have converged
   * by comparing the change in the values over the iteration period.
   */
  virtual bool areConstitutiveStateVariablesConverged() override;

  // This function updates the backstress
  virtual void updateBackstress();

  ///@{Varibles used in the Kalidindi 1992 slip system resistance constiutive model
  const Real _r;
  const Real _h;
  const Real _tau_sat;
  const Real _gss_a;
  const Real _ao;
  const Real _xm;
  const Real _gss_initial;
  ///@}

  /**
   * Slip system interaction matrix used to calculate the hardening contributions
   * from the self and latent slip systems, from Kalidindi et al (1992).
   */
  std::vector<Real> _hb;

  /// Increment of increased resistance for each slip system
  std::vector<Real> _slip_resistance_increment;

  /**
   * Stores the values of the slip system resistance from the previous substep
   * In classes which use dislocation densities, analogous dislocation density
   * substep vectors will be required.
   */
  std::vector<Real> _previous_substep_slip_resistance;

  /**
   * Caches the value of the current slip system resistance immediately prior
   * to the update of the slip system resistance, and is used to calculate the
   * the slip system resistance increment for the current substep (or step if
   * only one substep is taken) for the convergence check tolerance comparison.
   * In classes which use dislocation densities, analogous dislocation density
   * caching vectors will also be required.
   */
  std::vector<Real> _slip_resistance_before_update;

  /**
   * Flag to include the total twin volume fraction in the plastic velocity
   * gradient calculation, per Kalidindi IJP (2001).
   */
  const bool _include_twinning_in_Lp;

  /**
   * User-defined material property name for the total volume fraction of twins
   * in a twinning propagation constitutive model, when this class is used in
   * conjunction with the twinning propagation model.
   * Note that this value is the OLD material property and thus lags the current
   * value by a single timestep.
   */
  const MaterialProperty<Real> * const _twin_volume_fraction_total;

  const Real _c_backstress;
  const Real _d_backstress;

  // Backstress for each slip system
  MaterialProperty<std::vector<Real>> & _tau_b;
  const MaterialProperty<std::vector<Real>> & _tau_b_old;

  // std::vector<Real> _previous_substep_tau_b;
};
