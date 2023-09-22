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

#include "CrystalPlasticityKalidindiUpdate.h"

class CrystalPlasticityKalidindiUpdate;

class CrystalPlasticityKalidindiBackstress : public CrystalPlasticityKalidindiUpdate
{
public:
  static InputParameters validParams();

  CrystalPlasticityKalidindiBackstress(const InputParameters & parameters);

protected:
  /**
   * initializes the stateful properties such as
   * stress, plastic deformation gradient, slip system resistances, back stress, etc.
   */
  virtual void initQpStatefulProperties() override;

  virtual void setInitialConstitutiveVariableValues() override;

  virtual void setSubstepConstitutiveVariableValues() override;

  virtual bool calculateSlipRate() override;

  virtual void calculateConstitutiveSlipDerivative(std::vector<Real> & dslip_dtau) override;

  virtual bool updateStateVariables() override;

  // This function updates the backstress
  virtual void updateBackstress();

  const Real _c_backstress;
  const Real _d_backstress;

  // Backstress for each slip system
  MaterialProperty<std::vector<Real>> & _tau_b;
  const MaterialProperty<std::vector<Real>> & _tau_b_old;
};
