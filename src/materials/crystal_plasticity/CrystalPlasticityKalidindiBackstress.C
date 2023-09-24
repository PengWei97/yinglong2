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

#include "CrystalPlasticityKalidindiBackstress.h"
#include "libmesh/int_range.h"

registerMooseObject("TensorMechanicsApp", CrystalPlasticityKalidindiBackstress);

InputParameters
CrystalPlasticityKalidindiBackstress::validParams()
{
  InputParameters params = CrystalPlasticityStressUpdateBase::validParams();
  params.addClassDescription("Kalidindi version of homogeneous crystal plasticity.");
  params.addParam<Real>("r", 1.0, "Latent hardening coefficient");
  params.addParam<Real>("h", 541.5, "hardening constants");
  params.addParam<Real>("t_sat", 109.8, "saturated slip system strength");
  params.addParam<Real>("gss_a", 2.5, "coefficient for hardening");
  params.addParam<Real>("ao", 0.001, "slip rate coefficient");
  params.addParam<Real>("xm", 0.1, "exponent for slip rate");
  params.addParam<Real>("gss_initial", 60.8, "initial lattice friction strength of the material");

  params.addParam<MaterialPropertyName>(
      "total_twin_volume_fraction",
      "Total twin volume fraction, if twinning is considered in the simulation");

  params.addParam<Real>("c_backstress",0.0,"c parameter in Armstrong-Frederick backstress evolution");
  params.addParam<Real>("d_backstress",0.0,"d parameter in Armstrong-Frederick backstress evolution");

  return params;
}

CrystalPlasticityKalidindiBackstress::CrystalPlasticityKalidindiBackstress(
    const InputParameters & parameters)
  : CrystalPlasticityStressUpdateBase(parameters),
    // Constitutive values
    _r(getParam<Real>("r")),
    _h(getParam<Real>("h")),
    _tau_sat(getParam<Real>("t_sat")),
    _gss_a(getParam<Real>("gss_a")),
    _ao(getParam<Real>("ao")),
    _xm(getParam<Real>("xm")),
    _gss_initial(getParam<Real>("gss_initial")),

    // resize vectors used in the consititutive slip hardening
    _hb(_number_slip_systems, 0.0),
    _slip_resistance_increment(_number_slip_systems, 0.0),

    // resize local caching vectors used for substepping
    _previous_substep_slip_resistance(_number_slip_systems, 0.0),
    _slip_resistance_before_update(_number_slip_systems, 0.0),

    // Twinning contributions, if used
    _include_twinning_in_Lp(parameters.isParamValid("total_twin_volume_fraction")),
    _twin_volume_fraction_total(_include_twinning_in_Lp
                                    ? &getMaterialPropertyOld<Real>("total_twin_volume_fraction")
                                    : nullptr),

    _c_backstress(getParam<Real>("c_backstress")), // c parameter in Armstrong-Frederick backstress evolution
    _d_backstress(getParam<Real>("d_backstress")), // d parameter in Armstrong-Frederick backstress evolution
    _tau_b(declareProperty<std::vector<Real>>(_base_name + "tau_b")), // Backstress for each slip system
    _tau_b_old(getMaterialPropertyOld<std::vector<Real>>(_base_name + "tau_b"))
{
}

void
CrystalPlasticityKalidindiBackstress::initQpStatefulProperties()
{
  // 初始化_tua, schmit tensor, g^\alpha, gamma^\alpha
  CrystalPlasticityStressUpdateBase::initQpStatefulProperties();

  for (const auto i : make_range(_number_slip_systems))
  {
    _slip_resistance[_qp][i] = _gss_initial; // g^\alpha
    _slip_increment[_qp][i] = 0.0; // gamma^\alpha
  }

  _tau_b[_qp].resize(_number_slip_systems);
  for (const auto i : make_range(_number_slip_systems))
    _tau_b[_qp][i] = 0.0; // initialise back stress

}

void
CrystalPlasticityKalidindiBackstress::setInitialConstitutiveVariableValues()
{
  // Would also set old dislocation densities here if included in this model
  _slip_resistance[_qp] = _slip_resistance_old[_qp];
  _previous_substep_slip_resistance = _slip_resistance_old[_qp];

  _tau_b[_qp] = _tau_b_old[_qp];
}

void
CrystalPlasticityKalidindiBackstress::setSubstepConstitutiveVariableValues() // A
{
  // Would also set substepped dislocation densities here if included in this model
  _slip_resistance[_qp] = _previous_substep_slip_resistance;

  _tau_b[_qp] = _tau_b_old[_qp];
}

bool
CrystalPlasticityKalidindiBackstress::calculateSlipRate()
{
  Real tau_b = 0.0; // Backstress: temporary variable

  for (const auto i : make_range(_number_slip_systems))
  {
    tau_b = _tau_b[_qp][i];

    _slip_increment[_qp][i] =
        _ao * std::pow(std::abs((_tau[_qp][i] - tau_b) / _slip_resistance[_qp][i]), 1.0 / _xm);
    if (_tau[_qp][i] < 0.0)
      _slip_increment[_qp][i] *= -1.0;

    if (std::abs(_slip_increment[_qp][i]) * _substep_dt > _slip_incr_tol)
    {
      if (_print_convergence_message)
        mooseWarning("Maximum allowable slip increment exceeded ",
                     std::abs(_slip_increment[_qp][i]) * _substep_dt);

      return false;
    }
  }
  return true;
}

void
CrystalPlasticityKalidindiBackstress::calculateEquivalentSlipIncrement(
    RankTwoTensor & equivalent_slip_increment)
{
  if (_include_twinning_in_Lp)
  {
    for (const auto i : make_range(_number_slip_systems))
      equivalent_slip_increment += (1.0 - (*_twin_volume_fraction_total)[_qp]) *
                                   _flow_direction[_qp][i] * _slip_increment[_qp][i] * _substep_dt;
  }
  else // if no twinning volume fraction material property supplied, use base class
    CrystalPlasticityStressUpdateBase::calculateEquivalentSlipIncrement(equivalent_slip_increment);
}

void
CrystalPlasticityKalidindiBackstress::calculateConstitutiveSlipDerivative(
    std::vector<Real> & dslip_dtau)
{
  Real tau_b = 0.0; // Backstress: temporary variable

  for (const auto i : make_range(_number_slip_systems))
  {
    tau_b = _tau_b[_qp][i];

    if (MooseUtils::absoluteFuzzyEqual(_tau[_qp][i], 0.0))
      dslip_dtau[i] = 0.0;
    else
      dslip_dtau[i] = _ao / _xm *
                      std::pow((std::abs(_tau[_qp][i] - tau_b) / _slip_resistance[_qp][i]), 1.0 / _xm - 1.0) /
                      _slip_resistance[_qp][i];
  }
}

bool
CrystalPlasticityKalidindiBackstress::areConstitutiveStateVariablesConverged()
{
  return isConstitutiveStateVariableConverged(_slip_resistance[_qp],
                                              _slip_resistance_before_update,
                                              _previous_substep_slip_resistance,
                                              _resistance_tol);
}

void
CrystalPlasticityKalidindiBackstress::updateSubstepConstitutiveVariableValues()
{
  // Would also set substepped dislocation densities here if included in this model
  _previous_substep_slip_resistance = _slip_resistance[_qp];
}

void
CrystalPlasticityKalidindiBackstress::cacheStateVariablesBeforeUpdate()
{
  _slip_resistance_before_update = _slip_resistance[_qp];
}

void
CrystalPlasticityKalidindiBackstress::calculateStateVariableEvolutionRateComponent()
{
  for (const auto i : make_range(_number_slip_systems))
  {
    // Clear out increment from the previous iteration
    _slip_resistance_increment[i] = 0.0;

    _hb[i] = _h * std::pow(std::abs(1.0 - _slip_resistance[_qp][i] / _tau_sat), _gss_a);
    const Real hsign = 1.0 - _slip_resistance[_qp][i] / _tau_sat;
    if (hsign < 0.0)
      _hb[i] *= -1.0;
  }

  for (const auto i : make_range(_number_slip_systems))
  {
    for (const auto j : make_range(_number_slip_systems))
    {
      unsigned int iplane, jplane;
      iplane = i / 3;
      jplane = j / 3;

      if (iplane == jplane) // self vs. latent hardening
        _slip_resistance_increment[i] +=
            std::abs(_slip_increment[_qp][j]) * _hb[j]; // q_{ab} = 1.0 for self hardening
      else
        _slip_resistance_increment[i] +=
            std::abs(_slip_increment[_qp][j]) * _r * _hb[j]; // latent hardenign
    }
  }
}

bool
CrystalPlasticityKalidindiBackstress::updateStateVariables()
{
  // Now perform the check to see if the slip system should be updated
  for (const auto i : make_range(_number_slip_systems))
  {
    _slip_resistance_increment[i] *= _substep_dt;
    if (_previous_substep_slip_resistance[i] < _zero_tol && _slip_resistance_increment[i] < 0.0)
      _slip_resistance[_qp][i] = _previous_substep_slip_resistance[i];
    else
      _slip_resistance[_qp][i] =
          _previous_substep_slip_resistance[i] + _slip_resistance_increment[i];

    if (_slip_resistance[_qp][i] < 0.0)
      return false;
  }

  updateBackstress();
  return true;
}

// Update of the Armstrong-Frederick backstress term:
// Armstrong, P.J., Frederick, C.O., 1966. 
// A Mathematical Representation of the Multiaxial Bauschinger Effect, 
// G.E.G.B. Report RD/B/N. Central Electricity Generating Board
// Backstress term can be changed by changing this function only
void
CrystalPlasticityKalidindiBackstress::updateBackstress()
{
  Real dtau_b; // increment of backstress, temporary variable
  _tau_b[_qp] = _tau_b_old[_qp];

  // _dt here is the substep
  for (const auto i : make_range(_number_slip_systems))
  {	  
    dtau_b = 0.0;
    dtau_b = _c_backstress * _slip_resistance_increment[i];
    dtau_b -= _d_backstress * _tau_b[_qp][i] * std::abs(_slip_resistance_increment[i]);
      
    _tau_b[_qp][i] += dtau_b;
  }
}