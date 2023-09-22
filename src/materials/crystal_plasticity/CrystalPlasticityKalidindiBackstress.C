
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

registerMooseObject("yinglongApp", CrystalPlasticityKalidindiBackstress);

InputParameters
CrystalPlasticityKalidindiBackstress::validParams()
{
  InputParameters params = CrystalPlasticityKalidindiUpdate::validParams();
  params.addClassDescription("Kalidindi version of homogeneous crystal plasticity considering back stress for fatigure loading.");
  params.addParam<Real>("c_backstress",0.0,"c parameter in Armstrong-Frederick backstress evolution");
  params.addParam<Real>("d_backstress",0.0,"d parameter in Armstrong-Frederick backstress evolution");

  return params;
}

CrystalPlasticityKalidindiBackstress::CrystalPlasticityKalidindiBackstress(
    const InputParameters & parameters)
  : CrystalPlasticityKalidindiUpdate(parameters),
    _c_backstress(getParam<Real>("c_backstress")), // c parameter in Armstrong-Frederick backstress evolution
    _d_backstress(getParam<Real>("d_backstress")), // d parameter in Armstrong-Frederick backstress evolution
    _tau_b(declareProperty<std::vector<Real>>("tau_b")), // Backstress for each slip system
    _tau_b_old(getMaterialPropertyOld<std::vector<Real>>("tau_b"))
{
}

void
CrystalPlasticityKalidindiBackstress::initQpStatefulProperties()
{
  // tau, flow_direction, slip_resistance(g), slip_increment(delta_gamma^alpha), _tau_b
  CrystalPlasticityKalidindiUpdate::initQpStatefulProperties();

  _tau_b[_qp].resize(_number_slip_systems);
  for (const auto i : make_range(_number_slip_systems))
    _tau_b[_qp][i] = 0.0; // initialise back stress
}

void
CrystalPlasticityKalidindiBackstress::setInitialConstitutiveVariableValues()
{
 CrystalPlasticityKalidindiUpdate::setInitialConstitutiveVariableValues();

  _tau_b[_qp] = _tau_b_old[_qp];
}

void
CrystalPlasticityKalidindiBackstress::setSubstepConstitutiveVariableValues()
{
  CrystalPlasticityKalidindiUpdate::setSubstepConstitutiveVariableValues();

  _tau_b[_qp] = _tau_b_old[_qp];
}

bool
CrystalPlasticityKalidindiBackstress::calculateSlipRate()
{
  Real tau_b; // Backstress: temporary variable

  for (const auto i : make_range(_number_slip_systems))
  {
    tau_b = _tau_b[_qp][i];

    _slip_increment[_qp][i] =
        _ao * std::pow(std::abs((_tau[_qp][i] - tau_b) / _slip_resistance[_qp][i]), 1.0 / _xm) * std::copysign(1.0, (_tau[_qp][i] - tau_b));

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
CrystalPlasticityKalidindiBackstress::calculateConstitutiveSlipDerivative(
    std::vector<Real> & dslip_dtau)
{
  Real tau_b; // Backstress: temporary variable

  for (const auto i : make_range(_number_slip_systems))
  {
    tau_b = _tau_b[_qp][i];

    if (MooseUtils::absoluteFuzzyEqual(_tau[_qp][i], 0.0))
      dslip_dtau[i] = 0.0;
    else
      dslip_dtau[i] = _ao / _xm * 
                      std::pow(std::abs((_tau[_qp][i] - tau_b) / _slip_resistance[_qp][i]), 1.0 / _xm - 1.0) /
                      _slip_resistance[_qp][i];
  }
}

// void
// CrystalPlasticityKalidindiBackstress::updateSubstepConstitutiveVariableValues()
// {
//   // Would also set substepped dislocation densities here if included in this model
//   _previous_substep_slip_resistance = _slip_resistance[_qp];
// }

// void
// CrystalPlasticityKalidindiBackstress::cacheStateVariablesBeforeUpdate()
// {
//   _slip_resistance_before_update = _slip_resistance[_qp];
// }

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

    updateBackstress();

    if (_slip_resistance[_qp][i] < 0.0)
      return false;
  }
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

  for (const auto i : make_range(_number_slip_systems))
    _tau_b[i] = _tau_b_old[i];

  // _dt here is the substep
  for (const auto i : make_range(_number_slip_systems))
  {	  
    dtau_b = 0.0;
    dtau_b = _c_backstress * _slip_resistance_increment[i];
    dtau_b -= _d_backstress * _tau_b[_qp][i] * std::abs(_slip_resistance_increment[i]);
      
    _tau_b[_qp][i] += dtau_b;
  }
}
