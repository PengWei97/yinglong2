//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ACBulk.h"
#include "Material.h"
#include "GrainTrackerInterface.h"
#include "EBSDReader.h"

/**
 * This kernel (based ACSEDGPolyEBSD) adds the contribution of stored energy associated 
 * with dislocations to grain growth. The formulation is based on:
 * S. Gentry and K. Thornton, IOP Conf. Series: Materials Science and
 * Engineering 89, 012024, (2015), and other works cited therein.
 * This kernel depends grain_index instead of op_index to take full advantage of Grain Tracker
 * So a grain_tracker UserObject must be used. 
 * Meanwhile, EBSDReader is used to provide GNDs.
 */

// Forward Declarations
class GrainTrackerInterface;

class ACSEDGPolyEBSD : public ACBulk<Real>
{
public:
  static InputParameters validParams();

  ACSEDGPolyEBSD(const InputParameters & parameters);

protected:
  virtual Real computeDFDOP(PFFunctionType type);

  virtual Real getGNDsFromEBSD(const unsigned int & grain_id);

  const unsigned int _op_num;

  const std::vector<const VariableValue *> _vals;
  const std::vector<unsigned int> _vals_var;

  const Real _length_scale;
  const bool _concurrent_recovery;
  /// the prefactor needed to calculate the deformation energy from dislocation density
  const MaterialProperty<Real> & _beta;

  /// the average/effective dislocation density
  const MaterialProperty<Real> & _rho_eff;

  /// Grain tracker object
  const GrainTrackerInterface & _grain_tracker;

  /// Grain tracker object
  const EBSDReader & _GNDs_provider;

  /// index of the OP the kernel is currently acting on
  unsigned int _op_index;
};
