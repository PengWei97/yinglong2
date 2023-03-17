//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"
#include "EulerAngleProvider.h"
#include "GrainTrackerBase.h"

// Forward Declarations
class GrainTrackerBase;
class EulerAngleProvider;

/**
 * Output euler angles from user object to an AuxVariable.
 */
class OutputEulerAnglesGG : public AuxKernel
{
public:
  static InputParameters validParams();

  OutputEulerAnglesGG(const InputParameters & parameters);

protected:
  virtual Real computeValue();
  virtual void precalculateValue();

  /// Object providing the Euler angles
  const EulerAngleProvider & _euler;

  /// Grain tracker object
  const GrainTrackerBase & _grain_tracker;

  /// Number of grains
  MooseEnum _output_euler_angle;

  /// precalculated element value
  Real _value;
};
