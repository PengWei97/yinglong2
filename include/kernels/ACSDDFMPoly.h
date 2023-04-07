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
#include "EBSDReader.h"

/** 
 * This kernel adds the contribution of stored energy associated with dislocations to grain growth.
 * This kernel is created based on ACSEDGPoly, which can extract the average 
 * geometrically necessary dislocation densities (GNDs) in EBSD. 
 * The formulation is based on: S. Gentry and K. Thornton, IOP Conf. Series: Materials Science and
 * Engineering 89, 012024, (2015), and other works cited therein.
 * ACSDDFMPoly -- AC - SD(stored) - DF(deformtaion)
 */

// Forward Declarations
class GrainTrackerInterface;

class ACSDDFMPoly : public ACBulk<Real>
{
public:
  static InputParameters validParams();

  ACSDDFMPoly(const InputParameters & parameters);

protected:
  virtual Real computeDFDOP(PFFunctionType type);

  const Real _length_scale;
  const unsigned int _op_num;
  const std::vector<const VariableValue *> _vals;
  const std::vector<unsigned int> _vals_var;
  unsigned int _op_index; // index of the OP the kernel is currently acting on

  const MaterialProperty<Real> & _beta;
  const MaterialProperty<Real> & _rho_eff; // the average/effective dislocation density

  const GrainTrackerInterface & _grain_tracker;
  const EBSDReader & _GNDs_provider;
};
