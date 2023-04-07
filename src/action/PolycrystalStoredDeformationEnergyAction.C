//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PolycrystalStoredDeformationEnergyAction.h"
#include "Factory.h"
#include "Conversion.h"
#include "FEProblem.h"

registerMooseAction("PhaseFieldApp", PolycrystalStoredDeformationEnergyAction, "add_kernel");

InputParameters
PolycrystalStoredDeformationEnergyAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addClassDescription("Action that adds the contribution of stored energy associated with "
                             "dislocations to grain growth models");
  params.addRequiredParam<unsigned int>("op_num",
                                        "specifies the total number of OPs representing "
                                        "all grains (deformed + undeformed "
                                        "(recrystallized)) to create");
  params.addRequiredParam<std::string>("var_name_base", "specifies the base name of the variables");
  params.addRequiredParam<Real>("length_scale", "Length scale in m, where default is nm");
  params.addParam<VariableName>("c", "Name of coupled concentration variable");
  params.addParam<VariableName>("T", "Name of temperature variable");
  params.addParam<bool>(
      "use_displaced_mesh", false, "Whether to use displaced mesh in the kernels");
  params.addRequiredParam<UserObjectName>("grain_tracker",
                                          "The GrainTracker UserObject to get values from.");
  params.addRequiredParam<UserObjectName>("GNDs_provider", "Geometric necessary dislocation densities (GNDs) provider for EBSD reader");                                              
  return params;
}

PolycrystalStoredDeformationEnergyAction::PolycrystalStoredDeformationEnergyAction(const InputParameters & params)
  : Action(params),
    _op_num(getParam<unsigned int>("op_num")),
    _length_scale(getParam<Real>("length_scale")),
    _var_name_base(getParam<std::string>("var_name_base"))
{
}

void
PolycrystalStoredDeformationEnergyAction::act()
{
  for (unsigned int op = 0; op < _op_num; ++op)
  {
    //
    // Create variable names
    //

    std::string var_name = _var_name_base + Moose::stringify(op);
    std::vector<VariableName> v;
    v.resize(_op_num - 1);

    unsigned int ind = 0;
    for (unsigned int j = 0; j < _op_num; ++j)
      if (j != op)
        v[ind++] = _var_name_base + Moose::stringify(j);

    //
    // Set up ACSDDFMPoly Stored Energy in Deformed Grains kernels
    //

    InputParameters params = _factory.getValidParams("ACSDDFMPoly");
    params.set<NonlinearVariableName>("variable") = var_name;
    params.set<std::vector<VariableName>>("v") = v;
    params.set<Real>("length_scale") = _length_scale;
    params.set<UserObjectName>("grain_tracker") = getParam<UserObjectName>("grain_tracker");
    params.set<UserObjectName>("GNDs_provider") = getParam<UserObjectName>("GNDs_provider");
    params.set<bool>("use_displaced_mesh") = getParam<bool>("use_displaced_mesh");
    params.set<unsigned int>("op_index") = op;

    std::string kernel_name = "ACStoredEnergy_" + var_name;
    _problem->addKernel("ACSDDFMPoly", kernel_name, params);
  }
}
