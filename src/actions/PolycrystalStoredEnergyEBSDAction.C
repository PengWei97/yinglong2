//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PolycrystalStoredEnergyEBSDAction.h"
#include "Factory.h"
#include "Conversion.h"
#include "FEProblem.h"

registerMooseAction("yinglongApp", PolycrystalStoredEnergyEBSDAction, "add_kernel");

InputParameters
PolycrystalStoredEnergyEBSDAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addClassDescription("Action that adds the contribution of stored energy associated with dislocations to grain growth models");
  params.addRequiredParam<unsigned int>("op_num",
                                        "specifies the total number of OPs representing");
  params.addRequiredParam<std::string>("var_name_base", "specifies the base name of the variables");
  params.addParam<Real>("length_scale", 1.0e-9, "Length scale in m, where default is nm");
  params.addParam<bool>(
      "use_displaced_mesh", false, "Whether to use displaced mesh in the kernels");
  params.addRequiredParam<UserObjectName>("grain_tracker",
                                          "The GrainTracker UserObject to get values from.");
  params.addRequiredParam<UserObjectName>("GNDs_provider", "GNDs provider for EBSD reader.");  
  return params;
}

PolycrystalStoredEnergyEBSDAction::PolycrystalStoredEnergyEBSDAction(const InputParameters & params)
  : Action(params),
    _op_num(getParam<unsigned int>("op_num")),
    _var_name_base(getParam<std::string>("var_name_base")),
    _length_scale(getParam<Real>("length_scale"))
{
}

void
PolycrystalStoredEnergyEBSDAction::act()
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
    // Step up ACSEDGPolyEBSD -- Allen-Cahn Stored Energy in Deformed Grains kernels using EBSD
    //

    InputParameters params = _factory.getValidParams("ACSEDGPolyEBSD");
    params.set<NonlinearVariableName>("variable") = var_name;
    params.set<std::vector<VariableName>>("v") = v;
    params.set<Real>("length_scale") = _length_scale;
    params.set<bool>("use_displaced_mesh") = getParam<bool>("use_displaced_mesh");
    params.set<UserObjectName>("grain_tracker") = getParam<UserObjectName>("grain_tracker");
    params.set<UserObjectName>("GNDs_provider") = getParam<UserObjectName>("GNDs_provider");
    params.set<unsigned int>("op_index") = op;

    std::string kernel_name = "ACStoredEnergy_" + var_name;
    _problem->addKernel("ACSEDGPolyEBSD", kernel_name, params);
  }
}
