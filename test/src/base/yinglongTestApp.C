//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "yinglongTestApp.h"
#include "yinglongApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

InputParameters
yinglongTestApp::validParams()
{
  InputParameters params = yinglongApp::validParams();
  return params;
}

yinglongTestApp::yinglongTestApp(InputParameters parameters) : MooseApp(parameters)
{
  yinglongTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

yinglongTestApp::~yinglongTestApp() {}

void
yinglongTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  yinglongApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"yinglongTestApp"});
    Registry::registerActionsTo(af, {"yinglongTestApp"});
  }
}

void
yinglongTestApp::registerApps()
{
  registerApp(yinglongApp);
  registerApp(yinglongTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
yinglongTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  yinglongTestApp::registerAll(f, af, s);
}
extern "C" void
yinglongTestApp__registerApps()
{
  yinglongTestApp::registerApps();
}
