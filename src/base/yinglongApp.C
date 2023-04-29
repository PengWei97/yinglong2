#include "yinglongApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
yinglongApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  
  // Do not use legacy material output, i.e., output properties on INITIAL as well as TIMESTEP_END
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  params.set<bool>("use_legacy_output_syntax") = false;
  
  return params;
}

yinglongApp::yinglongApp(InputParameters parameters) : MooseApp(parameters)
{
  yinglongApp::registerAll(_factory, _action_factory, _syntax);
}

yinglongApp::~yinglongApp() {}

void
yinglongApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAll(f, af, syntax);
  Registry::registerObjectsTo(f, {"yinglongApp"});
  Registry::registerActionsTo(af, {"yinglongApp"});

  /* register custom execute flags, action syntax, etc. here */
  syntax.registerActionSyntax("PolycrystalKernelMOP2Action", "Kernels/PolycrystalKernelMOP2");
  syntax.registerActionSyntax("PolycrystalStoredEnergyEBSDAction", "Kernels/PolycrystalStoredEnergyEBSD");

}

void
yinglongApp::registerApps()
{
  registerApp(yinglongApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
yinglongApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  yinglongApp::registerAll(f, af, s);
}
extern "C" void
yinglongApp__registerApps()
{
  yinglongApp::registerApps();
}
