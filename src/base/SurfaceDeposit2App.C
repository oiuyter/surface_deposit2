#include "SurfaceDeposit2App.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

template <>
InputParameters
validParams<SurfaceDeposit2App>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

SurfaceDeposit2App::SurfaceDeposit2App(InputParameters parameters) : MooseApp(parameters)
{
  SurfaceDeposit2App::registerAll(_factory, _action_factory, _syntax);
}

SurfaceDeposit2App::~SurfaceDeposit2App() {}

void
SurfaceDeposit2App::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"SurfaceDeposit2App"});
  Registry::registerActionsTo(af, {"SurfaceDeposit2App"});

  /* register custom execute flags, action syntax, etc. here */
}

void
SurfaceDeposit2App::registerApps()
{
  registerApp(SurfaceDeposit2App);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
SurfaceDeposit2App__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  SurfaceDeposit2App::registerAll(f, af, s);
}
extern "C" void
SurfaceDeposit2App__registerApps()
{
  SurfaceDeposit2App::registerApps();
}
