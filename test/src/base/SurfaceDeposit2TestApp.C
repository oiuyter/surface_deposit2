//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "SurfaceDeposit2TestApp.h"
#include "SurfaceDeposit2App.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

template <>
InputParameters
validParams<SurfaceDeposit2TestApp>()
{
  InputParameters params = validParams<SurfaceDeposit2App>();
  return params;
}

SurfaceDeposit2TestApp::SurfaceDeposit2TestApp(InputParameters parameters) : MooseApp(parameters)
{
  SurfaceDeposit2TestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

SurfaceDeposit2TestApp::~SurfaceDeposit2TestApp() {}

void
SurfaceDeposit2TestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  SurfaceDeposit2App::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"SurfaceDeposit2TestApp"});
    Registry::registerActionsTo(af, {"SurfaceDeposit2TestApp"});
  }
}

void
SurfaceDeposit2TestApp::registerApps()
{
  registerApp(SurfaceDeposit2App);
  registerApp(SurfaceDeposit2TestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
SurfaceDeposit2TestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  SurfaceDeposit2TestApp::registerAll(f, af, s);
}
extern "C" void
SurfaceDeposit2TestApp__registerApps()
{
  SurfaceDeposit2TestApp::registerApps();
}
