//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef SURFACE_DEPOSIT2TESTAPP_H
#define SURFACE_DEPOSIT2TESTAPP_H

#include "MooseApp.h"

class SurfaceDeposit2TestApp;

template <>
InputParameters validParams<SurfaceDeposit2TestApp>();

class SurfaceDeposit2TestApp : public MooseApp
{
public:
  SurfaceDeposit2TestApp(InputParameters parameters);
  virtual ~SurfaceDeposit2TestApp();

  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs = false);
};

#endif /* SURFACE_DEPOSIT2TESTAPP_H */
