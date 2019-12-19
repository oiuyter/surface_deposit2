//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef SURFACE_DEPOSIT2APP_H
#define SURFACE_DEPOSIT2APP_H

#include "MooseApp.h"

class SurfaceDeposit2App;

template <>
InputParameters validParams<SurfaceDeposit2App>();

class SurfaceDeposit2App : public MooseApp
{
public:
  SurfaceDeposit2App(InputParameters parameters);
  virtual ~SurfaceDeposit2App();

  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s);
};

#endif /* SURFACE_DEPOSIT2APP_H */
