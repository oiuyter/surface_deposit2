//revised from ConstantRate

#ifndef UDEPOSIT_H
#define UDEPOSIT_H

#include "NodalKernel.h"

// Forward Declarations
class Udeposit;
class Function;
//class NodalKernel;
template <>
InputParameters validParams<Udeposit>();

class Udeposit : public NodalKernel
{
public:
  /**
   * Constructor initializes the rate
   */
  Udeposit(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
  const VariableValue & _couple_var;
  const unsigned _couple_var_offjac;
    Real _K_1;
  //const MaterialProperty<Real> & _K_O;

  Function & _func_3;
 
  Function & _func_4;

};

#endif
