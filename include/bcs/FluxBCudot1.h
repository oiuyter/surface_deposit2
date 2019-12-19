#ifndef FLUXBCUDOT1_H
#define FLUXBCUDOT1_H

#include "IntegratedBC.h"

// Forward Declarations
class FluxBCudot1;
class Function;

template <>
InputParameters validParams<FluxBCudot1>();

/**
 * Implements Neumann BC where grad(u)=udot-something on the boundary.
 */

class FluxBCudot1 : public IntegratedBC

{
public:

  FluxBCudot1(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:

  const VariableValue & _couple_var;
  const unsigned _couple_var_offjac;
  const VariableValue & _couple_var1;
  const unsigned _couple_var1_offjac;
  const VariableValue & _u_dot;
  const VariableValue & _du_dot_du;
  Real _K_O;
  Real _K_1;
  Real _sigma1;
  Function & _func_1;
  Function & _func_2;
  Function & _func_3;
  Function & _func_4;
};

#endif
