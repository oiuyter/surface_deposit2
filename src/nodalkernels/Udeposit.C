//doing the RHS of du/dt=v-u

#include "Udeposit.h"
#include "Function.h"
#include "NodalKernel.h"
registerMooseObject("SurfaceDeposit2App", Udeposit);

template <>
InputParameters
validParams<Udeposit>()
{
  InputParameters params = validParams<NodalKernel>();
  params.addCoupledVar("coupled_var", "	Choose the variable you want to couple");
  params.addParam<Real>("K_1", 1, "standard reaction rate");
  params.addRequiredParam<FunctionName>("Exp_3", "Exp_3 term in Butler Volmer Eq");
  params.addRequiredParam<FunctionName>("Exp_4", "Exp_4 term in Butler Volmer Eq");
  return params;
}

Udeposit::Udeposit(const InputParameters & parameters)
  : NodalKernel(parameters), 
    _couple_var(coupledValue("coupled_var")),
    _couple_var_offjac(coupled("coupled_var")),
    _K_1(getParam<Real>("K_1")),
    _func_3(getFunction("Exp_3")),
    _func_4(getFunction("Exp_4"))
{
}


Real
Udeposit::computeQpResidual()
{
    return  -_K_1*(_couple_var[_qp]*_func_3.value(_t, *_current_node) - _u[_qp]*_func_4.value(_t, *_current_node));
}


Real
Udeposit::computeQpJacobian()
{
  return _K_1*_func_4.value(_t, *_current_node);
}

Real
Udeposit::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _couple_var_offjac)
    return -_K_1*_func_3.value(_t, *_current_node);
  else
    return 0;
}
