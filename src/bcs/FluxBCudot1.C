// operate on C_R or C_O, in coupling term they have the same form.
#include "FluxBCudot1.h"
#include "Function.h"

registerMooseObject("SurfaceDeposit2App", FluxBCudot1);

template <>
InputParameters
validParams<FluxBCudot1>()
{
  InputParameters params = validParams<IntegratedBC>();

  params.addCoupledVar("coupled_var", "	Choose the variable you want to couple");
  params.addCoupledVar("coupled_var1", " Choose the variable you want to couple");
  params.addParam<Real>("K_O", 1, "standard reaction rate");
  params.addParam<Real>("K_1", 1, "standard reaction rate");
  params.addParam<Real>("sigma1", 0.01, "deposition layer");
  params.addRequiredParam<FunctionName>("Exp_1", "Exp_1");
  params.addRequiredParam<FunctionName>("Exp_2", "Exp_2");
  params.addRequiredParam<FunctionName>("Exp_3", "Exp_3");
  params.addRequiredParam<FunctionName>("Exp_4", "Exp_4 ");

  return params;
}

FluxBCudot1::FluxBCudot1(const InputParameters & parameters)
  : IntegratedBC(parameters),
    _couple_var(coupledValue("coupled_var")),
    _couple_var_offjac(coupled("coupled_var")),
    _couple_var1(coupledValue("coupled_var1")),
    _couple_var1_offjac(coupled("coupled_var1")),
    _u_dot(dot()),
    _du_dot_du(dotDu()),
    _K_O(getParam<Real>("K_O")),
    _K_1(getParam<Real>("K_1")),
    _sigma1(getParam<Real>("sigma1")),
    _func_1(getFunction("Exp_1")),
    _func_2(getFunction("Exp_2")),
    _func_3(getFunction("Exp_3")),
    _func_4(getFunction("Exp_4"))
{
}

Real
FluxBCudot1::computeQpResidual()
{
    return  _test[_i][_qp] * (_u_dot[_qp]*_sigma1 + _K_1*(_u[_qp]*_func_3.value(_t, _q_point[_qp])-_couple_var[_qp] * _func_4.value(_t, _q_point[_qp]))-_K_O*(_couple_var1[_qp]*_func_1.value(_t, _q_point[_qp])-_u[_qp]* _func_2.value(_t, _q_point[_qp])));
}


Real
FluxBCudot1::computeQpJacobian()
{
  return _test[_i][_qp] * (_du_dot_du[_qp]*_sigma1 + _K_O*_func_2.value(_t, _q_point[_qp])+_K_1*_func_3.value(_t, _q_point[_qp])) * _phi[_j][_qp];
}

Real
FluxBCudot1::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _couple_var_offjac)
    return _test[_i][_qp] * (-_K_1*_func_4.value(_t, _q_point[_qp])*_phi[_j][_qp]);

  if (jvar == _couple_var1_offjac)
    return _test[_i][_qp] * (-_K_O*_func_1.value(_t, _q_point[_qp])*_phi[_j][_qp]);
  else
    return 0;
}


