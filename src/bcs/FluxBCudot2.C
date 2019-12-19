//operate on C_R or C_O, in coupling term they have the same form.
#include "FluxBCudot2.h"
#include "Function.h"

registerMooseObject("SurfaceDeposit2App", FluxBCudot2);

template <>
InputParameters
validParams<FluxBCudot2>()
{
  InputParameters params = validParams<IntegratedBC>();

  params.addCoupledVar("coupled_var", "	Choose the variable you want to couple");
  params.addParam<Real>("K_L", 1,".");
  params.addParam<Real>("K_D", 1,".");
  params.addParam<Real>("G_S", 1,".");
  params.addParam<Real>("kf_1", 1,".");
  params.addParam<Real>("kb_1", 1,".");
  params.addParam<Real>("kf_2", 1,".");
  params.addParam<Real>("kb_2", 1,".");
  params.addParam<Real>("sigma", 0.01, "deposition layer");
  params.addRequiredParam<FunctionName>("Exp_1", "Exp_1 ");
  params.addRequiredParam<FunctionName>("Exp_2", "Exp_2 ");
  params.addRequiredParam<FunctionName>("Exp_3", "Exp_3 ");
  params.addRequiredParam<FunctionName>("Exp_4", "Exp_4 ");

  return params;
}

FluxBCudot2::FluxBCudot2(const InputParameters & parameters)
  : IntegratedBC(parameters),
    _couple_var(coupledValue("coupled_var")),
    _couple_var_offjac(coupled("coupled_var")),
    _u_dot(dot()),
    _du_dot_du(dotDu()),
    _K_L(getParam<Real>("K_L")),
    _K_D(getParam<Real>("K_D")),
    _G_S(getParam<Real>("G_S")),
    _kf_1(getParam<Real>("kf_1")),
    _kb_1(getParam<Real>("kb_1")),
    _kf_2(getParam<Real>("kf_2")),
    _kb_2(getParam<Real>("kb_2")),
    _func_1(getFunction("Exp_1")),
    _func_2(getFunction("Exp_2")),
    _func_3(getFunction("Exp_3")),
    _func_4(getFunction("Exp_4"))
{
}

Real
FluxBCudot2::theta()
{
  return (_K_L * _u[_qp]) / (( - _K_D * _K_L + _K_D * _K_D) * _u[_qp] * _u[_qp] + (_K_L - 2 * _K_D) * _u[_qp] + 1);
}

Real
FluxBCudot2::gamma()
{
  return (_u_dot[_qp] * _K_L * _G_S) * ((( - _K_D * _K_L + _K_D * _K_D) * _u[_qp] * _u[_qp] + (_K_L - 2 * _K_D) * _u[_qp] + 1) - _u[_qp] * (2 * _K_D * (_K_L - _K_D) * _u[_qp] + _K_L - 2 * _K_D)) / ((( - _K_D * _K_L + _K_D * _K_D) * _u[_qp] * _u[_qp] + (_K_L - 2 * _K_D) * _u[_qp] + 1) * (( - _K_D * _K_L + _K_D * _K_D) * _u[_qp] * _u[_qp] + (_K_L - 2 * _K_D) * _u[_qp] + 1));
}

Real
FluxBCudot2::diftheta()
{
  return (_K_L*(-(_K_D * _K_L - _K_D * _K_D)*_u[_qp]*_u[_qp]+(_K_L-2*_K_D)*_u[_qp]+1) + (_K_L*_u[_qp])*(-2*(_K_D * _K_L - _K_D * _K_D)*_u[_qp]+(_K_L-2*_K_D)))/((-(_K_D * _K_L - _K_D * _K_D)*_u[_qp]*_u[_qp]+(_K_L-2*_K_D)*_u[_qp]+1)*(-(_K_D * _K_L - _K_D * _K_D)*_u[_qp]*_u[_qp]+(_K_L-2*_K_D)*_u[_qp]+1));
}

Real
FluxBCudot2::difgamma()
{
  return ((_du_dot_du[_qp] * _K_L * _G_S)  *  ((_K_D * _K_D - _K_D * _K_L * _u[_qp] * _u[_qp] + (_K_L + 2 * _K_D) * _u[_qp] + 1) - _u[_qp] * (2 * _K_D * _K_D - _K_D * _K_L * _u[_qp] + _K_L - 2 * _K_D)) + ((_u_dot[_qp] * _K_L * _G_S) * - _u[_qp]  *  (2  *  _K_D  *  (_K_D - _K_L))) * ((_K_D * _K_D - _K_D * _K_L * _u[_qp] * _u[_qp] + (_K_L + 2 * _K_D) * _u[_qp] + 1) * (_K_D * _K_D - _K_D * _K_L * _u[_qp] * _u[_qp] + (_K_L + 2 * _K_D) * _u[_qp] + 1)) + 2 * (_K_D * _K_D - _K_D * _K_L * _u[_qp] * _u[_qp] - (_K_L + 2 * _K_D) * _u[_qp] + 1) * (2 * _K_D * _K_D - _K_D * _K_L * _u[_qp] + _K_L - 2 * _K_D) * (_u_dot[_qp] * _K_L * _G_S) * (( - (_K_D * _K_L - _K_D * _K_D) * _u[_qp] * _u[_qp] + (_K_L - 2 * _K_D) * _u[_qp] + 1) - _u[_qp] * (2 * _K_D * (_K_L - _K_D) * _u[_qp] + _K_L - 2 * _K_D))/(( - (_K_D * _K_L - _K_D * _K_D) * _u[_qp] * _u[_qp] + (_K_L - 2 * _K_D) * _u[_qp] + 1) * ( - (_K_D * _K_L - _K_D * _K_D) * _u[_qp] * _u[_qp] + (_K_L - 2 * _K_D) * _u[_qp] + 1))) / ((_K_D * _K_D - _K_D * _K_L * _u[_qp] * _u[_qp] + (_K_L + 2 * _K_D) * _u[_qp] + 1) * (_K_D * _K_D - _K_D * _K_L * _u[_qp] * _u[_qp] + (_K_L + 2 * _K_D) * _u[_qp] + 1) * (_K_D * _K_D - _K_D * _K_L * _u[_qp] * _u[_qp] + (_K_L + 2 * _K_D) * _u[_qp] + 1) * (_K_D * _K_D - _K_D * _K_L * _u[_qp] * _u[_qp] + (_K_L + 2 * _K_D) * _u[_qp] + 1));
}

Real
FluxBCudot2::computeQpResidual()
{
	if ( theta() <= 1)
	{
		return  _test[_i][_qp]  *  gamma() + ((1 - theta()) * ( - _couple_var[_qp] * _kf_1 * _func_1.value(_t, _q_point[_qp]) + _u[_qp] * _kb_1 *  _func_2.value(_t, _q_point[_qp])) + theta() * (- _couple_var[_qp] * _kf_2 * _func_3.value(_t, _q_point[_qp]) + _u[_qp] * _kb_2 *  _func_4.value(_t, _q_point[_qp])));
	}
	else
	{
		return  _test[_i][_qp]  * gamma() - (_couple_var[_qp] * _kf_2 * _func_3.value(_t, _q_point[_qp])) + (_u[_qp] * _kb_2 *  _func_4.value(_t, _q_point[_qp]));
	}
}










