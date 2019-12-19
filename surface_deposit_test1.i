inactive = 'Bounds'
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  xmax = 1
  ny = 5
[]

[Variables]
  [C_O]
    order = FIRST
    family = LAGRANGE
    scaling = 1
    [InitialCondition]
      type = FunctionIC
      function = C_O_IC_function
      variable = C_O
    []
  []
  [C_R]
    family = LAGRANGE
    order = FIRST
    [InitialCondition]
      type = ConstantIC
      variable = C_R
      value = 0
    []
  []
  [C_O_1]
    family = LAGRANGE
    scaling = 1E5
    order = FIRST
    [InitialCondition]
      type = ConstantIC
      value = 0
      variable = C_O_1
    []
  []
[]

[Kernels]
  [C_O_diff]
    type = MatDiffusion
    variable = C_O
    D_name = D_O
  []
  [C_O_dot]
    type = TimeDerivative
    variable = C_O
  []
  [C_R]
    type = NullKernel
    variable = C_R
  []
  [C_O_1_diff]
    type = MatDiffusion
    D_name = D_O_1
    variable = C_O_1
  []
  [C_O_1_dot]
    type = TimeDerivative
    variable = C_O_1
  []
[]

[BCs]
  [C_O_right]
    type = DirichletBC
    variable = C_O
    boundary = 'right'
    value = 1
  []
  [C_O_left_coupled_Flux]
    # C_O (primary var) will couple the value of C_R (coupled var)
    type = FluxBCudot
    variable = C_O
    boundary = 'left'
    coupled_var = 'C_O_1'
    Exp_1 = Exp_1
    Exp_2 = Exp_2
    sigma = 0.001
  []
  [C_O_1_right]
    type = DirichletBC
    variable = C_O_1
    boundary = 'right'
    value = 0
  []
  [C_O_1_left_coupled_Flux]
    # C_O (primary var) will couple the value of C_R (coupled var)
    type = FluxBCudot1
    variable = C_O_1
    boundary = 'left'
    coupled_var1 = 'C_O'
    Exp_4 = Exp_4
    Exp_1 = Exp_1
    Exp_3 = Exp_3
    Exp_2 = Exp_2
    coupled_var = 'C_R'
    sigma1 = 0.001
  []
[]

[Materials]
  [Diffusivity_of_C_O]
    type = GenericConstantMaterial
    prop_names = 'D_O'
    prop_values = '0.01'
  []
  [Diffusivity_of_C_O_1]
    type = GenericConstantMaterial
    prop_values = '0.01'
    prop_names = 'D_O_1'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  # solve_type = NEWTON
  # petsc_options_iname = '-snes_type'
  # petsc_options_value = 'test'
  type = Transient
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  petsc_options = '-ksp_converged_reason -snes_converged_reason -snes_test_display'
  solve_type = PJFNK
  num_steps = 1000
  end_time = 4
  dtmax = 0.5e-2
  line_search = basic
  dt = 1e-5
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.5e-10
  []
[]

[Outputs]
  # execute_on = 'TIMESTEP_END'
  # csv = true
  exodus = true
[]

[Debug]
  show_var_residual_norms = true
[]

[Postprocessors]
  [C_O]
    type = NodalVariableValue
    nodeid = 0
    variable = C_O
  []
  [Flux_C_O]
    type = SideFluxIntegral
    diffusivity = D_O
    variable = 'C_O'
    boundary = 'left'
  []
  [C_R]
    type = NodalVariableValue
    nodeid = 0
    variable = C_R
  []
  [E]
    type = FunctionValuePostprocessor
    function = E
  []
  [Flux_C_O_1]
    type = SideFluxIntegral
    diffusivity = D_O_1
    variable = 'C_O_1'
    boundary = 'left'
  []
  [C_O_1]
    type = NodalVariableValue
    nodeid = 0
    variable = C_O_1
  []
  [i_1]
    type = FunctionValuePostprocessor
    function = A_1
  []
  [i_2]
    type = FunctionValuePostprocessor
    function = A_2
  []
  [overall_i]
    type = FunctionValuePostprocessor
    function = oevrall
  []
[]

[NodalKernels]
  [Udeposit_dot]
    type = TimeDerivativeNodalKernel
    variable = C_R
    boundary = 'left'
  []
  [Udeposit_other]
    type = Udeposit
    variable = C_R
    coupled_var = 'C_O_1'
    boundary = 'left'
    Exp_3 = Exp_3
    Exp_4 = Exp_4
  []
[]

[Functions]
  Ef1 = '1.5'
  Ef2 = '0.5'
  E1 = '2'
  v = '1'
  [Exp_1]
    type = ParsedFunction
    vars = 'n F R T alpha'
    value = 'if(t<=2, exp(-alpha*n*F*(${E1}-${v}*t-${Ef1})/R/T), exp(-alpha*n*F*(${E1}+${v}*t-${Ef1}-2*${v}*2)/R/T))'
    vals = '1 96485 8.314 300 0.5'
  []
  [Exp_2]
    type = ParsedFunction
    vars = 'n F R T alpha'
    value = 'if(t<=2, exp((1-alpha)*n*F*(${E1}-${v}*t-${Ef1})/R/T), exp((1-alpha)*n*F*(${E1}+${v}*t-${Ef1}-2*${v}*2)/R/T))'
    vals = '1 96485 8.314 300 0.5'
  []
  [E]
    type = ParsedFunction
    value = 'if(t<=2, ${E1}-${v}*t, ${E1}+${v}*t-2*${v}*2)'
  []
  [C_O_IC_function]
    type = ParsedFunction
    vars = 'a b'
    value = 'a*x+b'
    vals = '1E-9 1'
  []
  [Exp_3]
    type = ParsedFunction
    vars = 'n F R T alpha'
    value = 'if(t<=2, exp(-alpha*n*F*(${E1}-${v}*t-${Ef2})/R/T), exp(-alpha*n*F*(${E1}+${v}*t-${Ef2}-2*${v}*2)/R/T))'
    vals = '2 96485 8.314 300 0.5'
  []
  [Exp_4]
    type = ParsedFunction
    vars = 'n F R T alpha'
    value = 'if(t<=2, exp((1-alpha)*n*F*(${E1}-${v}*t-${Ef2})/R/T), exp((1-alpha)*n*F*(${E1}+${v}*t-${Ef2}-2*${v}*2)/R/T))'
    vals = '2 96485 8.314 300 0.5'
  []
  [A_1]
    type = ParsedFunction
    vars = 'Exp_1 C_O Exp_2 C_O_1'
    value = 'Exp_1*C_O-C_O_1*Exp_2'
    vals = 'Exp_1 C_O Exp_2 C_O_1'
  []
  [A_2]
    type = ParsedFunction
    vars = 'Exp_3 C_O_1 Exp_4 C_R'
    value = 'Exp_3*C_O_1-C_R*Exp_4'
    vals = 'Exp_3 C_O_1 Exp_4 C_R'
  []
  [oevrall]
    type = ParsedFunction
    vars = 'A_1 A_2'
    value = '-96485*A_1-2*96485*A_2'
    vals = 'A_1 A_2'
  []
[]

[Bounds]
[]
