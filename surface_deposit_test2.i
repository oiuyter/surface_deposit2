[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 100
  ny = 10
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
      value = 0
      variable = C_R
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
  [C_R_diff]
    type = MatDiffusion
    variable = C_R
    D_name = D_R
  []
  [C_R_dot]
    type = TimeDerivative
    variable = C_R
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
    coupled_var = 'C_R'
    Exp_1 = Exp_1
    Exp_2 = Exp_2
    Exp_3 = Exp_3
    Exp_4 = Exp_4
    K_L = 855
    K_D = 918.5
    G_S = .00000001
    kf_1 = 32.965
    kb_1 = 32.965
    kf_2 = 31.05
    kb_2 = 31.05
    sigma = .01
  []
  [C_R_right]
    type = DirichletBC
    variable = C_R
    boundary = 'right'
    value = 0
  []
  [C_R_left_coupled_Flux]
    # C_R (primary var) will couple the value of C_0 (coupled var)
    type = FluxBCudot2
    variable = C_R
    boundary = 'left'
    coupled_var = 'C_O'
    Exp_1 = Exp_1
    Exp_2 = Exp_2
    Exp_3 = Exp_3
    Exp_4 = Exp_4
    K_L = 855
    K_D = 918.5
    G_S = .00000001
    kf_1 = 32.965
    kb_1 = 32.965
    kf_2 = 31.05
    kb_2 = 31.05
    sigma = .01
  []
[]

[Materials]
  [Diffusivity_of_C_O]
    type = GenericConstantMaterial
    prop_names = 'D_O'
    prop_values = '.00000602'
  []
  [Diffusivity_of_C_R]
    type = GenericConstantMaterial
    prop_names = 'D_R'
    prop_values = '.000000000000000000001'
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
  num_steps = 2000
  end_time = 1.8
  dtmax = 0.2e-2
  line_search = basic
  dt = 1
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.3e-2
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
  [Flux_C_R]
    type = SideFluxIntegral
    diffusivity = D_R
    variable = 'C_R'
    boundary = 'left'
  []
  [E]
    type = FunctionValuePostprocessor
    function = E
  []
  [Exp_1]
    type = FunctionValuePostprocessor
    function = Exp_1
  []
  [Exp_2]
    type = FunctionValuePostprocessor
    function = Exp_2
  []
  [Exp_3]
    type = FunctionValuePostprocessor
    function = Exp_3
  []
  [Exp_4]
    type = FunctionValuePostprocessor
    function = Exp_4
  []
[]

[NodalKernels]
  [Udeposit_dot]
    type = TimeDerivativeNodalKernel
    variable = C_R
    boundary = 'left'
  []
[]

[Functions]
  Ef1 = '2.0755'
  Ef2 = '2.0755'
  E1 = '.5'
  v = '.555555556'
  end_time = '1.8'
  [Exp_1]
    type = ParsedFunction
    vars = 'n F R T alpha1'
    value = 'if(t<=${end_time}/2, exp(-alpha1*n*F*(${E1}-${v}*t-${Ef1})/R/T), exp(-alpha1*n*F*(${E1}+${v}*t-${Ef1}-${v}*${end_time})/R/T))'
    vals = '1 96485 8.314 723 0.685'
  []
  [Exp_2]
    type = ParsedFunction
    vars = 'n F R T alpha1'
    value = 'if(t<=${end_time}/2, exp((1-alpha1)*n*F*(${E1}-${v}*t-${Ef1})/R/T), exp((1-alpha1)*n*F*(${E1}+${v}*t-${Ef1}-${v}*${end_time})/R/T))'
    vals = '1 96485 8.314 723 0.685'
  []
  [E]
    type = ParsedFunction
    value = 'if(t<=(${end_time}/2), ${E1}-${v}*t, ${E1}+${v}*t-2*${v}*(${end_time}/2))'
  []
  [Exp_3]
    type = ParsedFunction
    vars = 'n F R T alpha2'
    value = 'if(t<=${end_time}/2, exp(-alpha2*n*F*(${E1}-${v}*t-${Ef2})/R/T), exp(-alpha2*n*F*(${E1}+${v}*t-${Ef2}-${v}*${end_time})/R/T))'
    vals = '2 96485 8.314 723 0.685'
  []
  [Exp_4]
    type = ParsedFunction
    vars = 'n F R T alpha2'
    value = 'if(t<=${end_time}/2, exp((1-alpha2)*n*F*(${E1}-${v}*t-${Ef2})/R/T), exp((1-alpha2)*n*F*(${E1}+${v}*t-${Ef2}-${v}*${end_time})/R/T))'
    vals = '2 96485 8.314 723 0.685'
  []
  [C_O_IC_function]
    type = ParsedFunction
    vars = 'a b'
    value = 'a*x+b'
    vals = '1E-9 1'
  []
[]
