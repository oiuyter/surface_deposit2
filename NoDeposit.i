[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 10
  xmax = 1
  ny = 10
[]

[Variables]
  [C_O]
    order = FIRST
    family = LAGRANGE
    [InitialCondition]
      type = ConstantIC
      value = 1
      variable = C_O
    []
  []
  [C_R]
    order = FIRST
    family = LAGRANGE
    scaling = 1E1
    [InitialCondition]
      type = FunctionIC
      function = C_R_IC_funtion
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
  [C_R_diff]
    type = MatDiffusion
    variable = C_R
    D_name = D_R
  []
  [C_O_dot]
    type = TimeDerivative
    variable = C_O
  []
  [C_R_dot]
    type = TimeDerivative
    variable = C_R
  []
[]

[BCs]
  [C_R_right]
    type = DirichletBC
    variable = C_R
    boundary = 'right'
    value = 0
  []
  [C_O_right]
    type = DirichletBC
    variable = C_O
    boundary = 'right'
    value = 1
  []
  [C_R_left_coupled_Flux]
    # C_R (primary var) taking the flux from C_O (coupled)
    type = FluxBCudot
    variable = C_R
    boundary = 'left'
    coupled_var = 'C_O'
  []
  [C_O_left_coupled_Flux]
    # C_O (primary var) will couple the value of C_R (coupled var)
    type = FluxBCudot
    variable = C_O
    boundary = 'left'
    coupled_var = 'C_R'
  []
[]

[Materials]
  [Diffusivity_of_C_O]
    type = GenericConstantMaterial
    prop_names = 'D_O'
    prop_values = '0.01'
  []
  [Diffusivity_of_C_R]
    type = GenericConstantMaterial
    prop_names = 'D_R'
    prop_values = '0.01'
  []
  [K_O]
    type = GenericConstantMaterial
    prop_values = '0.1'
    prop_names = 'K_O'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  # petsc_options_iname = '-pc_type -pc_hypre_type'
  # petsc_options_value = 'hypre boomeramg'
  # solve_type = NEWTON
  # petsc_options_iname = '-snes_type'
  # petsc_options_value = 'test'
  type = Transient
  petsc_options = '-ksp_converged_reason -snes_converged_reason -snes_test_display'
  solve_type = PJFNK
  num_steps = 1000
  end_time = 1.8
  dtmax = 0.5e-2
  line_search = basic
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.5e-2
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

[Functions]
  [C_R_IC_funtion]
    # Give C_R a small value to trigger the flux exchange
    type = ParsedFunction
    vars = 'a'
    value = 'a*x-a'
    vals = '-1E-9'
  []
[]

[Postprocessors]
  [C_O]
    type = NodalVariableValue
    nodeid = 0
    variable = C_O
  []
  [C_R]
    type = NodalVariableValue
    nodeid = 0
    variable = C_R
  []
  [Flux_C_O]
    type = SideFluxIntegral
    diffusivity = D_O
    variable = 'C_O'
    boundary = 'left'
  []
  [Flux_C_R]
    type = SideFluxIntegral
    diffusivity = D_R
    variable = 'C_R'
    boundary = 'left'
  []
[]
