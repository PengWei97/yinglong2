my_filename = "y_test1"

[Mesh]
  type = GeneratedMesh
  dim = 2 # Problem dimension
  nx = 20 # Number of elements in the x direction
  ny = 20 # Number of elements in the y direction
  xmax = 500 # Maximum x-coordinate of mesh
  xmin = 0 # Minimum x-coordinate of mesh
  ymax = 500 # Maximum y-coordinate of mesh
  ymin = 0 # Minimum y-coordinate of mesh
  elem_type = QUAD4 # Type of elements used in the mesh
[]

[GlobalParams]
  # Parameters used by several kernels that are defined globally to simplify input file
  op_num = 3 # Number of order parameters used
  var_name_base = gr # base name of grains
  wGB = 20
  v = 'gr0 gr1 gr2' # Names of the grains
  theta1 = 90 # Angle the first grain makes at the triple junction
  theta2 = 90 # Angle the second grain makes at the triple junction
  length_scale = 1.0e-9 # Length scale in nm
  time_scale = 1.0e-9 # Time scale in ns
[]

# [UserObjects]
#   [grain_tracker]
#     type = GrainTracker
#     # to reduce the number of order parameters needed to model a large polycrystal system.
#     # halo_level = 4
#     connecting_threshold = 0.01 # 0.001
#     flood_entity_type = ELEMENTAL
#     execute_on = 'initial timestep_begin' #  TIMESTEP_BEGIN TIMESTEP_END
#     compute_halo_maps = true # For displaying HALO fields
#     # polycrystal_ic_uo = ebsd
#     compute_var_to_feature_map = true
#   []
# []

[Variables]
  [./gr0]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
       type = TricrystalTripleJunctionIC
       op_index = 1
    [../]
  [../]

  [./gr1]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
       type = TricrystalTripleJunctionIC
       op_index = 2
    [../]
  [../]

  [./gr2]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
       type = TricrystalTripleJunctionIC
       op_index = 3
    [../]
  [../]
[]

# [BCs]
#   [./Periodic]
#     [./All]
#       auto_direction = 'x'
#       variable = 'gr0 gr1 gr2'
#     [../]
#   [../]
# []

[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
  # [./unique_grains]
  #   order = FIRST
  #   family = LAGRANGE
  # [../]
  # [./var_indices]
  #   order = FIRST
  #   family = LAGRANGE
  # [../]
[]

[Kernels]
  # Kernels block where the kernels defining the residual equations are set up
  [./PolycrystalKernel]
    # Custom action creating all necessary kernels for grain growth.  All input parameters are up in GlobalParams
  [../]
[]

[AuxKernels]
  # AuxKernel block, defining the equations used to calculate the auxvars
  [./bnds_aux]
    # AuxKernel that calculates the GB term
    type = BndsCalcAux
    variable = bnds
    execute_on = 'initial timestep_end'
  [../]
[] 

[Materials]
  [./CuGrGranisotropic]
    type = GBAnisotropy
    T = 450 # K


    # molar_volume_value = 7.11e-6 #Units:m^3/mol
    Anisotropic_GB_file_name = anisotropy_mobility.txt
    inclination_anisotropy = false

    output_properties = 'kappa_op L mu gamma_asymm'
    outputs = my_exodus
  [../]
[]

[Postprocessors]
  [./dt]
    # Outputs the current time step
    type = TimestepSize
  [../]

  [./gr1_area]
    type = ElementIntegralVariablePostprocessor
    variable = gr1
  [../]
  [./gr2_area]
    type = ElementIntegralVariablePostprocessor
    variable = gr2
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 31'

  l_max_its = 30
  l_tol = 1e-4
  nl_max_its = 40
  nl_rel_tol = 1e-9
  start_time = 0.0
  num_steps =  20

  [TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.8
    dt = 0.1
    growth_factor = 1.1
    optimal_iterations = 7
  []
  [Adaptivity]
    initial_adaptivity = 4
    refine_fraction = 0.8
    coarsen_fraction = 0.3 #0.05
    max_h_level = 3
  []
[]

[Outputs]
  file_base = ./${my_filename}/out_${my_filename}
  [my_exodus]
    type = Exodus # Nemesis Exodus
    interval = 2
  [../]
  print_linear_residuals = false
  csv = true
[]

[Problem]
  solve = false
[]