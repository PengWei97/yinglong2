# GBAnisotropyMisori + EBSD data

# tw1: scale_factor_matrix = 1.0, my_tt1_mob = 2.4e-12
# tw2: scale_factor_matrix = 1.0, my_tt1_mob = 3.0e-13
# tw3: scale_factor_matrix = 1.0, my_tt1_mob = 1.2e-12
# tw3: scale_factor_matrix = 1.0, my_tt1_mob = 1.0e-12 + stored energy

my_filename = "gbAnisotropyMisori_tw3_1"

my_tt1_mob = 1.2e-12 #  3.0e-14
my_ct1_mob = 6.4e-13 # 6.4e-13

[Materials]
  [./CuGrGranisotropic]
    type = GBAnisotropyMisori
    T = 973.15 # K
    wGB = 1.0
  
    GBsigma_HAGB = 0.9565
    GBmob_HAGB = 2.4e-12

    TT1_sigma = 0.1019
    CT1_sigma = 0.0616
    TT1_mob = ${my_tt1_mob}
    CT1_mob = ${my_ct1_mob}

    euler_angle_provider = ebsd_reader

    gb_energy_anisotropy = true
    gb_mobility_anisotropy = true

    output_properties = 'kappa_op L mu misori_angle twinning_type'
    outputs = my_exodus
  [../]
[]

[Mesh]
  [ebsd_mesh]
    type = EBSDMeshGenerator
    # EBSD Data with GNDs during isothermal annealing with GNDs at 700℃ 
    filename = debug_loc3_tw1_10min_Ti700_rho.inl
  []
  parallel_type = distributed
[]

[GlobalParams]
  op_num = 10
  var_name_base = gr

  length_scale = 1.0e-6
  time_scale = 1.0

  grain_tracker = grain_tracker
[]

[UserObjects]
  [ebsd_reader]
    # Get Euler angles, coordinates, grain ID, phase ID, symmetry, GNDs from EBSD file
    type = EBSDReader
    custom_columns = 1
  []
  [ebsd]
    type = PolycrystalEBSD
    coloring_algorithm = jp
    ebsd_reader = ebsd_reader
    enable_var_coloring = true
    compute_var_to_feature_map = false
  []
  [./grain_tracker]
    type = GrainTrackerMerge
    threshold = 0.5
    connecting_threshold = 1.0e-2
    halo_level = 3
    flood_entity_type = ELEMENTAL
    polycrystal_ic_uo = ebsd
    compute_var_to_feature_map = true

    execute_on = 'INITIAL TIMESTEP_BEGIN'

    merge_grains_based_misorientaion = true
    euler_angle_provider = ebsd_reader
    
  [../]
  # [./term]
  #   type = Terminator
  #   expression = 'gr1_area < 500'
  # [../]
[]

[ICs]
  [./PolycrystalICs]
    [./PolycrystalColoringIC]
      polycrystal_ic_uo = ebsd
    [../]
  [../]
[]

[Variables]
  [./PolycrystalVariables]
  [../]
[]

# [Bounds]
#   [./gr0_upper_bound]
#     type = ConstantBoundsAux
#     variable = bounds_dummy
#     bounded_variable = gr0
#     bound_type = upper
#     bound_value = 1.2
#   [../]
#   [./gr0_lower_bound]
#     type = ConstantBoundsAux
#     variable = bounds_dummy
#     bounded_variable = gr0
#     bound_type = lower
#     bound_value = -0.01
#   [../]
#   [./gr1_upper_bound]
#     type = ConstantBoundsAux
#     variable = bounds_dummy
#     bounded_variable = gr1
#     bound_type = upper
#     bound_value = 1.2
#   [../]
#   [./gr1_lower_bound]
#     type = ConstantBoundsAux
#     variable = bounds_dummy
#     bounded_variable = gr1
#     bound_type = lower
#     bound_value = -0.01
#   [../]
#   [./gr2_upper_bound]
#     type = ConstantBoundsAux
#     variable = bounds_dummy
#     bounded_variable = gr2
#     bound_type = upper
#     bound_value = 1.2
#   [../]
#   [./gr2_lower_bound]
#     type = ConstantBoundsAux
#     variable = bounds_dummy
#     bounded_variable = gr2
#     bound_type = lower
#     bound_value = -0.01
#   [../]
#   [./gr3_upper_bound]
#     type = ConstantBoundsAux
#     variable = bounds_dummy
#     bounded_variable = gr3
#     bound_type = upper
#     bound_value = 1.2
#   [../]
#   [./gr3_lower_bound]
#     type = ConstantBoundsAux
#     variable = bounds_dummy
#     bounded_variable = gr3
#     bound_type = lower
#     bound_value = -0.01
#   [../]
#   [./gr4_upper_bound]
#     type = ConstantBoundsAux
#     variable = bounds_dummy
#     bounded_variable = gr4
#     bound_type = upper
#     bound_value = 1.2
#   [../]
#   [./gr4_lower_bound]
#     type = ConstantBoundsAux
#     variable = bounds_dummy
#     bounded_variable = gr4
#     bound_type = lower
#     bound_value = -0.01
#   [../]
#   [./gr5_upper_bound]
#     type = ConstantBoundsAux
#     variable = bounds_dummy
#     bounded_variable = gr5
#     bound_type = upper
#     bound_value = 1.2
#   [../]
#   [./gr5_lower_bound]
#     type = ConstantBoundsAux
#     variable = bounds_dummy
#     bounded_variable = gr5
#     bound_type = lower
#     bound_value = -0.01
#   [../]
#   [./gr6_upper_bound]
#     type = ConstantBoundsAux
#     variable = bounds_dummy
#     bounded_variable = gr6
#     bound_type = upper
#     bound_value = 1.2
#   [../]
#   [./gr6_lower_bound]
#     type = ConstantBoundsAux
#     variable = bounds_dummy
#     bounded_variable = gr6
#     bound_type = lower
#     bound_value = -0.01
#   [../]
#   [./gr7_upper_bound]
#     type = ConstantBoundsAux
#     variable = bounds_dummy
#     bounded_variable = gr7
#     bound_type = upper
#     bound_value = 1.2
#   [../]
#   [./gr7_lower_bound]
#     type = ConstantBoundsAux
#     variable = bounds_dummy
#     bounded_variable = gr7
#     bound_type = lower
#     bound_value = -0.01
#   [../]
#   [./gr8_upper_bound]
#     type = ConstantBoundsAux
#     variable = bounds_dummy
#     bounded_variable = gr8
#     bound_type = upper
#     bound_value = 1.2
#   [../]
#   [./gr8_lower_bound]
#     type = ConstantBoundsAux
#     variable = bounds_dummy
#     bounded_variable = gr8
#     bound_type = lower
#     bound_value = -0.01
#   [../]
#   [./gr9_upper_bound]
#     type = ConstantBoundsAux
#     variable = bounds_dummy
#     bounded_variable = gr9
#     bound_type = upper
#     bound_value = 1.2
#   [../]
#   [./gr9_lower_bound]
#     type = ConstantBoundsAux
#     variable = bounds_dummy
#     bounded_variable = gr9
#     bound_type = lower
#     bound_value = -0.01
#   [../]
# []

[AuxVariables]
  [bnds]
  []
  [unique_grains]
    order = CONSTANT
    family = MONOMIAL
  []
  [var_indices]
    order = CONSTANT
    family = MONOMIAL
  []
  [phi1]
    order = CONSTANT
    family = MONOMIAL
  []
  [Phi]
    order = CONSTANT
    family = MONOMIAL
  []
  [phi2]
    order = CONSTANT
    family = MONOMIAL
  []
  [./ebsd_rho]
    order = CONSTANT
    family = MONOMIAL
  [../]
  # [./bounds_dummy]
  #   order = FIRST
  #   family = LAGRANGE
  # [../]  
[]

[Kernels]
  [./PolycrystalKernel]
  [../]
[]

[AuxKernels]
  [BndsCalc]
    type = BndsCalcAux
    variable = bnds
  []
  [unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
  []
  [var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
  []
  [phi1]
    type = OutputEulerAngles
    variable = phi1
    euler_angle_provider = ebsd_reader
    
    output_euler_angle = 'phi1'
  []
  [Phi]
    type = OutputEulerAngles
    variable = Phi
    euler_angle_provider = ebsd_reader
    
    output_euler_angle = 'Phi'
  []
  [phi2]
    type = OutputEulerAngles
    variable = phi2
    euler_angle_provider = ebsd_reader
    
    output_euler_angle = 'phi2'
  []
  [grain_rho]
    type = EBSDReaderPointDataAux
    variable = ebsd_rho # GNDs
    ebsd_reader = ebsd_reader
    data_name = 'CUSTOM0'
  []
[]

[Modules]
  [PhaseField]
    [EulerAngles2RGB]
      crystal_structure = hexagonal # hexagonal cubic 
      euler_angle_provider = ebsd_reader
    []
  []
[]

[Postprocessors]
  [./dofs]
    type = NumDOFs
  [../]
  [./dt]
    type = TimestepSize
  [../]
  [./run_time]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  [../]
  # [./bnd_length]
  #   type = GrainBoundaryArea
  # [../]
  [./gr1_area]
    type = ElementIntegralVariablePostprocessor
    variable = gr5
  [../]
[]

[VectorPostprocessors]
  [./grain_volumes] 
    type = FeatureDataVectorPostprocessor
    flood_counter = grain_tracker
    output_centroids = true
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2 # Type of time integration (2nd order backward euler), defaults to 1st order backward euler

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  # Uses newton iteration to solve the problem.
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold' #  -snes_type
  petsc_options_value = 'hypre boomeramg 31 0.7' # vinewtonrsls

  l_tol = 1e-4 # Relative tolerance for linear solves
  nl_rel_tol = 1e-10 # Absolute tolerance for nonlienar solves
  l_max_its = 10 # Max number of linear iterations
  nl_max_its = 8 # Max number of nonlinear iterations
  dtmin = 1.0e-4
  # dtmax = 1.0

  start_time = 0.0
  end_time = 1e3
  # num_steps = 6

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]
  [./Adaptivity]
    initial_adaptivity = 3
    refine_fraction = 0.8
    coarsen_fraction = 0.05
    max_h_level = 2
  [../]
[]

[Outputs]
  [./my_checkpoint]
    file_base = ./${my_filename}/out_${my_filename}
    interval = 5
    type = Checkpoint
    additional_execute_on = 'FINAL'
  [../]  
  [my_exodus]
    file_base = ./ex_${my_filename}/out_${my_filename} 
    interval = 5
    type = Nemesis
  [../]
  [./csv]
    file_base = ./csv_${my_filename}/out_${my_filename}
    interval = 5
    type = CSV
  [../]
  print_linear_residuals = false
[]