# This test is used for grain growth considering stored energy driving 
# (refer to DeformedGrain.i), and the initial miscrosturcture and GNDs 
# (geometrically necessary dislocation densities) from EBSD data. 
# add 考虑晶界能和晶界迁移率各项异性
# TODO:使用MisorientationAngleCalculator,目前识别孪晶有问题

my_filename = "case1_gbAnisotropy"
# my_filename1 = "case1_gbAnisotropy_v1"

[Mesh]
  [ebsd_mesh]
    type = EBSDMeshGenerator
    # EBSD Data with GNDs during isothermal annealing with GNDs at 700℃ 
    filename = local_ebsd_Ti_Rho.inl
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
  [grain_tracker]
    type = GrainTrackerMerge
    threshold = 0.5
    connecting_threshold = 0.04
    halo_level = 3
    flood_entity_type = ELEMENTAL
    polycrystal_ic_uo = ebsd
    compute_var_to_feature_map = true

    merge_grains_based_misorientaion = true
    euler_angle_provider = ebsd_reader

    remap_grains = true
  []
  [./term]
    type = Terminator
    expression = 'grain_tracker < 10'
  [../]
[]

[ICs]
  [PolycrystalICs]
    [PolycrystalColoringIC]
      polycrystal_ic_uo = ebsd
    []
  []
[]

[Variables]
  [PolycrystalVariables]
  []
[]

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
[]

[Kernels]
  [./PolycrystalKernel]
  [../]
  [./PolycrystalStoredEnergyEBSD]
    
    GNDs_provider = ebsd_reader
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

[Materials]
  [./CuGrGr]
    # Material properties
    type = GBAnisotropyMisAng2
    T = 973.15 # K
  
    matrix_sigma = 0.9
    matrix_mob = 2.5e-12
    matrix_Q = 0.23

    inclination_anisotropy = false
    misorientation_anisotropy = true
    gb_mobility_anisotropy = true

    delta_theta_HAGB = 15.0
    GBsigma_HAGB = 0.9
    GBmob_HAGB = 2.5e-12
    TT1_sigma = 0.1019
    CT1_sigma = 0.0616
    wGB = 1.0

    euler_angle_provider = ebsd_reader

    output_properties = 'kappa_op gamma_asym L mu misori_angle twinning_type'
    outputs = my_exodus
  [../]
  [./deformed]
    type = DeformedGrainEBSDMaterial
    
    GNDs_provider = ebsd_reader
    output_properties = 'rho_eff'
    outputs = my_exodus
  [../]
[]

[VectorPostprocessors]
  [./grain_volumes] 
    type = FeatureDataVectorPostprocessor
    flood_counter = grain_tracker
    output_centroids = true
  [../]
[]

[Postprocessors]
  [dt]
    type = TimestepSize
  []
  [DOFs]
    type = NumDOFs
  []
  [./run_time]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  [../]
  [./ngrains]
    type = FeatureFloodCount
    variable = bnds
    threshold = 0.7
  [../]
  [./bnd_length]
    type = GrainBoundaryArea
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = PJFNK

  petsc_options_iname = '-pc_type -pc_hypre_type -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value = 'hypre    boomeramg      0.7'

  l_tol = 1.0e-4
  l_max_its = 20
  nl_max_its = 15
  nl_rel_tol = 1.0e-8
  dtmin = 1.0e-4

  start_time = 0.0
  num_steps = 1

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
    coarsen_fraction = 0.3
    max_h_level = 3
  []
[]

[Outputs]
  file_base = ./${my_filename}/out_${my_filename}
  # [./my_checkpoint]
  #   # interval = 5
  #   type = Checkpoint
  #   additional_execute_on = 'FINAL'
  # [../]
  [my_exodus]
    # file_base = ./${my_filename1}/out_${my_filename}    
    # interval = 5
    type = Nemesis
  [../]
  # [./csv]
  #   file_base = ./csv_${my_filename1}/out_${my_filename}
  #   # interval = 5
  #   type = CSV
  # [../]
  print_linear_residuals = false
[]