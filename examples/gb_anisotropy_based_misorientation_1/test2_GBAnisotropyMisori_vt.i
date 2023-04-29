# 本脚本用于对自定义 GBAnisotropyMisori 模型的测试
# 相比于原有的GBAnisotropy, GBAnisotropyMisori是基于GrainTracker构建
# 功能1: 基于欧拉角所计算的取向差来实时GB energe and mobiity；
# 功能2: 在每个正交点处，sigam_ij and mob_ij都不相同；
# 功能3：GBAnisotropy中考虑的孪晶的低能和低迁移率特性；


my_filename = "gbAnisotropyMisori_vt3_1"

my_connecting_threshold = 1.0e-6
my_mob_matrix = 0.0
my_interval = 5

[Materials]
  [./CuGrGranisotropic]
    type = GBAnisotropyMisori
    T = 973.15 # K
    wGB = 2.0

    sigma_matrix = 0.708
    GBsigma_HAGB = 0.708
    GBmob_HAGB = 2.4e-12
    mob_matrix = ${my_mob_matrix}

    grain_tracker = grain_tracker
    euler_angle_provider = euler_angle_file

    gb_energy_anisotropy = false
    gb_mobility_anisotropy = true

    output_properties = 'kappa_op gamma_asymm L mu misori_angle twinning_type mob_00' #  
    outputs = my_exodus
  [../]
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 50
  ny = 50
  nz = 0
  xmax = 200
  ymax = 200
  zmax = 0
  elem_type = QUAD4
[]

[GlobalParams]
  op_num = 5
  var_name_base = gr

  length_scale = 1.0e-6
  time_scale = 1.0
[]

[Variables]
  [./PolycrystalVariables]
  [../]
[]

[UserObjects]
  [./euler_angle_file]
    type = EulerAngleFileReader
    file_name = grn_5_test2_2D.tex
  [../]
  [./voronoi]
    type = PolycrystalVoronoi
    grain_num = 5 # Number of grains
    coloring_algorithm = bt
    rand_seed = 50
  [../]
  [./grain_tracker]
    type = GrainTracker
    threshold = 0.5
    connecting_threshold = ${my_connecting_threshold}
    halo_level = 3
    flood_entity_type = ELEMENTAL
    compute_var_to_feature_map = true
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  [../]
  [./term]
    type = Terminator
    expression = 'gr1_area < 1e3'
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
    [../]
  [../]
[]

[BCs]
  [./Periodic]
    [./top_bottom]
      auto_direction = 'x y' # 
    [../]
  [../]
[]

[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
  [./unique_grains]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./var_indices]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./euler_angle]
    order = CONSTANT
    family = MONOMIAL
  [../]  
  [mob_00]
    order = CONSTANT
    family = MONOMIAL
  []
  [mob_01]
    order = CONSTANT
    family = MONOMIAL
  []
  [mob_02]
    order = CONSTANT
    family = MONOMIAL
  []
  [mob_03]
    order = CONSTANT
    family = MONOMIAL
  []
  [mob_04]
    order = CONSTANT
    family = MONOMIAL
  []
  [mob_10]
    order = CONSTANT
    family = MONOMIAL
  []
  [mob_11]
    order = CONSTANT
    family = MONOMIAL
  []
  [mob_12]
    order = CONSTANT
    family = MONOMIAL
  []
  [mob_13]
    order = CONSTANT
    family = MONOMIAL
  []
  [mob_14]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Kernels]
  [./PolycrystalKernel]
  [../]
[]

[AuxKernels]
  [./BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = 'initial timestep_end'
  [../]
  [./unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    execute_on = timestep_end
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
  [../]
  [./var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    execute_on = timestep_end
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
  [../]
  [./euler_angle]
    type = OutputEulerAngles
    variable = euler_angle
    euler_angle_provider = euler_angle_file
    grain_tracker = grain_tracker
    output_euler_angle = 'phi1'
  [../]
  [mob_00]
    type = MaterialStdVectorAux
    property = mob_ij
    index = 0
    variable = mob_00
  [] 
  [mob_01]
    type = MaterialStdVectorAux
    property = mob_ij
    index = 1
    variable = mob_01
  []   
  [mob_02]
    type = MaterialStdVectorAux
    property = mob_ij
    index = 2
    variable = mob_02
  []   
  [mob_03]
    type = MaterialStdVectorAux
    property = mob_ij
    index = 3
    variable = mob_03
  []   
  [mob_04]
    type = MaterialStdVectorAux
    property = mob_ij
    index = 4
    variable = mob_04
  []  
  [mob_10]
    type = MaterialStdVectorAux
    property = mob_ij
    index = 5
    variable = mob_10
  [] 
  [mob_11]
    type = MaterialStdVectorAux
    property = mob_ij
    index = 6
    variable = mob_11
  []   
  [mob_12]
    type = MaterialStdVectorAux
    property = mob_ij
    index = 7
    variable = mob_12
  []   
  [mob_13]
    type = MaterialStdVectorAux
    property = mob_ij
    index = 8
    variable = mob_13
  []   
  [mob_14]
    type = MaterialStdVectorAux
    property = mob_ij
    index = 9
    variable = mob_14
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
  [./bnd_length]
    type = GrainBoundaryArea
  [../]
  [./gr1_area]
    type = ElementIntegralVariablePostprocessor
    variable = gr1
  [../] 
[]

# [VectorPostprocessors]
#   [./grain_volumes] 
#     type = FeatureDataVectorPostprocessor
#     flood_counter = grain_tracker
#     # output_centroids = true
#   [../]
# []

[Executioner]
  type = Transient
  scheme = bdf2 # Type of time integration (2nd order backward euler), defaults to 1st order backward euler

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  # Uses newton iteration to solve the problem.
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -mat_mffd_type'
  petsc_options_value = 'hypre boomeramg 101 ds'

  l_max_its = 30 # Max number of linear iterations
  l_tol = 1e-4 # Relative tolerance for linear solves
  nl_max_its = 40 # Max number of nonlinear iterations
  nl_rel_tol = 1e-10 # Absolute tolerance for nonlienar solves

  start_time = 0.0
  end_time = 7e3
  # num_steps = 6

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]
  [./Adaptivity]
    initial_adaptivity = 3
    refine_fraction = 0.8
    coarsen_fraction = 0.05
    max_h_level = 3
  [../]
[]

[Outputs]
  [my_exodus]
    file_base = ./ex_${my_filename}/out_${my_filename} 
    interval = ${my_interval}
    type = Nemesis
  [../]
  [./csv]
    file_base = ./csv_${my_filename}/out_${my_filename}
    # interval = 5
    type = CSV
  [../]
  print_linear_residuals = false
[]