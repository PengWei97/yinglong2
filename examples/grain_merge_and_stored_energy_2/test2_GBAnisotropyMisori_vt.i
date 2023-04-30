# 本脚本用于对自定义 GBAnisotropyMisori 模型的测试
# 相比于原有的GBAnisotropy, GBAnisotropyMisori是基于GrainTracker构建
# 功能1: 基于欧拉角所计算的取向差来实时GB energe and mobiity；
# 功能2: 在每个正交点处，sigam_ij and mob_ij都不相同；
# 功能3：GBAnisotropy中考虑的孪晶的低能和低迁移率特性；

# test
# vt1: matrix * 0.01
# vt2: matrix * 0.1

# vt5: matrix * 0.1
# vt6: matrix * 0.5
# vt7: matrix * 2.0

# vt9: matrix * 0.5
# vt10: matrix * 0.1
# vt11: matrix * 1.5

my_filename = "gbAnisotropyMisori_vt11_1"

my_tt1_mob = 3.0e-13
my_ct1_mob = 6.4e-13

[Materials]
  [./CuGrGranisotropic]
    type = GBAnisotropyMisori
    T = 973.15 # K
    wGB = 1.0
    # scale_factor_matrix = 0.1
  
    GBsigma_HAGB = 0.9565
    GBmob_HAGB = 2.4e-12

    TT1_sigma = 0.1019
    CT1_sigma = 0.0616
    TT1_mob = ${my_tt1_mob}
    CT1_mob = ${my_ct1_mob}

    euler_angle_provider = euler_angle_file

    gb_energy_anisotropy = true
    gb_mobility_anisotropy = true

    output_properties = 'kappa_op L mu gamma_asymm misori_angle twinning_type'
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

  grain_tracker = grain_tracker
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
    connecting_threshold = 1.0e-3
    halo_level = 3
    flood_entity_type = ELEMENTAL
    compute_var_to_feature_map = true
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  [../]
  [./term]
    type = Terminator
    expression = 'gr1_area < 2.5e3'
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
    variable = gr1
  [../] 
[]

[VectorPostprocessors]
  [./grain_volumes] 
    type = FeatureDataVectorPostprocessor
    flood_counter = grain_tracker
    # output_centroids = true
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2 # Type of time integration (2nd order backward euler), defaults to 1st order backward euler

  #Preconditioned JFNK (default)
  solve_type = PJFNK

  # Uses newton iteration to solve the problem.
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value = 'hypre boomeramg 31 0.7'

  l_tol = 1e-4 # Relative tolerance for linear solves
  nl_rel_tol = 1e-10 # Absolute tolerance for nonlienar solves
  l_max_its = 10 # Max number of linear iterations
  nl_max_its = 8 # Max number of nonlinear iterations
  dtmin = 1.0e-4

  start_time = 0.0
  end_time = 15000
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
  [./my_checkpoint]
    file_base = ./${my_filename}/out_${my_filename}
    interval = 5
    type = Checkpoint
    additional_execute_on = 'FINAL'
  [../] 
  [my_exodus]
    file_base = ./ex_${my_filename}/out_${my_filename} 
    # interval = 5
    type = Nemesis
  [../]
  [./csv]
    file_base = ./csv_${my_filename}/out_${my_filename}
    # interval = 5
    type = CSV
  [../]
  print_linear_residuals = false
[]