my_interval = 10
my_filename = 'a1_isotropy_gg'
my_filename2 = 'a1_isotropy_gg_2'

# TODO: 小尺度晶粒长大模拟，大尺度模拟造成内存溢出，下一步需要解决

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 20
  ny = 20
  nz = 20
  xmin = 0
  xmax = 32
  ymin = 0
  ymax = 32
  zmin = 0
  zmax = 32
  elem_type = HEX8

  parallel_type = distributed # Periodic BCs distributed replicated
[]

[GlobalParams]
  op_num = 26
  var_name_base = gr
[]

[Variables]
  [./PolycrystalVariables]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[UserObjects]
  [./voronoi]
    type = PolycrystalVoronoi
    grain_num = 100 # Number of grains
    rand_seed = 200 # 301
    coloring_algorithm = jp
  [../]
  [./grain_tracker]
    type = GrainTracker
    threshold = 0.1
    compute_halo_maps = true
    compute_var_to_feature_map = true
  [../]
  [./term]
    type = Terminator
    expression = 'grain_tracker < 25'
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
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
    field_display = UNIQUE_REGION
    execute_on = 'initial timestep_end'
    flood_counter = grain_tracker
  [../]
  [./var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    field_display = VARIABLE_COLORING
    execute_on = 'initial timestep_end'
    flood_counter = grain_tracker
  [../]
[]

# [BCs]
#  [./Periodic]
#    [./All]
#      auto_direction = 'x y z'
#    [../]
#  [../]
# []

[Materials]
  [./Copper]
    type = GBEvolution
    T = 673.15 # K
    wGB = 0.5 # um
    GBmob0 = 5.0e-13 #m^4/(Js)
    Q = 0.23 #Migration energy in eV
    GBenergy = 0.708 #GB energy in J/m^2
    length_scale = 1.0e-6
    time_scale = 1.0
  [../]
[]

[Postprocessors]
  [./dt]
    type = TimestepSize
  [../]
  [./run_time]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  [../]
  [./DOFs]
    type = NumDOFs
  [../]
  [./bnd_length]
    type = GrainBoundaryArea
  [../]
  [./ngrains]
    type = FeatureFloodCount
    variable = bnds
    threshold = 0.7
  [../]
[]

[VectorPostprocessors]
  [./grain_volumes] 
    type = FeatureDataVectorPostprocessor
    flood_counter = grain_tracker
    output_centroids = true
  [../]
[]

# [Preconditioning]
#  [./SMP]
#    type = SMP
#    full = true
#  [../]
# []

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
  dtmin = 1.0e-8

  start_time = 0.0
  end_time = 5.0e3 # 5.0e3

  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.9
    dt = 0.0002
    growth_factor = 1.1
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
    interval = ${my_interval}
    type = Checkpoint
    additional_execute_on = 'FINAL'
  [../]  
  [my_exodus]
    file_base = ./ex_${my_filename2}/out_${my_filename} 
    interval = ${my_interval}
    type = Nemesis
    additional_execute_on = 'FINAL'
  [../]
  [./csv]
    file_base = ./csv_${my_filename}/out_${my_filename}
    interval = ${my_interval}
    type = CSV
  [../]
  print_linear_residuals = false
[]
