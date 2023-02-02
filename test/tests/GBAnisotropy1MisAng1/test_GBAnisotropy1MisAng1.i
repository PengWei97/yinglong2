# This simulation was used to test the material class GBAnisotropy1MisAng1.

my_filename = "case1_GBAnisotropy1MisAng1"

[Mesh]
  # Mesh block.  Meshes can be read in or automatically generated
  type = GeneratedMesh
  dim = 2 # Problem dimension
  nx = 11 # Number of elements in the x-direction
  ny = 11 # Number of elements in the y-direction
  xmin = 0    # minimum x-coordinate of the mesh
  xmax = 600 # maximum x-coordinate of the mesh
  ymin = 0    # minimum y-coordinate of the mesh
  ymax = 600 # maximum y-coordinate of the mesh
  elem_type = QUAD4  # Type of elements used in the mesh
  uniform_refine = 3 # Initial uniform refinement of the mesh

  parallel_type = replicated # Periodic BCs
[]

[GlobalParams]
  # Parameters used by several kernels that are defined globally to simplify input file
  op_num = 8 # Number of order parameters used
  var_name_base = gr # Base name of grains
  grain_num = 10 # Number of grains

  length_scale = 1.0e-9
  time_scale = 1.0e-9
[]

[Variables]
  # Variable block, where all variables in the simulation are declared
  [./PolycrystalVariables]
  [../]
[]

[UserObjects]
  [./euler_angle_file]
    type = EulerAngleFileReader
    file_name = grn_10_rand_2D.tex
  [../]
  [./voronoi]
    type = PolycrystalVoronoi
    coloring_algorithm = jp
    rand_seed = 20
  [../]
  [./grain_tracker]
    type = GrainTrackerGG
    threshold = 0.2
    connecting_threshold = 0.08
    compute_halo_maps = true # Only necessary for displaying HALOS
    compute_var_to_feature_map = true
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
  # Dependent variables
  [./bnds]
    # Variable used to visualize the grain boundaries in the simulation
  [../]
  [./unique_grains]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./var_indices]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./ghost_regions]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./halos]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./phi1]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  # Kernel block, where the kernels defining the residual equations are set up.
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
  [./unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
    execute_on = 'initial timestep_end'
  [../]
  [./var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
    execute_on = 'initial timestep_end'
  [../]
  [./ghosted_entities]
    type = FeatureFloodCountAux
    variable = ghost_regions
    flood_counter = grain_tracker
    field_display = GHOSTED_ENTITIES
    execute_on = 'initial timestep_end'
  [../]
  [./halos]
    type = FeatureFloodCountAux
    variable = halos
    flood_counter = grain_tracker
    field_display = HALOS
    execute_on = 'initial timestep_end'
  [../]
  [./phi1]
    type = OutputEulerAnglesGG
    variable = phi1
    euler_angle_provider = euler_angle_file
    grain_tracker = grain_tracker
    output_euler_angle = 'phi1'
    execute_on = 'initial timestep_end'
  [../]
[]

[BCs]
  # Boundary Condition block
  [./Periodic]
    [./top_bottom]
      auto_direction = 'x y' # Makes problem periodic in the x and y directions
    [../]
  [../]
[]

[Materials]
  [./CuGrGranisotropic]
    type = GBAnisotropy1MisAng1
    T = 600 # K
    wGB = 15
    
    grain_tracker = grain_tracker
    euler_angle_provider = euler_angle_file

    misorientation_anisotropy = true
    gbMobility_anisotropy = false
    tb_anisotropy = false
    inclination_anisotropy = false

    delta_theta_HAB = 15
    GBsigma_HAB = 0.9
    GBmob_HAB = 2.5e-6

    output_properties = 'kappa_op L mu gamma_asymm misAngle'
    outputs = my_exodus
  [../]

  # [./CuGrGr]
  #   # Material properties
  #   type = GBEvolution
  #   T = 450 # Constant temperature of the simulation (for mobility calculation)
  #   wGB = 14 # Width of the diffuse GB
  #   GBmob0 = 2.5e-6 #m^4(Js) for copper from Schoenfelder1997
  #   Q = 0.23 #eV for copper from Schoenfelder1997
  #   GBenergy = 0.708 #J/m^2 from Schoenfelder1997
  # [../]
[]

[Postprocessors]
  # Scalar postprocessors
  [./ngrains]
    type = FeatureFloodCount
    variable = bnds
    threshold = 0.7
  [../]    
  [DOFs]
    type = NumDOFs
  []
  [./dt]
    # Outputs the current time step
    type = TimestepSize
  [../]
  [./run_time]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  [../]
[]

[VectorPostprocessors]
  [./grain_volumes] 
    type = FeatureDataVectorPostprocessor
    flood_counter = grain_tracker # The FeatureFloodCount UserObject to get values from.
    execute_on = 'initial timestep_end'
    output_centroids = true
  [../]
[]


[Executioner]
  type = Transient # Type of executioner, here it is transient with an adaptive time step
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
  # end_time = 1.0e3
  num_steps = 3

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.5 # Initial time step.  In this simulation it changes.
    optimal_iterations = 6 # Time step will adapt to maintain this number of nonlinear iterations
  [../]

  [./Adaptivity]
    # Block that turns on mesh adaptivity. Note that mesh will never coarsen beyond initial mesh (before uniform refinement)
    initial_adaptivity = 2 # Number of times mesh is adapted to initial condition
    refine_fraction = 0.7 # Fraction of high error that will be refined
    coarsen_fraction = 0.1 # Fraction of low error that will coarsened
    max_h_level = 4 # Max number of refinements used, starting from initial mesh (before uniform refinement)
  [../]
[]

[Outputs]
  file_base = ./${my_filename}/out_${my_filename}
  [my_exodus]
    type = Exodus # Nemesis Exodus
  [../]
  print_linear_residuals = false
  csv = true
  [./console]
    type = Console
    max_rows = 20 # Will print the 20 most recent postprocessor values to the screen
  [../]
[]