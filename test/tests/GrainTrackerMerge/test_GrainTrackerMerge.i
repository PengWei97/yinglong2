# This simulation was used to test the material class GrainTrackerMerge.

my_filename = "case1_GrainTrackerMerge"

[Mesh]
  type = GeneratedMesh
  dim = 2
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
    file_name = grn_10_test_2D.tex
  [../]
  [./voronoi]
    type = PolycrystalVoronoi
    coloring_algorithm = jp
    rand_seed = 20
  [../]
  [./grain_tracker]
    type = GrainTrackerMerge # GrainTrackerMerge
    threshold = 0.2
    connecting_threshold = 0.08
    compute_halo_maps = true # Only necessary for displaying HALOS
    compute_var_to_feature_map = true

    merge_grains_basedMisorAngle = true
    euler_angle_provider = euler_angle_file
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
  [./PolycrystalKernel]
  [../]
[]

[AuxKernels]
  [./bnds_aux]
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
    type = OutputEulerAngles
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
  [./CuGrGr]
    # Material properties
    type = GBEvolution
    T = 450 # Constant temperature of the simulation (for mobility calculation)
    wGB = 14 # Width of the diffuse GB
    GBmob0 = 2.5e-6 #m^4(Js) for copper from Schoenfelder1997
    Q = 0.23 #eV for copper from Schoenfelder1997
    GBenergy = 0.708 #J/m^2 from Schoenfelder1997
  [../]
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
  num_steps = 10

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
    type = Exodus
  [../]
  print_linear_residuals = false
  csv = true
  [./console]
    type = Console
    max_rows = 20 # Will print the 20 most recent postprocessor values to the screen
  [../]
[]

# # Here we'll load the wrong library and check for the correct error condition
# [Problem]
#   register_objects_from = 'TensorMechanicsApp'
#   library_path = '../../../../../tensor_mechanics/lib'
# []

[33mGrain #22 (op: 0) and Grain #28 (op: 5) was merged (misor angle: 0.454124).
[39m[33mGrain #101 (op: 1) and Grain #108 (op: 0) was merged (misor angle: 0.428499).
[39m[33mGrain #226 (op: 1) and Grain #234 (op: 2) was merged (misor angle: 0.920237).
[39m[33mGrain #68 (op: 2) and Grain #71 (op: 6) was merged (misor angle: 0.681307).
[39m[33mGrain #61 (op: 3) and Grain #62 (op: 4) was merged (misor angle: 0).
[39m[33mGrain #17 (op: 4) and Grain #41 (op: 2) was merged (misor angle: 0.601591).
[39m[33mSplit Grain (#22) detected on unmatched OPs (0, 5) attempting to remap to 0.
[39m[33mSplit Grain (#101) detected on unmatched OPs (0, 1) attempting to remap to 0.
[39m[33mSplit Grain (#226) detected on unmatched OPs (1, 2) attempting to remap to 1.
[39m[33mSplit Grain (#17) detected on unmatched OPs (2, 4) attempting to remap to 2.
[39m[33mSplit Grain (#68) detected on unmatched OPs (2, 6) attempting to remap to 2.
[39m[33mSplit Grain (#61) detected on unmatched OPs (3, 4) attempting to remap to 3.
[39m[33mGrain #17 intersects Grain #63 (variable index: 4)
[39m[32m- Depth 0: Remapping grain #17 from variable index 4 to 9 whose closest grain (#176) is at a distance of 63.6905

