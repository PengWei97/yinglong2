# heterogeneous stored energy distribution

my_filename = "test1_gbIsotropy"

my_filename2 = "test1_gbIsotropy1"

[Mesh]
  [ebsd_mesh]
    type = EBSDMeshGenerator
    filename = local3_10min_withTwin_Ti700_rho.inl
  []
  # using Distributed Mesh
[]

[GlobalParams]
  op_num = 10
  var_name_base = gr
  length_scale = 1.0e-6
  time_scale = 1.0
[]

[UserObjects]
  [ebsd_reader]
    # Get Euler angles, coordinates, grain ID, phase ID, symmetry, GNDs from EBSD file
    type = EBSDReader
    ebsd_meshgenerator = ebsd_mesh
    custom_columns = 1
    # execute_on = 'initial'
  []
  [ebsd]
    type = PolycrystalEBSD
    # Object for setting up a polycrystal structure from an EBSD Datafile
    coloring_algorithm = jp
    ebsd_reader = ebsd_reader
    enable_var_coloring = true
    compute_var_to_feature_map = false
    # execute_on = 'initial'
    # Instruct the Postprocessor to compute the active vars to features map
  []
  [grain_tracker]
    type = GrainTrackerGG2
    threshold = 0.5
    connecting_threshold = 0.04
    halo_level = 3
    flood_entity_type = ELEMENTAL
    polycrystal_ic_uo = ebsd
    compute_var_to_feature_map = true

    execute_on = 'initial TIMESTEP_BEGIN' #  TIMESTEP_BEGIN TIMESTEP_END

    remerge_grains = true
    euler_angle_provider = ebsd_reader   
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
  # [./ebsd_rho]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
[]

[Kernels]
  [PolycrystalKernel]
  []
[]

[AuxKernels]
  [BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = 'initial timestep_end'
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
    execute_on = 'initial timestep_end'
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
  []
  [phi1]
    type = OutputEulerAnglesGG
    variable = phi1
    euler_angle_provider = ebsd_reader
    grain_tracker = grain_tracker
    output_euler_angle = 'phi1'
    execute_on = 'initial TIMESTEP_BEGIN'
  []
  [Phi]
    type = OutputEulerAnglesGG
    variable = Phi
    euler_angle_provider = ebsd_reader
    grain_tracker = grain_tracker
    output_euler_angle = 'Phi'
    execute_on = 'initial TIMESTEP_BEGIN'
  []
  [phi2]
    type = OutputEulerAnglesGG
    variable = phi2
    euler_angle_provider = ebsd_reader
    grain_tracker = grain_tracker
    output_euler_angle = 'phi2'
    execute_on = 'initial TIMESTEP_BEGIN'
  []
  # [grain_rho]
  #   type = EBSDReaderPointDataAux
  #   variable = ebsd_rho # GNDs
  #   ebsd_reader = ebsd_reader
  #   data_name = 'CUSTOM0' 
  #   execute_on = 'initial TIMESTEP_BEGIN' 
  # []
[]

[Modules]
  [PhaseField]
    [EulerAngles2RGB]
      # EulerAngle2RGBAction
      # Set up auxvariables and auxkernels to output Euler angles as RGB values interpolated across inverse pole figure
      crystal_structure = hexagonal # hexagonal cubic 
      euler_angle_provider = ebsd_reader
      grain_tracker = grain_tracker
    []
  []
[]

[Materials]
  [./CuGrGr]
    # Material properties
    type = GBEvolutionGG # Quantitative material properties for copper grain growth.  Dimensions are nm and ns
    GBEnergy = 0.96 
    GBmob0 = 2.0e-13 # Mobility prefactor
    Q = 0.23
    T = 973.15 # K
    wGB = 1.0 # mum

    grain_tracker = grain_tracker
    euler_angle_provider = ebsd_reader

    output_properties = 'kappa_op L misAngle twType'
    outputs = my_exodus
  [../]
  # [./CuGrGranisotropic]
  #   type = GBAnisotropy1MisAng5
  #   T = 973.15 # K
  #   wGB = 1.0
    
  #   grain_tracker = grain_tracker
  #   euler_angle_provider = ebsd_reader

  #   misorientation_anisotropy = false
  #   gbMobility_anisotropy = false
  #   tb_anisotropy = false
  #   inclination_anisotropy = false
   
  #   matrix_sigma = 0.96 # 0.1 0.25 0.5 0.96
  #   matrix_mob = 2.0e-13 # 0.5 1.25 2.50
  #   matrix_Q = 0.23
    
  #   delta_theta_HAB = 15.0
  #   GBsigma_HAGB = 0.96
  #   GBmob_HAGB = 2.0e-13
  #   TT1_sigma = 0.1019
  #   CT1_sigma = 0.0616

  #   # output_properties = 'kappa_op L mu gamma_asymm misAngle twType'
  #   output_properties = 'misAngle twType'
  #   outputs = my_exodus

  #   misori_angle = '0.0	6.92	7.57	8.21	9.02	9.98	11.270	12.720	13.690	14.810	16.1	17.710	19.480	21.740	23.510	24.470	25.440	28.010	29.950	30.750	32.520	34.290	35.260	38.480	41.380	42.340	43.630	44.760	48.620	50.070	51.840	52.970	53.770	57.8	61.660	62.310	64.720	65.850	67.3	72.610	74.870	75.990	76.8	78.250	79.860	81.140	82.920	84.040	85.010	85.490	90	96.920'
    
  #   sigma_energy = '0.10	0.64382	0.73515	0.75521	0.78419	0.82430	0.86665	0.88005	0.88898	0.90905	0.91131	0.91581	0.88244	0.88027	0.87586	0.86698	0.86032	0.75570	0.85820	0.92726	0.91617	0.93403	0.92069	0.81610	0.95648	0.92087	0.92313	0.92984	0.84976	0.84088	0.76075	0.72069	0.66948	0.52259	0.60954	0.64074	0.70539	0.71655	0.70545	0.68554	0.52078	0.58540	0.63887	0.68568	0.66790	0.61003	0.58557	0.57224	0.56781	0.48318	0.10	0.10'
  # [../]
  [./deformedGrain]
    type = HeterStoredEnergy
    Elas_Mod = 2.50e10
    Burg_vec = 2.95e-10

    grain_tracker = grain_tracker
    GNDs_provider = ebsd_reader
    dataGNDs_name = 'CUSTOM0'

    output_properties = 'rho_eff GNDs'
    outputs = my_exodus
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
  # end_time = 1e6
  num_steps = 10

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
  [./my_checkpoint]
    file_base = ./${my_filename}/out_${my_filename}
    type = Checkpoint
    # interval = 5
    additional_execute_on = 'FINAL' # seems to be a necessary to avoid a Checkpoint bug
  [../]
  [my_exodus]
    file_base = ./${my_filename2}/out_${my_filename}
    type = Nemesis # Nemesis Exodus
    # interval = 5
  [../]
  [./csv]
    file_base = ./csv_${my_filename}/${my_filename}
    type = CSV
    # interval = 5
  [../]
  print_linear_residuals = false
[]