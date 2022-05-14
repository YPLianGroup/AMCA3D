simulations:
  - name: test
    time_integrator: ti_1

solvers:
   - cellular_automata
   - finite_element_method

transfers:
   - name: fe_ca
     bin_size: [5, 5, 5]   
     tolerance: 1.0e-16
     max_iterations: 20
     void_temperature: 293.0001

realms:

  - name: realm0
    type: finite_element
    dimension: 3
    mesh: ./PBF_x.txt

    solution_options:
       name: my_options

       options:
         - load_data_from_file:
             theta: theta
             for_whole_time: yes
             length_scale: 10.0
             time_scale: 1.0
             lines_for_title: 4
             lines_for_subtitle: 5
             x_offset: 0
             z_offset: 0

    output:
      output_data_base_name: ./Results/FEM
      output_frequency: 200000
      output_time_interval: 0.0002
      output_variables:
        - temperature

  - name: realm1
    type: cellular_automata_remelting
    dimension: 3

    domain:
      type: cubic
      original_point: [0.402,-0.134,-0.148]
      lateral_sizes: [0.448, 0.016, 0.15]

    discretization:
      cell_size: 0.001
      number_cells: [448, 16, 150]

    nucleation_rules:
      - surface:
          type: Gaussian
          site_density: 0.0
          mean: 0.5
          standard_deviation: 0.1

      - bulk:
          type: Gaussian
          site_density: 0
          mean: 10
          standard_deviation: 1.0

    problem_physics:
        type: RappazGandin
        # both
        melting_temperature: 1715.5
        solidus_temperature: 1598.5
        montecarlo_temperature: 1.0
        # ca
        a1: 0.0
        a2: 2.49e-3
        a3: 6.2e-4

    solution_options:
       name: my_options

       options:

         - random_type:
             rejection_sampling: no
         - phase_state:
             has_been_melted: no
         - cell_grain_length_cutoff_factor:
             factor_value: 4.0
         - output_microstructure_information:
             file_name: grainInfo.txt
             output_seeds: yes
         - load_microstructure_information:
             file_name: random
         - criterion_to_activate_grain:
             using_nucleation_seed: yes
         - parallel_specification:
             direction_with_all_processes: y
         - remove_melting:
             remove_melting: yes

    output:
      output_data_base_name: ./Results/Grains
      output_frequency: -1
      output_time_interval: 0.0002
      output_variables:
        - temperature
        - grain_orientation


time_integrators:
  - standard_time_integrator:
      name: ti_1
      termination_time: 100
      termination_step_count: 10
      time_step: 0.001
      time_step_factor: 0.4

      realms:
        - realm0
        - realm1
