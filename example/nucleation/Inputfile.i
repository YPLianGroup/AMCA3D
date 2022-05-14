simulations:
  - name: test
    time_integrator: ti_1

solvers:
   - cellular_automata

realms:
  - name: realm1
    type: cellular_automata
    dimension: 3

    domain:
      type: cubic
      original_point: [0,0,0]
      lateral_sizes: [0.03, 0.03, 0.03]

    discretization:
      cell_size: 0.0005
      number_cells: [800, 720, 252]

    nucleation_rules:
      - surface:
          type: Gaussian
          site_density: 0.0
          mean: 2
          standard_deviation: 0.5

      - bulk:
          type: Gaussian
          site_density: 10e6
          mean: 3
          standard_deviation: 1

    problem_physics:
        type: RappazGandin
        initial_temperature: -0.2
        melting_temperature: 0.0
        t_dot: -20
        a1: -0.544e-4 
        a2: 2.03e-4
        a3: 0.0

    solution_options:
       name: my_options

       options:

         - random_type:
             rejection_sampling: no
         - cell_grain_length_cutoff_factor:
             factor_value: 4.0
         - output_microstructure_information:
             file_name: grainInfo.txt
         - parallel_specification:
             direction_with_all_processes: y


    output:
      output_data_base_name: ./Results/Grains
      output_frequency: 100
      output_variables:
        - grain_ID

time_integrators:
  - standard_time_integrator:
      name: ti_1
      termination_time: 17.81
      termination_step_count: 10000
      time_step: 0.001

      realms:
        - realm1
