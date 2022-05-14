#ifndef Cellular_Automata_Manager_3DRemelting_H
#define Cellular_Automata_Manager_3DRemelting_H

#include "CellularAutomataManager_3D.h"
// forward class
class Grain;
class Orientation;

class CellularAutomataManager_3DRemelting: public CellularAutomataManager_3D
{
public:
  // Constructor and destructor
  CellularAutomataManager_3DRemelting();
  virtual ~CellularAutomataManager_3DRemelting();
    
  double nullMatTemp_;  // temperature value stands for null material
  // Load
  void load(const YAML::Node& node);

  // setup data member pointers associated with output and parallel manager.
  void setup_several_data_member_pointer();

  // Initialization
  void initialize();

  // Load microstructure information from file
  void load_microstructure_information(int &numberGlobalExpectedGrains);

  // Retrieve the grain information for cells without out material in previous simulation
  void retrieve_grain_information_before_simulation();

  // Integrate up to a given time
  void time_integration(double finalTime, int maxSteps);

  // update cell grain information
  void checkup_current_temperature_to_update_phase_state();
  void checkup_current_temperature();
  // grain growth if a liquid cell nearing it
  void checkup_cell_neigh_to_growth(double time);
  void checkup_liquid_cell_neigh_to_growth(double time);

  // Activate the solid-liquid interface cells
  void maintain_active_grains_on_interface_cells();
  // Calculate time step
  double calculate_time_step();

  // Take a single timestep and update time variable
  void single_integration_step();

  // capture cell
  void capture_cell(double* grainCenter, double length,
                    int orientationID, int iCell, 
                    bool resetCoordMatrix);

  // create new grain
  Grain * create_new_grain(int cellID, Orientation * orientation,
                           double * cellCenter, double length);

  void global_mapped_grain_ID_plus(double numberDensitySurface,
                                   double numberDensityBulk,
                                   int numberExistingGrains);

  // Nucleate grains at nucleation sites
  void nucleate_at_nucleation_sites(double time);
};

#endif
