#ifndef CELLULAR_AUTOMATA_MANAGER_3D_H
#define CELLULAR_AUTOMATA_MANAGER_3D_H

#include <mpi.h>

// STL include
#include <vector>
#include <list>

// Forward declarations
class Octahedron;
class Grain;
class Orientation;
class ParallelCommManager;
class ParallelCommManager_3D;
class ProblemPhysics;
class VtrManager;

// Base class
#include "CellularAutomataManager.h"

class CellularAutomataManager_3D : public CellularAutomataManager
{
public:

  // Constructor and destructor
  CellularAutomataManager_3D( );
  virtual ~CellularAutomataManager_3D( );

  // Load
  void load(const YAML::Node& node);
  // setup data member pointers associated with output and parallel manager.
  void setup_several_data_member_pointer();
  // Initialization
  void initialize();

  // Set up problem
  void setup_problem(int nDim, double* geometry, double h, int* meshInfo);
  void setup_cells();

  // Create population of nucleation sites
  void create_nucleation_sites_rejection_random_sampling
        ( double numberDensity, double dTcMean, double dTcSigma,
          bool surfacePopulation = false);
  // calculate time step
  double calculate_time_step();

  // Take a single timestep and update time variable
   void single_integration_step();

   friend class ParallelCommManager;
   friend class ParallelCommManager_2D;
   friend class ParallelCommManager_3D;
   friend class VtrManager;
   ParallelCommManager_3D * commManager_;
   VtrManager * vtrManager_;

public:
 
  // Compute number of cells, and offset, in given dimension
  void get_local_size(MPI_Comm cartComm, int dim, int n, int & nLocal, int & offset);

  // Set initial values at boundary cells
  void initialize_boundaries();

  // Compute random angle between 0 and pi/4
  void get_random_orientation_angle(double* angle);

  // Create new grain with a given cell (local id) and orientation
  // (note that this automatically pushes the grain onto grainVector_
  virtual Grain * create_new_grain(int cellID,
                                   Orientation * orientation,
                                   double * cellCenter, 
                                   double length);

  // Capture cell by creating a new grain with the appropriate orientation
  virtual void capture_cell(double * cellCenter,
                            double length,
                            int orientationID,
                            int iCell, 
                            bool resetCoordMatrix);

  // Get (x,y) coords at the center of a cell with a given ID
  void get_cell_center(int cellLocGhostID, double* cellCenter);

  // Get vector of cell neighbors based on locGhost ID
  void get_cell_neighbors(int cellLocGhostID, std::vector<int> & neighborVec);

  // Get vector of bulk or surface cell IDs (locGhost ID)
  void get_bulk_cells(std::vector<int> & cellIDs);
  void get_surface_cells(std::vector<int> & cellIDs);

  // Test whether a given cell (locGhost ID) is a ghost cell
  bool cell_is_ghost(int locGhostID);
  bool cell_is_global_ghost(int locGhostID);

  // Utilities to map between various numbering and indexing
  void convert_id_local_to_locghost(int localID, int & locGhostID);
  void convert_id_locghost_to_local(int locGhostID, int & localID);
  void convert_id_locghost_to_global(int locGhostID, int & globalID);
  void convert_id_global_to_locghost(int globalID, int & locGhost);
  bool convert_id_global_to_locghost_and_check(int globalID, int & locGhost);

  // Communicate new orientations with other procs and update global orientation list
  void update_global_orientation_list(const std::vector<Orientation *> newOrientations,
                                      const std::vector<int> updateCells);

  // Output cell data to files
  virtual void outputToFile(const std::string & fileName, int * u);
  virtual void writeToFile(const std::string & fileName, const int * u, const int * meltTag);

  /*----------------------------------------------
     output the final microstructure information
    ----------------------------------------------*/
  void output_microstructure_information();
  void load_microstructure_information();
  virtual void random_microstructure_information();
  void output_last_step_result();

public:
  // for mapvoxel manager testing
  void checkup_current_temperature_to_update_phase_state();
  bool get_amount_of_liquid_cells(int &numberLiquidCell);
  void get_amount_of_final_orientations(int &numberOrientation);
};

#endif