#ifndef CELLULAR_AUTOMATA_MANAGER_2D_H
#define CELLULAR_AUTOMATA_MANAGER_2D_H

#include <mpi.h>
// STL include
#include <vector>
#include <list>

// Forward declarations
class Octahedron;
class Grain;
class Orientation;
class ParallelCommManager;
class ParallelCommManager_2D;
class ProblemPhysics;

// Base class
#include "CellularAutomataManager.h"

class CellularAutomataManager_2D: public CellularAutomataManager
{
public:

  // Constructor and destructor
  CellularAutomataManager_2D();
  ~CellularAutomataManager_2D();
  // load
  void load(const YAML::Node& node);
  // setup data member pointers associated with output and parallel manager.
  void setup_several_data_member_pointer();
  // initialize
  void initialize();
  // Set up problem
  virtual void setup_problem(int nDim, double* geometry, double h, int* meshInfo);
  virtual void setup_cells();

  // Create population of nucleation sites by rejection random sampling
  void create_nucleation_sites_rejection_random_sampling(double numberDensity, 
                                                               double dTcMean, 
                                                              double dTcSigma,
                                               bool surfacePopulation = false);
  // calculate time step
  double calculate_time_step();

  // Take a single timestep and update time variable
  virtual void single_integration_step();

  friend class ParallelCommManager;
  friend class ParallelCommManager_2D;
  ParallelCommManager_2D * commManager_;
  
protected:

  // Set initial values at boundary cells
  virtual void initialize_boundaries();

  // Compute number of cells, and offset, in given dimension
  virtual void get_local_size(MPI_Comm cartComm, int dim, int n, int & nLocal, int & offset);

  // Compute random angle between 0 and pi/4
  virtual void get_random_orientation_angle(double* angle);

  // Create new grain with a given cell (local id) and orientation
  // (note that this automatically pushes the grain onto grainVector_
  virtual Grain * create_new_grain(int cellID,
    Orientation * orientation,
    double * cellCenter, double length = 0);

  // Check whether growing grain captures the cell (i.e., cell center
  // in inside grain square)
  virtual bool cell_is_captured(const Grain * grain, int iCell);

  // Capture cell by creating a new grain with the appropriate orientation
  virtual void capture_cell(double xcgrain,
    double ycgrain,
    double length,
    int orientationID,
    int iCell);

  // Get (x,y) coords at the center of a cell with a given ID
  virtual void get_cell_center(int cellLocGhostID, double* cellCenter);

  // Get vector of cell neighbors based on locGhost ID
  virtual void get_cell_neighbors(int cellLocGhostID, std::vector<int> & neighborVec);

  // Get vector of bulk or surface cell IDs (locGhost ID)
  virtual void get_bulk_cells(std::vector<int> & cellIDs);
  virtual void get_surface_cells(std::vector<int> & cellIDs);

  // Test whether a given cell (locGhost ID) is a ghost cell
  virtual bool cell_is_ghost(int locGhostID);
  virtual bool cell_is_global_ghost(int locGhostID);

  // Utilities to map between various numbering and indexing
  virtual void convert_id_local_to_locghost(int localID, int & locGhostID);
  virtual void convert_id_locghost_to_local(int locGhostID, int & localID);
  virtual void convert_id_locghost_to_global(int locGhostID, int & globalID);
  virtual void convert_id_global_to_locghost(int globalID, int & locGhost);
  bool convert_id_global_to_locghost_and_check(int globalID, int & locGhost);

  // Communicate new orientations with other procs and update global orientation list
  virtual void update_global_orientation_list(const std::vector<Orientation *> newOrientations,
    const std::vector<int> updateCells);


  // Output cell data to files
  virtual void outputToFile(const std::string & fileName, int * u);
  virtual void writeToFile(const std::string & fileName, const int * u, const int * meltTag);
};

#endif