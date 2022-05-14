#ifndef CELLULAR_AUTOMATA_MANAGER_H
#define CELLULAR_AUTOMATA_MANAGER_H

#include <mpi.h>
// STL include
#include <vector>
#include <list>
#include "CafeEnv.h"
#include "Simulation.h"
#include "ProblemPhysics.h"
// Forward declarations
class Grain;
class Orientation;
class Envelope;
class ParallelCommManager;
class ProblemPhysics;
class Timer;
class MapVoxelManager;

class CellularAutomataManager
{
public:
  // Constructor and destructor
  CellularAutomataManager();
  virtual ~CellularAutomataManager();
  // load
  virtual void load(const YAML::Node& node)=0;
  // setup data member pointers associated with output and parallel manager.
  virtual void setup_several_data_member_pointer()=0;
  // initialize
  virtual void initialize()=0;

  // Set up problem
  virtual void setup_problem(int nDim, double* geometry, double h, int* meshInfo) = 0;
  virtual void setup_cells()=0;


  // Setup integration parameters
  virtual void setup_integration(double maxDt);

  // Setup time step for CAFE
  virtual void setup_time_step_for_coupling(double dT);

  // Integrate up to a given time
  virtual void time_integration(double finalTime, int maxSteps);

  // Calculate time step
  virtual double calculate_time_step()=0;
  
  // Take a single timestep and update time variable
  virtual void single_integration_step()=0;
  // Set problem physics
  void set_problem_physics(ProblemPhysics * p) { problemPhysics_ = p; };
  virtual double get_lower_temperature(){return problemPhysics_->get_solidus_temperature();};
  // Create population of nucleation sites
  virtual void create_nucleation_sites(double numberDensity,
                                       double dTcMean,
                                       double dTcSigma,
                                       bool surfacePopulation = false);


  virtual void global_mapped_grain_ID(double numberDensitySurface,
                                      double numberDensityBulk);

  // Nucleate grains at nucleation sites
  virtual void nucleate_at_nucleation_sites(double time);

  // Output cell data to files
  virtual void output_orientations(const std::string & fileName);

  // Domain geometry data
  int nDim_;          // dimension size of the domain
  double h_;         // size of each cubical cell
  int nx_, ny_, nz_, nxy_; // number of cells in each direction
  double x0_, y0_, z0_;    // coordinates of corner of cell

  int xStart_, yStart_, zStart_;     // starting grid indices
  int nxLocal_, nyLocal_, nzLocal_, nxyLocal_;   // number of cells in local grid in each dimension
  int nxLocGhost_, nyLocGhost_, nzLocGhost_, nxyLocGhost_; // number of cells in grid, including ghosts, in each dimension
  int numLocalCells_;       // number of cells on this processor
  int numLocGhostCells_;    // number of cells on proc, including ghosts

  // Integration data
  double maxTimestep_;
  double time_;
  int iStep_;
  double dT_;
  double timeStepFactor_;
  double timeStep_;

  // time cost record
  Timer * CAtimer_;
  void set_timer_pointer(Timer *timer);
  
  // Cell data which should be written to a cell class in the future
  int numSitesSurface_;       // number of the nucleation sites of the overall surface
  int numSitesBulk_;          // number of the nucleation sites of the overall bulk
  int * mappedGrainID;        // Random map from orientation ID to "shuffled" IDs
  bool shafferVtrID_;            // output original orientation ID, or "shuffled" IDs
  bool shafferVtkID_;            // output original orientation ID, or "shuffled" IDs
  int * cellOrientationID_;   // GLOBAL orientation ID for cell
  int * cellGrainID_;         // LOCAL grain ID for cell
  double * cellNucleationDT_; // critical undercooling for nucleation at each cell (0 if not a nucleation site)
  double *cellTemperature_;   // cell temperature
  bool *cooling_;
  std::vector<bool> hasBeenMelted_;     // locally defined list of cells to show whether this cell has been melted
  std::list<int> nucleationCellIDList_; // List of cell IDs that are nucleation sites
  
  std::vector<bool> cellLiquidNeigh; // a cell has liquid Neigh or not
  // grain growth if a liquid cell nearing it
  virtual void checkup_liquid_cell_neigh_to_growth(double time){};

  // Grain and Orientation vectors
  Envelope * envelope_;
  std::vector<Grain *> grainVector_;             // *locally* defined list of grains
  std::list<Grain *> activeGrainList_;           // *locally* defined list of *active* grains
  std::vector<Orientation *> orientationVector_; // *globally* consistent list of all orientations
  double cellGrainLegnthCutoff_;                 //  length limitation of grain associated with a cell

  // Parallel info
  int numProcs_;
  int procID_;
  MPI_Comm cartComm_;
  int cartRank_;
  int procXP_, procXM_, procYP_, procYM_, procZP_, procZM_; // neighboring procs for communication in 4 directions

  // interpolate temperature from FE to CA
  MapVoxelManager* mapVoxelManager_;
  virtual void setup_map_voxel_manager_and_accessory(MapVoxelManager* p);
  virtual void checkup_current_temperature_to_update_phase_state() {};

  // some control variables
  bool setInitialStatedAsMelted_;
  bool setNucleationSite_;    // for remelting problem
  bool activateExistingGrain_; // for remelting problem
  bool capture_mushy_only_;   // only mushy zone cell where T_s < T < T_l can be caputare

  // model type
  double  couplingTimeStepFactor_;

protected:

  // Problem info
  ProblemPhysics * problemPhysics_;

    // Compute number of cells, and offset, in given dimension
  virtual void get_local_size(MPI_Comm cartComm, int dim, int n, int & nLocal, int & offset)=0;

  // Allocate memory for cell data
  virtual void allocate_arrays();
  
  // Set initial cell values
  virtual void initialize_cells();

  // Set initial values at boundary cells
  virtual void initialize_boundaries()=0;

  // Compute random angle between 0 and pi/4
  virtual void get_random_orientation_angle(double* angle)=0;

  // Create new grain with a given cell (local id) and orientation
  // (note that this automatically pushes the grain onto grainVector_
  virtual Grain * create_new_grain(int cellID,
                           Orientation * orientation,
                         double * cellCenter, double length=0);

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
  virtual void get_cell_center(int cellLocGhostID, double* cellCenter)=0;

  // Get vector of cell neighbors based on locGhost ID
  virtual void get_cell_neighbors(int cellLocGhostID, std::vector<int> & neighborVec)=0;

  // Get vector of bulk or surface cell IDs (locGhost ID)
  virtual void get_bulk_cells(std::vector<int> & cellIDs)=0;
  virtual void get_surface_cells(std::vector<int> & cellIDs)=0;
  
  // Test whether a given cell (locGhost ID) is a ghost cell
  virtual bool cell_is_ghost(int locGhostID)=0;
  virtual bool cell_is_global_ghost(int locGhostID)=0;

  // Utilities to map between various numbering and indexing
  virtual void convert_id_local_to_locghost(int localID, int & locGhostID)=0;
  virtual void convert_id_locghost_to_local(int locGhostID, int & localID)=0;
  virtual void convert_id_locghost_to_global(int locGhostID, int & globalID)=0;
  virtual void convert_id_global_to_locghost(int globalID, int & locGhost)=0;
  virtual bool convert_id_global_to_locghost_and_check(int globalID, int & locGhost)=0;

  // Communicate new orientations with other procs and update global orientation list
  virtual void update_global_orientation_list(const std::vector<Orientation *> newOrientations,
                                      const std::vector<int> updateCells)=0;

  // Output cell data to files
  virtual void outputToFile(const std::string & fileName, int * u)=0;
  virtual void writeToFile(const std::string & fileName, const int * u, const int *meltTag)=0;

  /*----------------------------------
     nucleation sites density and
     probability density functions
   -----------------------------------*/
  double ns_;           // surface nucleation site density
  double deltaTs_max_;  // mean undercooling for surface sites
  double deltaTs_sigma_;// standard deviation of undercooling for surface sites

  double nv_;           // bulk nucleation site density
  double deltaTv_max_;  // mean undercooling for bulk sites
  double deltaTv_sigma_;// standard deviation of undercooling for bulk sites
  /*----------------------------------
     solution control parameters
    ----------------------------------*/
  bool first_nearnest_neighbours;  // specify neighbours

public:
  /*----------------------------------------------
     output/load the final microstructure information
    ----------------------------------------------*/
  bool output_nucleation_seeds_;
  bool load_nucleation_seeds_;
  int mappedGrainIDSizeFromFile_;
  std::string microstructureInformationFileName_load_;
  std::string microstructureInformationFileName_output_;
  bool outputMicrostructureInfo_;
  bool loadMicrostructureInfo_;  
  bool randomMicrostructureInfo_;   // random generate microstructure info first
  virtual void output_microstructure_information();
  virtual void load_microstructure_information();
  virtual void random_microstructure_information() {};
  virtual void output_last_step_result();

  /*----------------------------------------------------------------
      statistics methods
    ----------------------------------------------------------------*/
  bool solidificationDone_;
  virtual void get_amount_of_final_orientations(int &numberOrientation) {};
  virtual bool get_amount_of_liquid_cells(int &numberLiquidCell) { return false; };
};

#endif
