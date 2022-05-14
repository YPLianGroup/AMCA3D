#ifndef MANAGER
#define MANAGER

#include "mpi.h"


// Local includes

class Manager
{
public:

	// Constructor and destructor
	Manager();
	~Manager();

	// Set up problem
	virtual void setup_cells() =0;

	// Nucleate initial grains and set up nucleation sites
	virtual void initialize_nucleation()=0;

	// Set up integration parameters
	virtual void setup_integration(double maxDt)=0;

	// Integrate up to a given time
	virtual void time_integration(double finalTime, int maxSteps) =0;

	// Take a single timestep and update time variable
	virtual void single_integration_step(double & time) =0;

	// Create population of nucleation sites
	virtual void create_nucleation_sites(double numberDensity,
		double dTcMean,
		double dTcSigma,
		bool surfacePopulation = false)=0;

	// Nucleate grains at nucleation sites
	virtual void nucleate_at_nucleation_sites(double time)=0;

	// Output cell data to files
	virtual void output_orientations(const std::string & fileName)=0;

	// Domain geometry data
	double h_;             // size of each cubical cell
	int nx_, ny_, nz;      // number of cells in each direction
	double x0_, y0_, z0_;       // coordinates of corner of cell

					   //int numProcX_, numProcY_; // number of procs in each dimension
	int xStart_, yStart_;     // starting grid indices
	int nxLocal_, nyLocal_;   // number of cells in local grid in each dimension
	int nxLocGhost_, nyLocGhost_; // number of cells in grid, including ghosts, in each dimension
	int numLocalCells_;       // number of cells on this processor
	int numLocGhostCells_;    // number of cells on proc, including ghosts

							  // Integration data
	double maxTimestep_;

	// Cell data
	int * cellOrientationID_; // GLOBAL orientation ID for cell
	int * cellGrainID_;       // LOCAL grain ID for cell
	double * cellNucleationDT_; // critical undercooling for nucleation at each cell (0 if not a nucleation site)

								// PETSc objects
								//DM distributedArray_; // Structured grid array info
								//Vec cellDataVec_;     // PETSc vector to hold cell data

								// Grain and Orientation vectors
	std::vector<Grain *> grainVector_;             // *locally* defined list of grains
	std::list<Grain *> activeGrainList_;           // *locally* defined list of *active* grains
	std::vector<Orientation *> orientationVector_; // *globally* consistent list of all orientations

												   // List of cell IDs that are nucleation sites
	std::list<int> nucleationCellIDList_;

	// Parallel info
	int numProcs_;
	int procID_;
	MPI_Comm cartComm_;
	int cartRank_;
	int procXP_, procXM_, procYP_, procYM_; // neighboring procs for communication in 4 directions
	friend class ParallelCommManager;
	ParallelCommManager * commManager_;


protected:

	// Problem info
	ProblemPhysics * problemPhysics_;

	// Compute number of cells, and offset, in given dimension
	void get_local_size(MPI_Comm cartComm, int dim, int n, int & nLocal, int & offset);

	// Allocate memory for cell data
	void allocate_arrays();

	// Set initial cell values
	void initialize_cells();

	// Set initial values at boundary cells
	void initialize_boundaries();

	// Nucleate grains that are present from the initial time
	void nucleate_initial_grains();

	// Compute random angle between 0 and pi/4
	void get_random_orientation_angle(double & angle);

	// Create new grain with a given cell (local id) and orientation
	// (note that this automatically pushes the grain onto grainVector_
	Grain * create_new_grain(int cellID,
		Orientation * orientation,
		double xc, double yc, double length = 0);

	// Check whether growing grain captures the cell (i.e., cell center
	// in inside grain square)
	bool cell_is_captured(const Grain * grain, int iCell);

	// Capture cell by creating a new grain with the appropriate orientation
	void capture_cell(double xcgrain,
		double ycgrain,
		double length,
		int orientationID,
		int iCell);

	// Get (x,y) coords at the center of a cell with a given ID
	void get_cell_center(int cellLocGhostID, double & xc, double & yc);

	// Get vector of cell neighbors based on locGhost ID
	void get_cell_neighbors(int cellLocGhostID, std::vector<int> & neighborVec);

	// Get vector of bulk or surface cell IDs (locGhost ID)
	void get_bulk_cells(std::vector<int> & cellIDs);
	void get_surface_cells(std::vector<int> & cellIDs);

	// Test whether a given cell (locGhost ID) is a ghost cell
	bool cell_is_ghost(int locGhostID);

	// Utilities to map between various numbering and indexing
	void convert_id_local_to_locghost(int localID, int & locGhostID);
	void convert_id_locghost_to_local(int locGhostID, int & localID);
	void convert_id_locghost_to_global(int locGhostID, int & globalID);
	void convert_id_global_to_locghost(int globalID, int & locGhost);

	// Communicate new orientations with other procs and update global orientation list
	void update_global_orientation_list(const std::vector<Orientation *> newOrientations,
		const std::vector<int> updateCells);


	// Output cell data to files
	void outputToFile(const std::string & fileName, int * u);
	void writeToFile(const std::string & fileName, const int * u, int nx, int ny);

};

#endif

