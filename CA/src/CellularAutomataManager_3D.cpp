// General includes
#include <iostream>
#include <fstream>
#include <random>
#include <mpi.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <Timer.h>
#include <stdlib.h>
#include <algorithm>

// Local includes
#include "CellularAutomataManager_3D.h"
#include "Octahedron.h"
#include "Grain.h"
#include "Orientation.h"
#include "ParallelCommManager_3D.h"
#include "ProblemPhysics.h"
#include "VtrManager.h"
#include "CafeParsing.h"
#include "MapVoxelManager.h"

#define M_PI 3.14159265358979323846 
/*=======================================================================
    Class Definition
      CellularAutomataManager_3D - manage 3D CA
  =======================================================================*/
/*-------------------
    Constructor
  -------------------*/
CellularAutomataManager_3D::CellularAutomataManager_3D()
  : CellularAutomataManager()
{
  // Parallel setup
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs_);
  MPI_Comm_rank(MPI_COMM_WORLD, &procID_);  

  // Random seed
  srand(procID_ * 10000);
  envelope_ = new Octahedron;
}

/*------------------
   Destructor
  ------------------*/ 
CellularAutomataManager_3D::~CellularAutomataManager_3D()
{
  if (commManager_)
    delete commManager_;
  if (vtrManager_)
    delete vtrManager_;
  if (envelope_)
    delete envelope_;
}

/*--------------------------------------
   setup_several_data_member_pointer
  --------------------------------------*/
void
CellularAutomataManager_3D::setup_several_data_member_pointer()
{
  commManager_ = new ParallelCommManager_3D(this);
  vtrManager_ = new VtrManager(this);
}

/*-----------------
   Load
  -----------------*/
void
CellularAutomataManager_3D::load(const YAML::Node & node)
{
  double original_point[3] = { -1 - 1 - 1 };
  double lateral_size[3] = { -1 - 1 - 1 };
  double dCell;
  int    numbCells[3] = { -1 - 1 - 1 };
  
  // get domain size
  const YAML::Node * domain_node = node.FindValue("domain");
  if (domain_node)
  {
    CafeEnv::self().caOutputP0() << "\n" << "Domain Size Review" << "\n";
    CafeEnv::self().caOutputP0() << "=============================" << "\n";
    std::string fileName;
    get_if_present_no_default(*domain_node, "load_microstructure_information", fileName);
    if(fileName=="")
    {
      loadMicrostructureInfo_ = false;
    }
    else
    {
      microstructureInformationFileName_load_ = fileName;
      loadMicrostructureInfo_ = true;
      // Prepare the output file name
      std::size_t dotPos = microstructureInformationFileName_load_.find_last_of(".");
      std::string extension = std::to_string(procID_) + ".txt";
      fileName = microstructureInformationFileName_load_.substr(0, dotPos);
      fileName += extension;
      std::fstream inputStream;
      inputStream.open(fileName.c_str(), std::ios::in);
      if (!inputStream) {
        std::cout << "fail to open the file" << fileName << std::endl;
        throw std::runtime_error("fail to open the specified file");
      }
      double X=-1, Y=-1, Z=-1, nx=-1, ny=-1, nz=-1;
      int numCells, numOrientations, numLocalGrains, numGlobalExpectedGrains;
      std::string line;
      std::getline(inputStream, line); // skip the description line
      std::getline(inputStream, line); // get the value of the header
      std::istringstream iss(line);
      iss >> numCells >> numOrientations >> numLocalGrains >> numGlobalExpectedGrains >> X >> Y >> Z >> nx >> ny >> nz >> dCell;
      if(nx<0.5)
      {
        CafeEnv::self().caOutputP0() << "parser error: realm-domain size from file: " << fileName << "\n";
        loadMicrostructureInfo_ = false;
      }
      else
      {
        CafeEnv::self().caOutputP0() << "load domain size from file: " << fileName << "\n";
        original_point[0] = X;
        original_point[1] = Y;
        original_point[2] = Z;
        lateral_size[0] = nx*dCell;
        lateral_size[1] = ny*dCell;
        lateral_size[2] = nz*dCell;
        numbCells[0] = nx;
        numbCells[1] = ny;
        numbCells[2] = nz;
      }
    }
  }
  else
    throw std::runtime_error("parser error: realm-domain");

  if (domain_node && !loadMicrostructureInfo_)
  {
    // Make some things pretty in the log file
    get_required(*domain_node, "original_point", original_point);
    //get_if_present(*domain_node, "lateral_sizes", lateral_size, lateral_size);
    get_if_present_no_default(*domain_node, "lateral_sizes", lateral_size);
  }
  else if(!loadMicrostructureInfo_)
    throw std::runtime_error("parser error: realm-domain");
  CafeEnv::self().caOutputP0() << "Original Point: " << "\n";
  CafeEnv::self().caOutputP0() << "    x0=" << original_point[0] 
                                << "  y0=" << original_point[1] 
                                << "  z0=" << original_point[2] << "\n";
  if (lateral_size[0] > 0)
  {
    CafeEnv::self().caOutputP0() << "Lateral Size: " << "\n";
    CafeEnv::self().caOutputP0() << "    Lx=" << lateral_size[0]
                                  << "    Ly=" << lateral_size[1]
                                  << "    Lz=" << lateral_size[2] << "\n";
  }
  // get discretization info
  const YAML::Node * discretization_node = node.FindValue("discretization");
  if (discretization_node && !loadMicrostructureInfo_)
  {
    // Make some things pretty in the log file
    CafeEnv::self().caOutputP0() << "\n" << "Discretization Review" << "\n";
    CafeEnv::self().caOutputP0() << "=============================" << "\n";

    get_required(*discretization_node, "cell_size", dCell);
    if (lateral_size[0] < 0)
       get_required(*discretization_node, "number_cells", numbCells);
    else
    {
      for (int i = 0; i < 3; i++)
      {
        numbCells[i] = lateral_size[i] / dCell;
        lateral_size[i] = numbCells[i] * dCell;
      }
    }
  }
  else if(!loadMicrostructureInfo_)
    throw std::runtime_error("parser error: realm-discretization");
  CafeEnv::self().caOutputP0() << " Cell Size: " << dCell << "\n";
  CafeEnv::self().caOutputP0() << " Number Cells: x=" << numbCells[0]
                                << " y=" << numbCells[1]
                                << " z=" << numbCells[2] << "\n";
  // original point coordinates
  x0_ = original_point[0];
  y0_ = original_point[1];
  z0_ = original_point[2];

  // mesh information
  h_ = dCell;
  nx_ = numbCells[0];
  ny_ = numbCells[1];
  nz_ = numbCells[2];

  // nucleation sites density and pdf
  const YAML::Node * nucleation = node.FindValue("nucleation_rules");
  if (nucleation)
  {
    for (size_t iterator = 0; iterator < nucleation->size(); iterator++)
    {
      const YAML::Node & nucleation_node = (*nucleation)[iterator];
      const YAML::Node * surface_node = nucleation_node.FindValue("surface");
      if (surface_node)
      {
        // Make some things pretty in the log file
        CafeEnv::self().caOutputP0() << "\n" << "Surface Nucleation Rule Review" << "\n";
        CafeEnv::self().caOutputP0() << "=============================" << "\n";

        std::string type = "Gaussian";
        get_required(*surface_node, "type", type);
        get_required(*surface_node, "site_density", ns_);
        get_required(*surface_node, "mean", deltaTs_max_);
        get_required(*surface_node, "standard_deviation", deltaTs_sigma_);

        CafeEnv::self().caOutputP0() << "PDF: " << type << "\n";
        CafeEnv::self().caOutputP0() << "Site Density: " << ns_ << "\n";
        CafeEnv::self().caOutputP0() << "Mean: " << deltaTs_max_ << "\n";
        CafeEnv::self().caOutputP0() << "Standard Deviation: " << deltaTs_sigma_ << "\n";
      }
      const YAML::Node * bulk_node = nucleation_node.FindValue("bulk");
      if (bulk_node)
      {
        // Make some things pretty in the log file
        CafeEnv::self().caOutputP0() << "\n" << "Bulk Nucleation Rule Review" << "\n";
        CafeEnv::self().caOutputP0() << "=============================" << "\n";

        std::string type = "Gaussian";
        get_if_present(*bulk_node, "type", type, type);
        get_required(*bulk_node, "site_density", nv_);
        get_required(*bulk_node, "mean", deltaTv_max_);
        get_required(*bulk_node, "standard_deviation", deltaTv_sigma_);

        CafeEnv::self().caOutputP0() << "PDF: " << type << "\n";
        CafeEnv::self().caOutputP0() << "Site Density: " << nv_ << "\n";
        CafeEnv::self().caOutputP0() << "Mean: " << deltaTv_max_ << "\n";
        CafeEnv::self().caOutputP0() << "Standard Deviation: " << deltaTv_sigma_ << "\n";
      }
    }
  }
  else
  {
      CafeEnv::self().caOutputP0() << "\n" << "Nucleation Information Review" << "\n";
      CafeEnv::self().caOutputP0() << "=============================" << "\n";
      CafeEnv::self().caOutputP0() << "\n" << "Nucleation as default values" << "\n";
      CafeEnv::self().caOutputP0() << "Surface Nucleation Site Density: " << ns_ << "\n";
      CafeEnv::self().caOutputP0() << "Bulk Nucleation Site Density: " << nv_ << "\n";
  }
  //throw std::runtime_error("parser error: realm-nucleation_rules");

  // solution options
  const YAML::Node * solutionOptions = node.FindValue("solution_options");
  if (solutionOptions)
  {
    // Make some things pretty in the log file
    CafeEnv::self().caOutputP0() << "\n" << "Solution Options Review" << "\n";
    CafeEnv::self().caOutputP0() << "=============================" << std::endl;
    std::string name_;
    get_required(*solutionOptions, "name", name_);

    const YAML::Node * options_node = expect_sequence(*solutionOptions, "options", true);
    if (options_node)
    {
      for (size_t ioption = 0; ioption < options_node->size(); ioption++)
      {
        const YAML::Node & option_node = (*options_node)[ioption];
        if (expect_map(option_node, "output_microstructure_information", true))
        {
          const YAML::Node& fileName = *option_node.FindValue("output_microstructure_information");
          get_required(fileName, "file_name", microstructureInformationFileName_output_);
          outputMicrostructureInfo_ = true;
          get_if_present(fileName, "output_seeds", output_nucleation_seeds_, output_nucleation_seeds_);
          CafeEnv::self().caOutputP0() << "Output final microstructure information to file: "
                                       << microstructureInformationFileName_output_ << std::endl
                                       << "Output nucleations seeds: "<< output_nucleation_seeds_<< std::endl;
        }
        else if (expect_map(option_node, "random_type", true))
        {
          const YAML::Node& randomType = *option_node.FindValue("random_type");
          get_if_present(randomType, "shuffer_vtr_grain_ID", shafferVtrID_,shafferVtrID_);
          get_if_present(randomType, "shuffer_vtk_grain_ID", shafferVtkID_,shafferVtrID_);
          get_if_present(randomType, "random_microstructure_information", randomMicrostructureInfo_,randomMicrostructureInfo_);
          CafeEnv::self().caOutputP0() << "shuffer vtr grain_ID: "<< shafferVtrID_ << std::endl;
          CafeEnv::self().caOutputP0() << "shuffer vtk grain_ID: "<< shafferVtkID_ << std::endl;
          CafeEnv::self().caOutputP0() << "random microstructure info: " << randomMicrostructureInfo_ << std::endl;
        }
        else if (expect_map(option_node, "phase_state", true))
        {
          const YAML::Node& phaseState = *option_node.FindValue("phase_state");
          get_if_present(phaseState, "has_been_melted", setInitialStatedAsMelted_, setInitialStatedAsMelted_);
          CafeEnv::self().caOutputP0() << "phase state: has been melted = " << setInitialStatedAsMelted_ << std::endl;
        }
        else if (expect_map(option_node, "cell_grain_length_cutoff_factor", true))
        {
          const YAML::Node& cutOff = *option_node.FindValue("cell_grain_length_cutoff_factor");
          get_if_present(cutOff, "factor_value", cellGrainLegnthCutoff_, cellGrainLegnthCutoff_);          
        }
        else if (expect_map(option_node, "parallel_specification", true)) {
          const YAML::Node& parallel = *option_node.FindValue("parallel_specification");
          std::string axis = "o";
          get_if_present(parallel, "direction_with_all_processes", axis, axis);
          if (axis != "o") {
            commManager_->controledByInputFile_ = true;            
            int specifiedDirection = 0;
            if (axis == "y") {
              specifiedDirection = 1;
            }
            else if (axis == "z") {
              specifiedDirection = 2;
            }
            commManager_->specifiedDirection_ = specifiedDirection;
          }
        }
        else if (expect_map(option_node, "envelope", true)) {
          const YAML::Node& envelopeInfor = *option_node.FindValue("envelope");
          double diagonalLengthFactor = -1;
          get_if_present(envelopeInfor, "max_half_diagonal_length", diagonalLengthFactor, diagonalLengthFactor);
          if (diagonalLengthFactor > 0) {
            envelope_->set_truncation_coefficient(diagonalLengthFactor);
          }
        }
        else if (expect_map(option_node, "load_microstructure_information", true)) {
            if(randomMicrostructureInfo_)
                continue;
            else
            {
                const YAML::Node& fileName = *option_node.FindValue("load_microstructure_information");
                get_required(fileName, "file_name", microstructureInformationFileName_load_);
                if (microstructureInformationFileName_load_ == "random")
                {
                    loadMicrostructureInfo_ = false;
                    randomMicrostructureInfo_ = true;
                }
                else {
                    loadMicrostructureInfo_ = true;
                    get_if_present(fileName, "load_seeds", load_nucleation_seeds_, load_nucleation_seeds_);
                    if (load_nucleation_seeds_) {
                        setNucleationSite_ = false;
                    }
                }
            }
        }
      }
    }
  }

  if(randomMicrostructureInfo_^loadMicrostructureInfo_) {
      if(randomMicrostructureInfo_)
          CafeEnv::self().caOutputP0() << "random microstructure info: " << randomMicrostructureInfo_ << std::endl;
      else
          CafeEnv::self().caOutputP0() << "Load final microstructure information to file: "
            << microstructureInformationFileName_load_
            << "Load nucleation seeds: " << load_nucleation_seeds_ << std::endl;
  }
  else{
      if(randomMicrostructureInfo_)
          loadMicrostructureInfo_ = false;
      else{
          randomMicrostructureInfo_ = true;
          CafeEnv::self().caOutputP0() << "NO microstructure info, so random generate are used to initialize microstructure!!! " << std::endl;
      }
  }
  // get output
  vtrManager_->load(node);
}

/*----------------------
   initialize
  ----------------------*/
void CellularAutomataManager_3D::initialize()
{
  nxy_ = nx_ * ny_;
  setup_cells();
  // Create nucleation sites on surface
  create_nucleation_sites(ns_, deltaTs_max_, deltaTs_sigma_, true);
  // Create nucleation sites in volume
  create_nucleation_sites(nv_, deltaTv_max_, deltaTv_sigma_, false);
  // generate global mapped grain ID
  global_mapped_grain_ID(ns_, nv_);
  // set up the length cutoff
  cellGrainLegnthCutoff_ *= h_;
  vtrManager_->initialize();
}

/*-----------------
   Set up problem
  -----------------*/
void CellularAutomataManager_3D::setup_problem(int nDim, double* geometry, double h, int* meshInfo)
{
  // geometry info
  nDim_ = nDim;
  x0_ = geometry[0];
  y0_ = geometry[1];
  z0_ = geometry[2];
  
  // mesh information
  h_ = h;
  nx_ = meshInfo[0];
  ny_ = meshInfo[1];
  nz_ = meshInfo[2];
  
  nxy_ = nx_ * ny_;
}

/*---------------
   setup_cells
  ---------------*/ 
void
CellularAutomataManager_3D::setup_cells()
{
  // Set up Cartesian topology
  int dims[3];
  commManager_->setup_cartesian_comm();
  cartComm_ = commManager_->get_cart_comm();

  // Compute size of local grid
  get_local_size(cartComm_, 0, nx_, nxLocal_, xStart_);
  get_local_size(cartComm_, 1, ny_, nyLocal_, yStart_);
  get_local_size(cartComm_, 2, nz_, nzLocal_, zStart_);

  nxyLocal_ = nxLocal_ * nyLocal_;
  numLocalCells_ = nxyLocal_ * nzLocal_;
  
  nxLocGhost_ = nxLocal_ + 2;
  nyLocGhost_ = nyLocal_ + 2;
  nzLocGhost_ = nzLocal_ + 2;
  nxyLocGhost_ = nxLocGhost_ * nyLocGhost_;
  numLocGhostCells_ = nxyLocGhost_ * nzLocGhost_;

  commManager_->setup_indexed_array();

  // Allocate data
  allocate_arrays();

  // Initialize cells
  initialize_cells();
  initialize_boundaries();

}

/*-------------------------
   calculate_time_step
  -------------------------*/
double
CellularAutomataManager_3D::calculate_time_step()
{
  CAtimer_->start_timing_nucleation();
  // Nucleate new grains based on current undercooling
  nucleate_at_nucleation_sites(time_);
  CAtimer_->stop_timing_nucleation();

  // Loop over grains to compute timestep
  double maxVelocity = 1e-10;
  for (std::list<Grain *>::iterator iter = activeGrainList_.begin();
    iter != activeGrainList_.end();
    ++iter)
  {
    Grain * grain = *iter;
    double vel = 0.0;

    // for CAFE solver
    if (mapVoxelManager_)
    {
      double T;
      assert(mapVoxelManager_->evaluate_temperature(grain->cell_, T));

      // check phase state
      if (T >= problemPhysics_->get_melting_temperature())
      {
        hasBeenMelted_[grain->cell_] = true;
        vel = 0.0;
      }
      else if (hasBeenMelted_[grain->cell_])
      {
        /*-----------------------------------------
          In the application to AM, there is a case that a cell with grain on is 
          surrouned by cells that has never been melted. Consequencely, the cell
          will never been set as inactive as it alwasys has neighbour cells that 
          are not captured by it. Moreworse, if the temperature of the cell keeps 
          decreasing, the growth rate of grain on this cell will increase all the way.
          On the other hand, you can not set the grain as inactive based on the unmelted
          neighbour cells because those cells also get the chance to be melted at some
          time point.

          In order to avoid such a numerical-causing fake phenomenon, we set the grain
          length cutoff. 
        */
        if (grain->length_ < cellGrainLegnthCutoff_)
        {
          vel = problemPhysics_->compute_growth_velocity_with_temperature_input(T);
        }
        else
        {
          // Do not consider the growth rate of the grain on this cell!
          grain->length_ = cellGrainLegnthCutoff_; 
          vel = 0.0;
        }
      }
    } // end for CAFE solver
    // for CA solver only
    else
    {
      vel = problemPhysics_->compute_growth_velocity(grain, time_);
    }

    if (vel > maxVelocity)
    {
      maxVelocity = vel;
    }
  }
  double dtLoc = std::min(timeStepFactor_*h_ / maxVelocity, maxTimestep_);
  MPI_Allreduce(&dtLoc, &dT_, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  
  return dT_;
}
/*-------------------------
   single_integration_step
  -------------------------*/
void
CellularAutomataManager_3D::single_integration_step()
{

  CAtimer_->start_timing_capture();
  // Grow grains
  for(std::list<Grain *>::iterator iter = activeGrainList_.begin();
      iter != activeGrainList_.end();
      ++iter)
  {
    Grain * grain = *iter;
    double vel;
    if (mapVoxelManager_)
    {
      double T;
      assert( mapVoxelManager_->evaluate_temperature(grain->cell_, T));
      vel = problemPhysics_->compute_growth_velocity_with_temperature_input(T);
    }
    else
    {
      vel = problemPhysics_->compute_growth_velocity(grain, time_);
    }
    grain->grow(vel, dT_);
  }

  // Check for capture of neighbors
  std::list<Grain *>::iterator iter = activeGrainList_.begin();

  double meltingTemperature = problemPhysics_->get_melting_temperature();
  while(iter != activeGrainList_.end())
  {
    Grain * grain = *iter;
    int iCell = grain->cell_;
    std::vector<int> neighborVec;
    get_cell_neighbors(iCell, neighborVec);
    bool hasNonCapturedNeighbors = false;
    for (int i = 0; i < neighborVec.size(); ++i)
    {
      int iNeigh = neighborVec[i];
      bool nonVoid = true;
      bool nonGrain = (cellGrainID_[iNeigh] == -1); // it is now not a grain

      // check up temperature first
      bool mushyTemperature = false;
      double temperature;

      // for CAFE solver
      if (mapVoxelManager_)
      {
        // interpolate temperature from FEM result
        if (mapVoxelManager_->evaluate_temperature(iNeigh, temperature))
        {
          if (temperature >= meltingTemperature)
          {
            hasBeenMelted_[iNeigh] = true;
            mushyTemperature = false;        // confirm this
            hasNonCapturedNeighbors = true;
          }
          else if (hasBeenMelted_[iNeigh])
          {
            mushyTemperature = true;
          }
          else
          {
            hasNonCapturedNeighbors = true;
          }
        }
        else
        {
          // void cell: without element tag and temperature
          nonVoid = false;
          hasNonCapturedNeighbors = true;
        }
      } // end for CAFE solver
      // for CA solver
      else
      {
        mushyTemperature = true;
      }

      if (nonGrain && nonVoid && mushyTemperature)
      {
        double cellCenter[3];
        get_cell_center(iNeigh, cellCenter);
        bool isCaptured = envelope_->cell_is_captured(grain, cellCenter);

        // This is where we'll need something to handle parallelism, eventually
        if (isCaptured)
        {
          bool cellIsGhost = cell_is_ghost(iNeigh);
          double grainCenter[3];
          grainCenter[0] = grain->xc_;
          grainCenter[1] = grain->yc_;
          grainCenter[2] = grain->zc_;
          if (cellIsGhost)
          {
            commManager_->pack_for_comm(grainCenter,
                                        grain->length_,
                                        grain->orientation_->id_,
                                        iNeigh);
          }
          else
          {
            capture_cell(grainCenter, 
                         grain->length_, 
                         grain->orientation_->id_, 
                         iNeigh, false);
          }
        }
        else
        {
          hasNonCapturedNeighbors = true;
        }
        
      } // end if (nonGrain nonVoid mushyTemperature)
    } // end neighbors loop

    // Deactivate grain if it has no more liquid neighbors; otherwise
    // update iterator
    if (!hasNonCapturedNeighbors)
    {
      grain->inactivate();
      iter = activeGrainList_.erase(iter);
    }
    else
    {
      ++iter;
    }
    
  }
 
  CAtimer_->stop_timing_capture();
  CAtimer_->start_timing_crossProcCapture();

  // Parallel communication, and handle captures across proc boundaries
  commManager_->cross_proc_captures();
  
  CAtimer_->stop_timing_crossProcCapture();
  // Update time
  time_ += dT_;
  iStep_++;
  CAtimer_->stop_timing_CA_simulation();
  // post-processing
    int numberLiquidCells;
    get_amount_of_liquid_cells(numberLiquidCells);
    // call vtrManager
    vtrManager_->execute(iStep_, time_);
    CafeEnv::self().caOutputP0() << "************************************************"
                                 << std::endl;
    CafeEnv::self().caOutputP0() << "Time Step Count: "<< iStep_
                                 <<"  Current Time: " << time_
                                 <<"  dtN:"<< dT_ << "\n"
                                 <<" Remaining Liquid Cells: "<<numberLiquidCells<< std::endl;
  
}

/*------------------------------------------------------------------
   Capture cell by creating a new grain with the appropriate 
     orientation
   capture_cell
  -----------------------------------------------------------------*/
void 
CellularAutomataManager_3D::capture_cell(double* grainCenter,
  double length,
  int orientationID,
  int iCell, bool resetCoordMatrix)
{
  Orientation * orientation = orientationVector_[orientationID];
  if (resetCoordMatrix)
    envelope_->set_coordinate_transformation_matrix(orientation->EulerAngle_);
  double cellCenter[3];   // the coordinates of cell center
  get_cell_center(iCell, cellCenter);
  
  double gcCellCenter[3]; // the coordinates of cell center relative to grain center
  for (int i = 0; i < 3; i++)
  {
    gcCellCenter[i] = cellCenter[i] - grainCenter[i];
  }
  
  double newGrainInfo[4] = { 0.0, 0.0, 0.0, 0.0 };
  envelope_->capture_cell(gcCellCenter, length, h_, newGrainInfo);

  // Create new grain with same orientation and assign it to this cell
  double newGrainCenter[3];
  for (int i=0; i<3; i++)
      newGrainCenter[i] = newGrainInfo[i] + grainCenter[i];

  double newGrainLength = newGrainInfo[3];

  Grain * newGrain = create_new_grain(iCell, orientation, newGrainCenter, newGrainLength);
  cellGrainID_[iCell] = newGrain->id_;
  cellOrientationID_[iCell] = orientationID;
}

/*------------------
   get_local_size
  ------------------*/
void
CellularAutomataManager_3D::get_local_size(MPI_Comm cartComm, int dim, int n, int & nLocal, int & offset)
{
  // Get my process coordinates
  int cartDims[3];
  int periods[3];
  int cartCoords[3];
  MPI_Cart_get(cartComm, 3, cartDims, periods, cartCoords);

  // Compute number of local interior nodes, and offset
  // global_I = local_i + offset
  nLocal = n/cartDims[dim];
  offset = cartCoords[dim] * nLocal;
  // Account for non-uniform distribution on last proc: offset + nlocal = n
  if (cartCoords[dim] == cartDims[dim] - 1)
  {
    nLocal = n - offset;
  }
}

/*------------------------
   initialize_boundaries
  ------------------------*/
void
CellularAutomataManager_3D::initialize_boundaries()
{
  // Set all ghost boundaries to have 0 grain and orientation id (so that
  // they are treated as non-liquid by their neighbors)
  
  // Get my process coordinates
  int cartDims[3];
  int periods[3];
  int cartCoords[3];
  MPI_Cart_get(cartComm_, 3, cartDims, periods, cartCoords);

  int nx = nxLocGhost_;
  int ny = nyLocGhost_;
  int nz = nzLocGhost_;
  int nxy = nxyLocGhost_;

  // Lower BC: coords of current proc is with z=0
  if (cartCoords[2] == 0)
  {
    int k = 0;
    for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++)
      {
        int ID = nxy * k + nx * j + i;
        cellOrientationID_[ID] = 0;
        cellGrainID_[ID] = 0;
      }
  }

  // Upper BC: coords of current proc is with z=Top
  if (cartCoords[2] == cartDims[2] - 1)
  {
    int k = (nz-1);
    for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++)
      {
        int ID = nxy * k + nx * j + i;
        cellOrientationID_[ID] = 0;
        cellGrainID_[ID] = 0;
      }
  }

  // Left BC: x0
  if (cartCoords[0] == 0)
  {
    for (int k = 0; k < nz; k++)
    {
      for (int j = 0; j < ny; j++)
      {
        int i = 0;
        int ID = nxy * k + nx * j + i;
        cellOrientationID_[ID] = 0;
        cellGrainID_[ID] = 0;
      }
    }
  }

  // Right BC: x=nx*h
  if (cartCoords[0] == cartDims[0] - 1)
  {
    for (int k = 0; k < nz; k++)
    {
      for (int j = 0; j < ny; j++)
      {
        int i = nx - 1; 
        int ID = nxy * k + nx * j + i;
        cellOrientationID_[ID] = 0;
        cellGrainID_[ID] = 0;
      }
    }
  }

  // Front: y0
  if (cartCoords[1] == 0)
  {
      for (int k=0; k<nz; k++)
      for (int j=0; j<1; j++)
        for (int i = 0; i < nx; i++)
        {
          int ID = nxy * k + nx * j + i;
          cellOrientationID_[ID] = 0;
          cellGrainID_[ID] = 0;
        }
  }

  // Rear: y=ny*h
  if (cartCoords[1] == cartDims[1] -1)
  {
    for (int k = 0; k<nz; k++)
      for (int j = ny - 1; j < ny; j++)
        for (int i = 0; i < nx; i++)
        {
          int ID = nxy * k + nx * j + i;
          cellOrientationID_[ID] = 0;
          cellGrainID_[ID] = 0;
        }
  } 
}

/*------------------------------
   get_random_orientation_angle
  ------------------------------*/ 
void
CellularAutomataManager_3D::get_random_orientation_angle(double* angle)
{
  //Reference: http://mathworld.wolfram.com/SpherePointPicking.html
  // Choose U and V to be random variables on (0, 1)
  double U, V;
  U = static_cast <double> (rand()) / RAND_MAX;
  V = static_cast <double> (rand()) / RAND_MAX;
  // Then calculate the random variables of \alpha and \beta
  angle[0] = U * 2 * M_PI;
  angle[1] = acos(2 * V - 1);
  angle[2] = static_cast <double> (rand()) / RAND_MAX * 0.5 * M_PI;  
}

/*--------------------
   create_new_grain
  --------------------*/
Grain *
CellularAutomataManager_3D::create_new_grain(int cellID,
  Orientation * orientation,
  double * cellCenter, double length)
{
  Grain * grain = new Grain(this);
  grainVector_.push_back(grain);
  grain->id_ = grainVector_.size() - 1;
  grain->cell_ = cellID;
  grain->xc_ = cellCenter[0];
  grain->yc_ = cellCenter[1];
  grain->zc_ = cellCenter[2];
  grain->orientation_ = orientation;
  grain->length_ = length;

  return grain;
}

/*-------------------
   get_cell_center
  -------------------*/
void
CellularAutomataManager_3D::get_cell_center(int cellLocGhostID, double* cellCenter)
{
  int localK = cellLocGhostID / nxyLocGhost_;
  int localI = (cellLocGhostID % nxyLocGhost_) % nxLocGhost_;
  int localJ = (cellLocGhostID % nxyLocGhost_) / nxLocGhost_;
  int globalI = xStart_ + localI;
  int globalJ = yStart_ + localJ;
  int globalK = zStart_ + localK;
  cellCenter[0] = x0_ + h_ * (globalI - 0.5);
  cellCenter[1] = y0_ + h_ * (globalJ - 0.5);
  cellCenter[2] = z0_ + h_ * (globalK - 0.5);
}

/*---------------------
   get_cell_neighbors
  ---------------------*/
void
CellularAutomataManager_3D::get_cell_neighbors(int cellLocGhostID, std::vector<int> & neighborVec)
{
  int nx = nxLocGhost_;
  int ny = nyLocGhost_;
  int nxy = nxyLocGhost_;
  int numNeighbors = 26;
  int nbr[26];
  if (first_nearnest_neighbours)
  {
    numNeighbors = 6;
    nbr[0] = cellLocGhostID - 1 * nxy + 0 * nx + 0;
    nbr[1] = cellLocGhostID + 0 * nxy - 1 * nx + 0;
    nbr[2] = cellLocGhostID + 0 * nxy + 0 * nx - 1;
    nbr[3] = cellLocGhostID + 0 * nxy + 0 * nx + 1;
    nbr[4] = cellLocGhostID + 0 * nxy + 1 * nx + 0;
    nbr[5] = cellLocGhostID + 1 * nxy + 0 * nx + 0;
  }
  else
  {
    int count = 0;
    for (int k = -1; k < 2; k++)
    {
      for (int j = -1; j < 2; j++)
      {
        for (int i = -1; i < 2; i++)
        {
          if (k == 0 && j == 0 && i == 0)
            continue;
          nbr[count++] = cellLocGhostID + k*nxy + j*nx + i;
        }
      }
    }
  }
  neighborVec.clear();
  neighborVec.assign(&nbr[0], &nbr[0] + numNeighbors);
}

/* ----------------
   get_bulk_cells
  ----------------*/
void
CellularAutomataManager_3D::get_bulk_cells(std::vector<int> & cellIDs)
{
  // Return vector of all cells, excluding ghosts
  cellIDs.clear();
  cellIDs.resize(numLocalCells_);
  int count = 0;
  for (int k = 1; k < nzLocGhost_ - 1; ++k)
  {
    for (int j = 1; j < nyLocGhost_ - 1; ++j)
    {
      for (int i = 1; i < nxLocGhost_ - 1; ++i)
      {
        int id = nxyLocGhost_ * k + nxLocGhost_ * j + i;
        cellIDs[count++] = id;
      }
    }
  }
  
}

/*-------------------
   get_surface_cells
  -------------------*/
void
CellularAutomataManager_3D::get_surface_cells(std::vector<int> & cellIDs)
{
  // Loop over cells and, for surface cells, push back to vector.
  // Corner cells get included twice, and that's the intent.
  cellIDs.clear();
  for (int i = 0; i < numLocGhostCells_; ++i)
  {
    int locGhostK = i / nxyLocGhost_;
    // skip the local ghost cells
    if (locGhostK == 0 || locGhostK == nzLocGhost_ - 1)
      continue;

    int locGhostI = (i % nxyLocGhost_) % nxLocGhost_;
    // skip the local ghost cells
    if (locGhostI == 0 || locGhostI == nxLocGhost_ - 1)
      continue;
      
    int locGhostJ = (i % nxyLocGhost_) / nxLocGhost_;
    // skip the local ghost cells
    if (locGhostJ == 0 || locGhostJ == nxLocGhost_ - 1)
      continue;

      int globalI = locGhostI + xStart_;
      int globalJ = locGhostJ + yStart_;
    int globalK = locGhostK + zStart_;

    // lower and upper
    if (globalK == 1 || globalK == nz_)
      if (globalI != 0 && globalI != nx_ + 1 && globalJ != 0 && globalJ != ny_ +1)
              cellIDs.push_back(i);

    // left and right
    if (globalI == 1 || globalI == nx_)
      if (globalJ != 0 && globalJ != ny_ + 1 && globalK != 0 && globalK != nz_ + 1)
        cellIDs.push_back(i);
    
    // front and rear
    if (globalJ == 1 || globalJ == ny_)
      if (globalI != 0 && globalI != nx_ + 1 && globalK != 0 && globalK != nz_ + 1)
        cellIDs.push_back(i);
  }
}
  
/*-----------------
   cell_is_ghost
  -----------------*/
bool
CellularAutomataManager_3D::cell_is_ghost(int locGhostID)
{
  int locGhostK = locGhostID / nxyLocGhost_;
    int locGhostI = (locGhostID % nxyLocGhost_) % nxLocGhost_;
    int locGhostJ = (locGhostID % nxyLocGhost_) / nxLocGhost_;
    bool isGhost = (locGhostI == 0               ||
                    locGhostI == nxLocGhost_ - 1 ||
                    locGhostJ == 0               ||
                    locGhostJ == nyLocGhost_ - 1 ||
                    locGhostK == 0               ||
                    locGhostK == nzLocGhost_ - 1);
    return isGhost;
}

/*----------------------
   cell_is_global_ghost
  ----------------------*/
bool
CellularAutomataManager_3D::cell_is_global_ghost(int locGhostID)
{
  int globalGhostK = locGhostID / nxyLocGhost_ + zStart_;
  int globalGhostI = (locGhostID % nxyLocGhost_) % nxLocGhost_ + xStart_;
  int globalGhostJ = (locGhostID % nxyLocGhost_) / nxLocGhost_ + yStart_;

  bool isGlobalGhost = (globalGhostI == 0       ||
                        globalGhostI == nx_ + 1 ||
                        globalGhostJ == 0       ||
                        globalGhostJ == ny_ + 1 ||
                        globalGhostK == 0       ||
                        globalGhostK == nz_ + 1);
  return isGlobalGhost;
}
/*-------------------------------
   convert_id_local_to_locghost
  -------------------------------*/
void
CellularAutomataManager_3D::convert_id_local_to_locghost(int localID, int & locGhostID)
{
  int localK = localID / nxyLocal_;
  int localI = (localID % nxyLocal_) % nxLocal_;
  int localJ = (localID % nxyLocal_) / nxLocal_;
  locGhostID = (localK + 1) * nxyLocGhost_ + (localJ + 1) * nxLocGhost_ + localI + 1;
}

/*-------------------------------
   convert_id_locghost_to_local
  -------------------------------*/
void
CellularAutomataManager_3D::convert_id_locghost_to_local(int locGhostID, int & localID)
{
  int locGhostK = locGhostID / nxyLocGhost_;
  int locGhostI = (locGhostID % nxyLocGhost_) % nxLocGhost_;
  int locGhostJ = (locGhostID % nxyLocGhost_) / nxLocGhost_;
  int localI = locGhostI - 1;
  int localJ = locGhostJ - 1;
  int localK = locGhostK - 1; 
  // Make sure this cell is within local boundaries
  assert((localI >= 0) && (localI < nxLocal_));
  assert((localJ >= 0) && (localJ < nyLocal_));
  assert((localK >= 0) && (localK < nzLocal_));
  localID = localK * nxyLocal_ + localJ * nxLocal_ + localI;
}


/*-------------------------------
   convert_id_locghost_to_global
  -------------------------------*/
void
CellularAutomataManager_3D::convert_id_locghost_to_global(int locGhostID, int & globalID)
{
  int locGhostK = locGhostID / nxyLocGhost_;
  int locGhostI = (locGhostID % nxyLocGhost_) % nxLocGhost_;
  int locGhostJ = (locGhostID % nxyLocGhost_) / nxLocGhost_;
  int globalI = locGhostI + xStart_;
  int globalJ = locGhostJ + yStart_;
  int globalK = locGhostK + zStart_;
  globalID = (nxy_ + 2*(nx_ + ny_) + 4)*globalK + (nx_+2)*globalJ + globalI;
}

/*-------------------------------
   convert_id_global_to_locghost
  -------------------------------*/
void
CellularAutomataManager_3D::convert_id_global_to_locghost(int globalID, int & locGhostID)
{
  int nxyGlobal = (nx_ + 2)*(ny_ + 2);
  int globalK = globalID / nxyGlobal;
  int globalI = (globalID % nxyGlobal) % (nx_ + 2);
  int globalJ = (globalID % nxyGlobal) / (nx_ + 2);
  int locGhostI = globalI - xStart_;
  int locGhostJ = globalJ - yStart_;
  int locGhostK = globalK - zStart_;
  locGhostID = nxyLocGhost_ * locGhostK + nxLocGhost_ * locGhostJ + locGhostI;
}

/*---------------------------------
   convert_id_global_to_locghost
   meanwhile check up the position
  ---------------------------------*/
bool
CellularAutomataManager_3D::convert_id_global_to_locghost_and_check
                                     (int globalID, int & locGhostID)
{
  int nxyGlobal = (nx_ + 2)*(ny_ + 2);
  int globalK = globalID / nxyGlobal;
  int globalI = (globalID % nxyGlobal) % (nx_ + 2);
  int globalJ = (globalID % nxyGlobal) / (nx_ + 2);
  int locGhostI = globalI - xStart_;
  if (locGhostI <= 0 || locGhostI >= nxLocGhost_ - 1)
     return false;
  int locGhostJ = globalJ - yStart_;
  if (locGhostJ <= 0 || locGhostJ >= nyLocGhost_ - 1)
    return false;
  int locGhostK = globalK - zStart_;
  if (locGhostK <= 0 || locGhostK >= nzLocGhost_ - 1)
    return false;
  locGhostID = nxyLocGhost_ * locGhostK + nxLocGhost_ * locGhostJ + locGhostI;
  return true;
}

/*---------------------------------
   update_global_orientation_list
  ---------------------------------*/
void
CellularAutomataManager_3D::update_global_orientation_list(const std::vector<Orientation *> newOrientations,
                                                        const std::vector<int> updateCells)
{
  // Pack orientation angles into an array
  int numLocalOrientations = newOrientations.size();
  std::vector<double> localAngleArray(numLocalOrientations*3);
  int count = 0;
  for (std::vector<Orientation *>::const_iterator iter = newOrientations.begin();
       iter != newOrientations.end();
       ++iter)
  {
    for (int i = 0; i < 3; i++)
      localAngleArray[count++] = (*iter)->EulerAngle_[i];
  }

  // Communicate with other procs through Allgatherv
  std::vector<int> numOrientationsPerProc(numProcs_);
  MPI_Allgather(&numLocalOrientations,
                1,
                MPI_INTEGER,
                &numOrientationsPerProc[0],
                1,
                MPI_INTEGER,
                MPI_COMM_WORLD);
  std::vector<int> displs(numProcs_);
  displs[0] = 0;
  for (int i = 0; i < numProcs_; ++i)
  {
    // Since now the numOrientationsPerProc is tripled due to 3 Euler angles for each orientation
    numOrientationsPerProc[i] = numOrientationsPerProc[i] * 3;
  }
  for (int i = 1; i < numProcs_; ++i)
  {
    displs[i] = displs[i-1] + numOrientationsPerProc[i-1];
  }
  // The _3 stands for considering there are 3 Euler angles for each orientation
  int numGlobalOrientations_3 = displs[numProcs_-1] + numOrientationsPerProc[numProcs_-1];
  std::vector<double> globalAngleArray(numGlobalOrientations_3);
  if (true)//numLocalOrientations > 0)
  {
	  MPI_Allgatherv(&localAngleArray[0],
		  numLocalOrientations * 3,
		  MPI_DOUBLE,
		  &globalAngleArray[0],
		  &numOrientationsPerProc[0],
		  &displs[0],
		  MPI_DOUBLE,
		  MPI_COMM_WORLD);
  }
                 
  // Create new Orientations and add to global vector (make sure to
  // just copy the Orientations from the local proc, not create new
  // ones)
  int myDisp = displs[procID_];
  int myDisp_third = myDisp / 3;
  int numCurrentOrientations = orientationVector_.size();
  int numGlobalOrientations = numGlobalOrientations_3 / 3;
  for (int i = 0; i < numGlobalOrientations; ++i)
  {
    // From current proc
    if ((myDisp_third <= i) && (i < myDisp_third + numLocalOrientations))
    {
      newOrientations[i-myDisp_third]->id_ = numCurrentOrientations + i;
      orientationVector_.push_back(newOrientations[i-myDisp_third]);
    }
    // From other procs
    else
    {
    double EulerAngles[3];
    for (int j = 0; j < 3; j++)
      EulerAngles[j] = globalAngleArray[i * 3 + j];
      orientationVector_.push_back(new Orientation(numCurrentOrientations + i,3, EulerAngles));
    }
  }
  
  // Update cells on local proc with correct global orientation ID
  for (std::vector<int>::const_iterator iter = updateCells.begin();
       iter != updateCells.end();
       ++iter)
  {
    int cellID = *iter;
    int oldOrientationID = cellOrientationID_[cellID];
    int newOrientationID = numCurrentOrientations + myDisp_third + oldOrientationID;
    cellOrientationID_[cellID] = newOrientationID;
  }
}

/*---------------
   outputToFile
  ---------------*/
void
CellularAutomataManager_3D::outputToFile(const std::string & fileName, int * u)
{
  // Communicate all cells (including ghosted boundary cells) to proc 0 for output

  int myProcID, numProcs;
  MPI_Comm_rank(cartComm_, &myProcID);
  MPI_Comm_size(cartComm_, &numProcs);

  int gridSizeX = nxLocGhost_;
  int gridSizeY = nyLocGhost_;
  int gridSizeZ = nzLocGhost_;
  int xOffset = xStart_;
  int yOffset = yStart_;
  int zOffset = zStart_;

  // On process 0, allocate data for entire u array, as well as arrays
  // to hold grid sizes and offsets gathered from individual procs
  int * uAll;
  int *hasBeenMeltedAll;
  int *gridSizeXArray, *gridSizeYArray, *gridSizeZArray;
  int *xOffsetArray, *yOffsetArray, *zOffsetArray;
  int fullSizeX = nx_ + 2; // Include boundary cells
  int fullSizeY = ny_ + 2; // Include boundary cells
  int fullSizeZ = nz_ + 2; // Include boundary cells
  int fullSizeXY = fullSizeX * fullSizeY;
  if (myProcID == 0)
  {
    uAll = new int[fullSizeX*fullSizeY*fullSizeZ];
    hasBeenMeltedAll = new int[fullSizeX*fullSizeY*fullSizeZ];
    gridSizeXArray = new int[numProcs];
    gridSizeYArray = new int[numProcs];
    gridSizeZArray = new int[numProcs];
    xOffsetArray = new int[numProcs];
    yOffsetArray = new int[numProcs];
    zOffsetArray = new int[numProcs];
  }
  // Gather grid sizes and offsets
  MPI_Gather(&gridSizeX, 1, MPI_INTEGER,
    gridSizeXArray, 1, MPI_INTEGER, 0, cartComm_);
  MPI_Gather(&gridSizeY, 1, MPI_INTEGER,
    gridSizeYArray, 1, MPI_INTEGER, 0, cartComm_);
  MPI_Gather(&gridSizeZ, 1, MPI_INTEGER,
    gridSizeZArray, 1, MPI_INTEGER, 0, cartComm_);
  MPI_Gather(&xOffset, 1, MPI_INTEGER,
    xOffsetArray, 1, MPI_INTEGER, 0, cartComm_);
  MPI_Gather(&yOffset, 1, MPI_INTEGER,
    yOffsetArray, 1, MPI_INTEGER, 0, cartComm_);
  MPI_Gather(&zOffset, 1, MPI_INTEGER,
    zOffsetArray, 1, MPI_INTEGER, 0, cartComm_);
  // On each processor, send grid data to process 0. Use non-blocking
  // communication to avoid deadlock.
  MPI_Request request;
  MPI_Isend(u, gridSizeX*gridSizeY*gridSizeZ, MPI_INTEGER, 0, 0, cartComm_, &request);
  // On process 0, loop over processes and receive the sub-block using
  // a derived data type
  if (myProcID == 0)
  {
    for (int proc = 0; proc < numProcs; ++proc)
    {
      MPI_Datatype subblockType;
      int count = gridSizeYArray[proc] * gridSizeZArray[proc];
      int length = gridSizeXArray[proc];
      int *blocksLength = new int [count];
      int *stride = new int[count];
      int i = 0;
      for (int k = 0; k < gridSizeZArray[proc]; k++)
      {
        for (int j = 0; j < gridSizeYArray[proc]; j++)
        {
          blocksLength[i] = length;
          stride[i++] = fullSizeXY * k + fullSizeX * j;
        }
      }
      MPI_Type_indexed(count, blocksLength, stride, MPI_INTEGER, &subblockType);
      MPI_Type_commit(&subblockType);
      int * recvPointer = &uAll[zOffsetArray[proc] * fullSizeXY + 
                              yOffsetArray[proc] * fullSizeX + 
                              xOffsetArray[proc]];
      MPI_Status status;
      MPI_Recv(recvPointer, 1, subblockType, proc,
        0, cartComm_, &status);
      MPI_Type_free(&subblockType);

      delete[] blocksLength;
      delete[] stride;
    }
  }
  MPI_Wait(&request, MPI_STATUS_IGNORE);
  // for hasBeenMelted variable
  int * localHasBeenMelted;
  localHasBeenMelted = new int[numLocGhostCells_];
  int sizeHasBeenMelted = hasBeenMelted_.size();
  for (int iM = 0; iM < numLocGhostCells_; iM++) {
    if (sizeHasBeenMelted) {
      localHasBeenMelted[iM] = hasBeenMelted_[iM];
    }
    else {
      localHasBeenMelted[iM] = 1;
    }    
  }
  MPI_Isend(localHasBeenMelted, gridSizeX*gridSizeY*gridSizeZ, MPI_INTEGER, 0, 0, cartComm_, &request);
  // On process 0, loop over processes and receive the sub-block using
  // a derived data type
  if (myProcID == 0)
  {
    for (int proc = 0; proc < numProcs; ++proc)
    {
      MPI_Datatype subblockType;
      int count = gridSizeYArray[proc] * gridSizeZArray[proc];
      int length = gridSizeXArray[proc];
      int *blocksLength = new int [count];
      int *stride = new int[count];
      int i = 0;
      for (int k = 0; k < gridSizeZArray[proc]; k++)
      {
        for (int j = 0; j < gridSizeYArray[proc]; j++)
        {
          blocksLength[i] = length;
          stride[i++] = fullSizeXY * k + fullSizeX * j;
        }
      }
      MPI_Type_indexed(count, blocksLength, stride, MPI_INTEGER, &subblockType);
      MPI_Type_commit(&subblockType);
      int * recvPointer = &hasBeenMeltedAll[zOffsetArray[proc] * fullSizeXY + 
                              yOffsetArray[proc] * fullSizeX + 
                              xOffsetArray[proc]];
      MPI_Status status;
      MPI_Recv(recvPointer, 1, subblockType, proc,
        0, cartComm_, &status);
      MPI_Type_free(&subblockType);

      delete[] blocksLength;
      delete[] stride;
    }
  }
  MPI_Wait(&request, MPI_STATUS_IGNORE);
  // Output to file from proc 0
  if (myProcID == 0)
  {
    writeToFile(fileName, uAll, hasBeenMeltedAll);
  }

  // Delete arrays on proc 0
  if (myProcID == 0)
  {
    delete[] uAll;
    delete[] gridSizeXArray;
    delete[] gridSizeYArray;
    delete[] gridSizeZArray;
    delete[] xOffsetArray;
    delete[] yOffsetArray;
    delete[] zOffsetArray;
  }

}

/*--------------
   writeToFile
  --------------*/
void
CellularAutomataManager_3D::writeToFile(const std::string & fileName, const int * u, const int * meltTag)
{
  int linecount = 5;
  int count = 0;
  int NumPerLine = 10; // one line no more than this number
  int nx = nx_ + 2;
  int ny = ny_ + 2;
  int nz = nz_ + 2;
  int nxy = nx*ny;
  CafeEnv::self().caOutputP0() << "\n" << "**************************************************************";
  CafeEnv::self().caOutputP0() << "\n" << "VTKfileHeader";
  CafeEnv::self().caOutputP0() << "\n" << "**************************************************************";
  std::ofstream vtkFile;
  vtkFile.open(fileName.c_str());

  // Prepare the output file name     
  std::string fileNameOne = vtrManager_->pvdPath_ + "/finalResult.txt";    

  std::ofstream oneFile;
  std::ofstream numCellsFile;
  oneFile.open(fileNameOne.c_str());  

  vtkFile << "# vtk DataFile Version 2.0" << std::endl;
  vtkFile << "Test Data" << std::endl;
  vtkFile << "ASCII" << std::endl;
  vtkFile << "DATASET RECTILINEAR_GRID" << std::endl;

  // Random map from orientation ID to "shuffled" IDs
  // Get max value
  int maxVal = -1;
  for (int k = 1; k < nz - 1; ++k)
    for (int j = 1; j < ny - 1; ++j)
    {
      for (int i = 1; i < nx - 1; ++i)
      {
        int val = u[k*nxy + j*nx + i];
        if (val > maxVal) maxVal = val;
      }
    }

  std::vector<int> mappedID(maxVal + 1);
  if(shafferVtkID_)
  {
      for (int i = 0; i < maxVal + 1; ++i)
      {
          mappedID[i] = i;
      }
      std::random_shuffle(mappedID.begin(), mappedID.end());
  }
  // Dimensions
  vtkFile << "DIMENSIONS "
          << nx_ + 1 << " "
          << ny_ + 1 << " "
          << nz_ + 1 << std::endl;

  // Coordinates
  vtkFile << "X_COORDINATES " << nx_ + 1 << " double" << std::endl;
  linecount++;
  vtkFile << x0_<< " ";
  for (int i = 0; i < nx_; ++i)
  {
    vtkFile << x0_ + (i + 1)*h_ << " ";
    count ++;
    if(count % NumPerLine == 0) {vtkFile << std::endl; count=0; linecount++;}  // each line has NumPerLine data
  }
  if(count!=0) {vtkFile << std::endl; linecount++;}
  else
      count = 0;
  vtkFile << std::endl;
  linecount++;
  vtkFile << "Y_COORDINATES " << ny_ + 1 << " double" << std::endl;
  linecount++;
  vtkFile << y0_ << " ";
  for (int i = 0; i < ny_; ++i)
  {
    vtkFile << y0_ + (i + 1)*h_ << " ";
    count ++;
    if(count % NumPerLine == 0) {vtkFile << std::endl; count=0; linecount++;}  // each line has NumPerLine data
  }
  if(count!=0) {vtkFile << std::endl; linecount++;}
  else
      count = 0;
  vtkFile << std::endl;
  linecount++;
  vtkFile << "Z_COORDINATES " << nz_ + 1 << " double" << std::endl;  linecount++;
  vtkFile << z0_ << " ";
  for (int i = 0; i < nz_; ++i)
  {
    vtkFile << z0_ + (i + 1)*h_ << " ";
    count ++;
    if(count % NumPerLine == 0) {vtkFile << std::endl; count=0; linecount++;}  // each line has NumPerLine data
  }
  if(count!=0) {vtkFile << std::endl; linecount++;count=0;}
  else
      count = 0;
  vtkFile << std::endl;
  linecount++;

  // Cell data
  vtkFile << "CELL_DATA " << nx_*ny_*nz_ << std::endl;
  vtkFile << "SCALARS grain_id double 1" << std::endl;
  vtkFile << "LOOKUP_TABLE default" << std::endl;
  linecount += 3;
  CafeEnv::self().caOutputP0() << "\n" << "grain_id," << linecount;
  for (int k = 1; k < nz -1; ++ k){
    for (int j = 1; j < ny - 1; ++j)
    {
      for (int i = 1; i < nx - 1; ++i)
      {
        int val = u[k*nxy + j*nx + i];
        if(shafferVtkID_)
        {
            if (val >= 0)
                val = mappedID[val];
        }
        else if(shafferVtrID_)
            if (val >= 0) val = mappedGrainID[val];

        vtkFile << val << " ";        
        count ++;
        if(count % NumPerLine == 0) {vtkFile << std::endl; count=0; linecount++;}  // each line has NumPerLine data
      }
    }
  }
  if(count!=0) {vtkFile << std::endl; linecount++;}
  else
      count = 0;

  // Cell data - number of cells per grain associated with the current cell
  vtkFile << "SCALARS number_cells int 1" << std::endl;
  vtkFile << "LOOKUP_TABLE default" << std::endl;
  linecount += 2;
  CafeEnv::self().caOutputP0() << "\n" << "numCell," << linecount;
  int *numCells;
  numCells = new int[maxVal + 1];
  for (int i = 0; i < maxVal + 1; i++)
  {
    numCells[i] = 0;
  }
  // get the number of cells for each grain
  for (int k = 1; k < nz - 1; ++k) {
    for (int j = 1; j < ny - 1; ++j)
    {
      for (int i = 1; i < nx - 1; ++i)
      {
        int val = u[k*nxy + j*nx + i];
        if (val >= 0) numCells[val]++;
      }
    }
  }

  for (int k = 1; k < nz - 1; ++k) {
    for (int j = 1; j < ny - 1; ++j)
    {
      for (int i = 1; i < nx - 1; ++i)
      {
        int val = u[k*nxy + j*nx + i];
        if (val >= 0) val = numCells[val];               
        vtkFile << val << " ";

        count ++;
        if(count%10==0) {vtkFile << std::endl; count=0; linecount++;} // each line has NumPerLine data
      }
    }
  }
  if(count!=0) {vtkFile << std::endl; linecount ++;}
  else
      count = 0;

  // Cell data - melted_flag
  vtkFile << "SCALARS melted_flag int 1" << std::endl;
  vtkFile << "LOOKUP_TABLE default" << std::endl;
  linecount+=2;
  CafeEnv::self().caOutputP0() << "\n" << "meltflag," << linecount;
  for (int k = 1; k < nz - 1; ++k) {
    for (int j = 1; j < ny - 1; ++j) {
      for (int i = 1; i < nx - 1; ++i) {
        vtkFile << meltTag[k*nxy + j*nx + i] << " ";
        count ++;
        if(count%10==0) {vtkFile << std::endl; count=0; linecount++;} // each line has NumPerLine data
      }
    }
  }
  if(count!=0) {vtkFile << std::endl;linecount++;}
  else
      count = 0;
  // Cell data - Euler Angle    
  vtkFile << "VECTORS EulerAngle float" << std::endl;  
  linecount++;
  CafeEnv::self().caOutputP0() << "\n" << "EulerAngle," << linecount;
  for (int k = 1; k < nz - 1; ++k) {
    for (int j = 1; j < ny - 1; ++j)
    {
      for (int i = 1; i < nx - 1; ++i)
      {
        int val = u[k*nxy + j*nx + i];
        if (val >= 0) {
          double * EulerAngle;
          EulerAngle = (orientationVector_[val]->EulerAngle_);            
          vtkFile << EulerAngle[0] << " " << EulerAngle[1] << " " << EulerAngle[2] << " ";          

          count ++;
          if(count%10==0) {vtkFile << std::endl; count=0; linecount++;} // each line has NumPerLine data
        }
        else {          
          vtkFile << 0.0 << " " << 0.0 << " " << 0.0 << " ";
        }
      }
    }
  }
  vtkFile << std::endl;
  vtkFile.close();

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  oneFile << "X0" << "  " << "Y0" << "  " << "Z0" << "  " << " Dcell" << "  "
    << "nx_cell" << "  " << "ny_cell" << "  " << "nz_cell" << std::endl;
  oneFile << x0_ << " " << y0_ << " " << z0_ << " " << h_ << " " 
          << nx_ << " " << ny_ << " " << nz_ << std::endl;
  oneFile << "grain_ID" << "  " << "number_cells" << "  " << "melted_flag" << "  "
    << "EulerAngle1" << " " << "EulerAngle2" << " " << "EulerAngle3" << std::endl;
  for (int k = 1; k < nz - 1; ++k) {
    for (int j = 1; j < ny - 1; ++j)
    {
      for (int i = 1; i < nx - 1; ++i)
      {
        int val = u[k*nxy + j*nx + i];
        // grain ID        
		if (shafferVtkID_) { if (val >= 0) val = mappedID[val]; }
        else if(shafferVtrID_) { if (val >= 0) val = mappedGrainID[val]; }
        oneFile << val << " ";

        // number_cells    
		val = u[k*nxy + j * nx + i];
        if (val >= 0) val = numCells[val];
        oneFile << val << " ";

        // melted_flag        
        oneFile << meltTag[k*nxy + j*nx + i] << " ";

        // EulerAngle
		val = u[k*nxy + j * nx + i];
        if (val >= 0) {
          double * EulerAngle;
          EulerAngle = (orientationVector_[val]->EulerAngle_);
          oneFile << EulerAngle[0] << " " << EulerAngle[1] << " " << EulerAngle[2] << " ";
        }
        else {
          oneFile << 0.0 << " " << 0.0 << " " << 0.0 << " ";
        }
        // new line
        oneFile << std::endl;
      }
    }
  }

  oneFile.close();

  delete[] numCells;
}

/*----------------------------------------
   create_nucleation_sites_rejection_random_sampling
  ----------------------------------------*/
void
CellularAutomataManager_3D::create_nucleation_sites_rejection_random_sampling(double numberDensity,
                                                           double dTcMean,
                                                           double dTcSigma,
                                                           bool surfacePopulation)
{
  /*----------------------------------------------------
     Step 1: randomly choose unique nulceation sites and 
             determine the critical undercoolings
    ----------------------------------------------------*/
  int nxyGlobal = (nx_ + 2)*(ny_ + 2);
  int nxGlobal = nx_ + 2;

  int numCells;
  int * cellID;      // store global cell ID
  int * cellIdIndex; // store the index of cellID
  if (surfacePopulation)
  {
    numCells = (nx_*ny_ + nx_*(nz_ - 2) + (ny_ - 2)*(nz_ - 2)) * 2;
  }
  else
  {
    numCells = nx_*ny_*nz_;
  }
    
  // Global ID which counts overall ghost cells
  cellID = new int[numCells];
  int count = 0;
  for ( int k=1; k<nz_+1; k++)
    for (int j = 1; j < ny_+1; j++)
    {
      for (int i = 1; i < nx_+1; i++)
      {
        if (surfacePopulation &&
            (k != 1 && k != nz_ ) &&
            (j != 1 && j != ny_ ) &&
            (i != 1 && i != nx_ ))
              continue;

        cellID[count++] = k*nxyGlobal + j*nxGlobal + i;
      }
    }      

  cellIdIndex = new int[numCells];
  for (int i = 0; i < numCells; i++)
  {
    cellIdIndex[i] = i;
  }

  double volumePerCell;
  volumePerCell = surfacePopulation ? h_*h_ : h_*h_*h_;
  int numSites = round(numberDensity * numCells * volumePerCell) + 0.5; 
    // Add 0.5 to avoid floating point error
  if (numSites > numCells)
      throw std::runtime_error("nucleation site is too much, check the density!!");

  int * siteID;
  double * cellNucleDT;
  siteID = new int[numSites];
  cellNucleDT = new double[numSites];

  for (int i = 0; i < numSites; i++)
  {
    cellNucleDT[i] = -1;
  }

  // Set up normal distribution (note that this requires C++11)
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(dTcMean, dTcSigma);
  
  // Loop over sites to be created
  if (procID_ == 0)
  {
    for (int iSite = 0; iSite < numSites; iSite++)
    {
      // Select cell at random from list of site ids
      int t = iSite + rand() % (numCells - iSite);
      std::swap(cellIdIndex[iSite], cellIdIndex[t]);
      // Add cells to nucleationCellIDList_
      siteID[iSite] = cellIdIndex[iSite];

      // Create undercoolings from normal distribution (clip negative values)
      double dT = distribution(generator);
      dT = (dT < 0) ? 0 : dT;
      //Add undercoolings to cellNucleDT;
      cellNucleDT[iSite] = dT;
    }
  }
  
  /*--------------------------------------------------------
     Step 2: broadcast the nucleation sites and corresponding
             critical undercoolings
    --------------------------------------------------------*/
  MPI_Bcast(siteID, numSites, MPI_INTEGER, 0, cartComm_);

  MPI_Bcast(cellNucleDT, numSites, MPI_DOUBLE, 0, cartComm_);

  /*---------------------------------------------------------
     Step 3: get the nucleate cells and corresponding 
             nucleation DT that belongs to current process
    ---------------------------------------------------------*/
  int numGeneratedSites=0;
  for (int iSite = 0; iSite < numSites; iSite++)
  {
    int iCellGlobal = cellID[siteID[iSite]];
                    
    int nucleationCell;
    if(convert_id_global_to_locghost_and_check(iCellGlobal, nucleationCell))
    {
      // Add cells to nucleationCellIDList_ and undercoolings to
      // cellNucleationDT_; if a site already exists, only add it if
      // undercooling is less than existing one
      double dT = cellNucleDT[iSite];
                                 
      double dTcurr = cellNucleationDT_[nucleationCell];
      if (dTcurr < 0)
      {
        nucleationCellIDList_.push_back(nucleationCell);
        cellNucleationDT_[nucleationCell] = dT;

        // count the number of nucleation sites that belong to the local process
        numGeneratedSites++;
      }
      else if (dT < dTcurr)
      {
        cellNucleationDT_[nucleationCell] = dT;
      }

    } // end selection   
  }
  // Sort list of sites for improved performance (maybe?)
  nucleationCellIDList_.sort();

  // Assign sites' number
  if (surfacePopulation)
  {
    numSitesSurface_ = numGeneratedSites;
  }
  else
  {
    numSitesBulk_ = numGeneratedSites;
  }
  /*---------------------------------------------------------
     Step 4: clear up "new-ed" data array
    ---------------------------------------------------------*/
  if (cellID)
    delete[] cellID;

  if (cellIdIndex)
    delete[] cellIdIndex;
  
  if (siteID)
    delete[] siteID;
  
  if (cellNucleDT)
    delete[] cellNucleDT;
  
}

/*-----------------------------------------------------
   ~ checkup_current_temperature_to_updata_phase_state    
  -----------------------------------------------------*/
void
CellularAutomataManager_3D::
checkup_current_temperature_to_update_phase_state()
{
  double meltingTemp = problemPhysics_->get_melting_temperature();
  for (int iCell=0; iCell<numLocGhostCells_; iCell++)
  {
    double Tcurr=0.0;
    if (!hasBeenMelted_[iCell])      
    {
      if (mapVoxelManager_->evaluate_temperature(iCell, Tcurr))
      {
        //FIXME: do we need store temperature on the cell of CA ?
        if (Tcurr >= meltingTemp)
        {
          hasBeenMelted_[iCell] = true;
        }
      }
      else
      {
        // this is a void cell
      }
    }    
  }
}

/*----------------------------------------
     output_microstructure_information
  ----------------------------------------*/
void
CellularAutomataManager_3D::
output_microstructure_information()
{  
  // exchange boundary data before output
  commManager_->exchange_boundary_data(cellOrientationID_);
  // The reason to do the gather operation is that we declare the mappedGrainID as a array pointer
  int numSitesSurfaceAllProcesses;
  int numSitesBulkAllProcesses;
  MPI_Allreduce(&numSitesSurface_, &numSitesSurfaceAllProcesses, 1, MPI_INTEGER, MPI_SUM, cartComm_);
  MPI_Allreduce(&numSitesBulk_, &numSitesBulkAllProcesses, 1, MPI_INTEGER, MPI_SUM, cartComm_);
    
  int numCells, numOrientations, numLocalGrains, numGlobalExpectedGrains;
  numCells = numLocGhostCells_;
  numOrientations = orientationVector_.size();
  numLocalGrains = grainVector_.size();
  numGlobalExpectedGrains = numSitesSurfaceAllProcesses + numSitesBulkAllProcesses;
  numGlobalExpectedGrains += mappedGrainIDSizeFromFile_;

  // Prepare the output file name
  std::size_t dotPos = microstructureInformationFileName_output_.find_last_of(".");
  std::string extension = std::to_string(procID_) + ".txt";
  std::string fileName = microstructureInformationFileName_output_.substr(0, dotPos);
  fileName += extension;

  std::ofstream microInfoFile;
  microInfoFile.open(fileName.c_str());

  // Head line
  microInfoFile << "numCell; numOrientation; numLocalGrain; numGlobalExpectedGrain; start Point X; Y; Z; number Point nx; ny; nz; dCell" << "\n";
  microInfoFile << numCells << "  "<< numOrientations << "  "
                << numLocalGrains << "  "<< numGlobalExpectedGrains << "  "
                << x0_  << "  "<<  y0_  << "  "<<  z0_  << "  "<< nxLocGhost_ << "  "<< nyLocGhost_  << "  "<< nzLocGhost_  << "  "<< h_
                << std::endl;

  // Cell information  
  microInfoFile << "Cell informatin: ID; grainID; orientationID; hasBeenMelted; undercooling of seed" << std::endl;

  for (int iCell = 0; iCell < numCells; iCell++) {
    bool melted = false;
    bool cellSolid = true;
    if (hasBeenMelted_.size()) {
      melted = hasBeenMelted_[iCell];
    }
    microInfoFile << iCell << "  "
                  << cellGrainID_[iCell] << "  "
                  << cellOrientationID_[iCell] << "  "
                  << melted << "  "
                  << cellNucleationDT_[iCell];
    microInfoFile << std::endl;
  }
  
  // Orientation information
  microInfoFile << "Orientation informatin: ID; EularAngle" << std::endl;
  for (int iOrien = 0; iOrien < numOrientations; iOrien++) {
    microInfoFile << orientationVector_[iOrien]->id_ << "  "; // ID
    for (int i = 0; i < 3; i++) { // Eular angles
      microInfoFile << orientationVector_[iOrien]->EulerAngle_[i] << "  ";
    }
    microInfoFile << std::endl; // change the line
  }

  // Grain information
  microInfoFile << "Grain information: ID; cellID; grainCenter; orientationID" << std::endl;
  for (int iGrain = 0; iGrain < numLocalGrains; iGrain++) {
    microInfoFile << grainVector_[iGrain]->id_ << "  ";
    microInfoFile << grainVector_[iGrain]->cell_ << "  ";
    
    microInfoFile << grainVector_[iGrain]->xc_ << "  ";
    microInfoFile << grainVector_[iGrain]->yc_ << "  ";
    microInfoFile << grainVector_[iGrain]->zc_ << "  ";

    //microInfoFile << grainVector_[iGrain]->length_ << "  ";
    int orientationID = grainVector_[iGrain]->orientation_->id_;
    microInfoFile << orientationID << "   " << std::endl;;
  }

  microInfoFile.close();
}

/*--------------------------------------
   load_microstructure_information
  --------------------------------------*/
void
CellularAutomataManager_3D::
load_microstructure_information()
{
  // do nothing here
}

/*--------------------------------------
   random_microstructure_information
  --------------------------------------*/
void
CellularAutomataManager_3D::
random_microstructure_information()
{
    std::vector<Orientation *> newOrientationVector;
    std::vector<int> nucleatedCells;

    for (int iCell = 0; iCell < numLocGhostCells_; iCell++){
        if(cell_is_ghost(iCell))
        {
            cellGrainID_[iCell] = -1;
            cellOrientationID_[iCell] = -1;
            continue;
        }

        // Nucleate with random orientation angle
        double angle[3];
        get_random_orientation_angle(angle);

        int iOrientation = newOrientationVector.size();
        Orientation * orientation = new Orientation(iOrientation, nDim_, angle);
        newOrientationVector.push_back(orientation);
        // Create new grain at cell center
        double cellCenter[3];
        get_cell_center(iCell, cellCenter);
        Grain * grain = create_new_grain(iCell, orientation, cellCenter,0);
        int grainID = grain->id_;
        int orientationID = orientation->id_;

        // Add to cell, and keep track of cells so we can correct the orientation IDs
        cellGrainID_[iCell] = grainID;
        cellOrientationID_[iCell] = orientationID;
        nucleatedCells.push_back(iCell);
    }
    // Add local orientations to global list and update IDs stored on cells
    update_global_orientation_list(newOrientationVector, nucleatedCells);
    commManager_->exchange_boundary_data(cellOrientationID_);
    int grainNum = orientationVector_.size();
    std::cout << "Maxium initial grain ID is: " << grainNum << std::endl;
    if (shafferVtrID_) {
        int numGlobalExpectedGrains = (nx_+2)*(ny_+2)*(nz_+2);
        mappedGrainID = new int[(numGlobalExpectedGrains)];
        if (procID_ == 0)
        {
            // Random map from orientation ID to "shuffled" IDs
            for (int i = 0; i < numGlobalExpectedGrains; i++)
            {
                mappedGrainID[i] = i;
            }
            if(shafferVtrID_)
                std::random_shuffle(mappedGrainID, mappedGrainID + numGlobalExpectedGrains);
        }
        MPI_Bcast(mappedGrainID, numGlobalExpectedGrains, MPI_INTEGER, 0, cartComm_);
    }
    else
        mappedGrainID = NULL;
}

/*-----------------------------------
output_last_step_result
-----------------------------------*/
void
CellularAutomataManager_3D::output_last_step_result()
{
  vtrManager_->execute(0, time_);
}


/*------------------------------------------
get_amount_of_final_orientations
------------------------------------------*/
void
CellularAutomataManager_3D::
get_amount_of_final_orientations(int &numberOrientation)
{
    // work for both normal and remelting problems
    int orienSize = orientationVector_.size();
    int * localOrienIDRecord;
    localOrienIDRecord = new int[orienSize];
    int * globalOrienIDRecord;
    globalOrienIDRecord = new int[orienSize];

    for (int i = 0; i < orienSize; i++) {
        localOrienIDRecord[i] = 0;
        globalOrienIDRecord[i] = 0;
    }

    for (int i = 0; i < numLocGhostCells_; i++)
    {
        int ID = cellOrientationID_[i];
        if (ID > -1) {
            localOrienIDRecord[ID] = 1;
        }
    }

    MPI_Allreduce(localOrienIDRecord, globalOrienIDRecord, orienSize, MPI_INTEGER, MPI_SUM, cartComm_);

    numberOrientation = 0;
    for (int i = 0; i < orienSize; i++) {
        if (globalOrienIDRecord[i] > 0) {
            numberOrientation++;
        }
    }

    delete[] localOrienIDRecord;
    delete[] globalOrienIDRecord;
}

/*---------------------------------
get_amount_of_liquid_cells
---------------------------------*/
bool
CellularAutomataManager_3D::
get_amount_of_liquid_cells(int &numberLiquidCell)
{
    numberLiquidCell = 0;
    int localNLC = 0;
    for (int i = 0; i < numLocGhostCells_; i++)
    {
        if (!cell_is_ghost(i)) {
            if (cellOrientationID_[i] < 0) {
                localNLC++;
            }
        }
    }

    MPI_Allreduce(&localNLC, &numberLiquidCell, 1, MPI_INTEGER, MPI_SUM, cartComm_);
    if (numberLiquidCell == 0) {
        solidificationDone_ = true;
        return false;
    }
    else
        return true;
}
