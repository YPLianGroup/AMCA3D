// General includes
#include <yaml-cpp/yaml.h>
#include <iostream>
#include <fstream>
#include <random>
#include <mpi.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <stdlib.h>
#include <algorithm>
// Local includes
#include "CellularAutomataManager.h"
#include "Grain.h"
#include "Orientation.h"
#include "ParallelCommManager.h"
#include "ProblemPhysics.h"
#include "Timer.h"
#include "CafeParsing.h"
#include "MapVoxelManager.h"

#define M_PI 3.14159265358979323846 
/*=======================================================================
   Class Definition
     CellularAutomataManager - manage flowchat of the calculation
     (abstract base class)
=======================================================================*/

/*------------------
   Constructor
  ------------------*/
CellularAutomataManager::CellularAutomataManager()
  : maxTimestep_(1.0),
  cellOrientationID_(NULL),
  cellGrainID_(NULL),
  cellNucleationDT_(NULL),
  cellTemperature_(NULL),
  cooling_(NULL),
  mapVoxelManager_(NULL),
  shafferVtrID_(false),
  shafferVtkID_(false),
  time_(0.0),
  iStep_(0),
  timeStepFactor_(0.2),
  numSitesSurface_(0),
  numSitesBulk_(0),
  first_nearnest_neighbours(false),
  ns_(0),//2.5e2),
  deltaTs_max_(0.5),
  deltaTs_sigma_(0.1),
  nv_(0),//5.5e1),
  deltaTv_max_(6.0),
  deltaTv_sigma_(0.1),
  setInitialStatedAsMelted_(true),
  cellGrainLegnthCutoff_(1.0e8),
  outputMicrostructureInfo_(false),
  loadMicrostructureInfo_(false),
  randomMicrostructureInfo_(false),
  setNucleationSite_(false),
  activateExistingGrain_(true),
  solidificationDone_(false),
  output_nucleation_seeds_(false),
  load_nucleation_seeds_(false),
  mappedGrainIDSizeFromFile_(0),
  couplingTimeStepFactor_(0.2),
  capture_mushy_only_(false)
{
}

/*------------------
   Destructor
  ------------------*/
CellularAutomataManager::~CellularAutomataManager()
{
  // Delete cell arrays
  if (cellOrientationID_) 
        delete[] cellOrientationID_;
  if (cellGrainID_) 
        delete[] cellGrainID_;
  if (cellNucleationDT_) 
        delete[] cellNucleationDT_;

  // Delete grain and orientation pointers
  for (int i = 0; i < grainVector_.size(); ++i){
    delete grainVector_[i];
  }

  for (int i = 0; i < orientationVector_.size(); ++i){
    delete orientationVector_[i];
  }

  if (mappedGrainID){
    delete[] mappedGrainID;
  }

  if (cellTemperature_) {
    delete[] cellTemperature_;
  }
  if (cooling_) {
    delete[] cooling_;
  }
}

/*----------------------
   set_timer_pointer
  ----------------------*/
void
CellularAutomataManager::set_timer_pointer(Timer *timer)
{
  CAtimer_ = timer;
}

/*-------------------------------
   setup_time_step_for_coupling
  -------------------------------*/
void
CellularAutomataManager::setup_time_step_for_coupling(double dT)
{
  dT_ = dT;
}
/*------------------
   setup_integration
  ------------------*/
void
CellularAutomataManager::setup_integration(double maxDt)
{
  maxTimestep_ = maxDt;
}

/*------------------
   time_integration
  ------------------*/
void
CellularAutomataManager::time_integration(double finalTime,
  int maxSteps)
{  
  while (time_ < finalTime && !solidificationDone_)
  {
    CAtimer_->start_timing_CA_simulation();
    //Nucleate new grains based on current undercooling and then compute time step
    calculate_time_step();
    single_integration_step();
  }

  if (outputMicrostructureInfo_)
    output_microstructure_information();
}

/*-----------------------------------------
   create_nucleation_sites_random_sampling
  -----------------------------------------*/
void
CellularAutomataManager::create_nucleation_sites(double numberDensity,
  double dTcMean,
  double dTcSigma,
  bool surfacePopulation)
{
  // Compute number of sites to create
  std::vector<int> cellIDs;
  if (surfacePopulation)
    get_surface_cells(cellIDs);
  else
    get_bulk_cells(cellIDs);

  double volumePerCell;
  if (nDim_ == 2)
    volumePerCell = surfacePopulation ? h_ : h_*h_;
  else
    volumePerCell = surfacePopulation ? h_*h_ : h_*h_*h_;
  int numCells = cellIDs.size();
  int numSites = round(volumePerCell * numberDensity * numCells ) + 0.5; // Add 0.5 to avoid floating point error

  // Set up normal distribution (note that this requires C++11)
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(dTcMean, dTcSigma);

  int numGeneratedSites = 0;
  // Loop over sites to be created
  for (int iSite = 0; iSite < numSites; ++iSite)
  {
    // Select cell at random from list of site ids
    int i = rand() % numCells;
    int nucleationCell = cellIDs[i];

    // Create undercoolings from normal distribution (clip negative values)
    double dT = distribution(generator);
    dT = (dT < 0) ? 0 : dT;

    // Add cells to nucleationCellIDList_ and undercoolings to
    // cellNucleationDT_; if a site already exists, only add it if
    // undercooling is less than existing one
    double dTcurr = cellNucleationDT_[nucleationCell];
    if (dTcurr < 0)
    {
      nucleationCellIDList_.push_back(nucleationCell);
      cellNucleationDT_[nucleationCell] = dT;
      numGeneratedSites++;
    }
    else if (dT < dTcurr)
    {
      cellNucleationDT_[nucleationCell] = dT;
    }
  }

  // Sort list of sites for improved performance (maybe?)
  nucleationCellIDList_.sort();

  // Assign sites' number
  if (surfacePopulation)
    numSitesSurface_ = numGeneratedSites;
  else
    numSitesBulk_ = numGeneratedSites;
}

/*--------------------------------
   nucleate_at_nucleation_sites
  --------------------------------*/
void
CellularAutomataManager::nucleate_at_nucleation_sites(double time)
{
  std::vector<Orientation *> newOrientationVector;
  std::vector<int> nucleatedCells;

  double meltingTemperature = problemPhysics_->get_melting_temperature();
  // Loop over nucleation sites
  std::list<int>::iterator iter = nucleationCellIDList_.begin();
  while (iter != nucleationCellIDList_.end())
  {
    int iCell = *iter;
    // If cell has already been captured, don't nucleate here; just remove nucleation site from list
    bool isLiquid = (cellGrainID_[iCell] == -1);
    if (isLiquid)
    {
      double dTcrit = cellNucleationDT_[iCell];
      double temperature;
        if (mapVoxelManager_)
        {
            // interpolate temperature from FEM result
            if (mapVoxelManager_->evaluate_temperature(iCell, temperature))
            {
                if (!hasBeenMelted_[iCell])
                {
                    // reset the temperature above melting temperature to skip the following work
                    // can not use continue here
                    temperature = meltingTemperature + dTcrit + 1.0;
                }
            }
            else
            {
                // no element tag right now corresponding to no material within the cell
                // reset the temperature above melting temperature to skip the following work
                // can not use continue here
                temperature = meltingTemperature + dTcrit + 1.0;
            }
        }
        else
        {
            temperature = problemPhysics_->compute_cell_temperature(iCell, time);
        }

      double dT = meltingTemperature - temperature;
      if (dT > dTcrit)
      {
        // Nucleate with random orientation angle
        double angle[3];
        get_random_orientation_angle(angle);
        int iOrientation = newOrientationVector.size();
        Orientation * orientation = new Orientation(iOrientation, nDim_, angle);
        newOrientationVector.push_back(orientation);

        // Create new grain at cell center
        double cellCenter[3];
        get_cell_center(iCell, cellCenter);
        Grain * grain = create_new_grain(iCell, orientation, cellCenter);
        int grainID = grain->id_;
        int orientationID = orientation->id_;

        // Add to cell, and keep track of cells so we can correct the orientation IDs
        cellGrainID_[iCell] = grainID;
        cellOrientationID_[iCell] = orientationID;
        nucleatedCells.push_back(iCell);

        isLiquid = false;
      }
    }
    // Remove cell from nucleation list if it is no longer liquid
    if (isLiquid)
      ++iter;
    else
    {
      cellNucleationDT_[iCell] = 0;
      iter = nucleationCellIDList_.erase(iter);
    }
  }
  // Add local orientations to global list and update IDs stored on cells
  update_global_orientation_list(newOrientationVector, nucleatedCells);
}

/*------------------
   Allocate
  ------------------*/
void
CellularAutomataManager::allocate_arrays()
{
  // Allocate for grid including ghosts
  cellOrientationID_ = new int[numLocGhostCells_];
  cellGrainID_ = new int[numLocGhostCells_];
  cellNucleationDT_ = new double[numLocGhostCells_];
  grainVector_.reserve(numLocalCells_);
}

/*-------------------
  initialize_cells
  -------------------*/
void
CellularAutomataManager::initialize_cells()
{
  // Initial grain and orientation id's of -1 indicate liquid cells
  // everywhere
  for (int i = 0; i < numLocGhostCells_; ++i)
  {
    cellOrientationID_[i] = -1;
    cellGrainID_[i] = -1;
    cellNucleationDT_[i] = -1.0;
  }
}

/*---------------------
   create_new_grain
  --------------------*/
Grain *
CellularAutomataManager::create_new_grain(int cellID,
  Orientation * orientation,
  double * cellCenter, double length)
{
  Grain * grain = new Grain(this);
  grainVector_.push_back(grain);
  grain->id_ = grainVector_.size() - 1;
  grain->cell_ = cellID;
  grain->xc_ = cellCenter[0];
  grain->yc_ = cellCenter[1];
  grain->orientation_ = orientation;
  grain->length_ = length;

  return grain;
}

/*------------------
   cell_is_captured
  ------------------*/
bool
CellularAutomataManager::cell_is_captured(const Grain * grain, int iCell)
{
  bool isCaptured = false;
  double length = grain->length_;
  if (length == 0) return isCaptured;
  double xcgrain = grain->xc_;
  double ycgrain = grain->yc_;
  double angle = grain->orientation_->angle_;
  double xc, yc;

  // Compute "grain coordinates" of cell center
  double cellCenter[2];
  get_cell_center(iCell, cellCenter);
  xc = cellCenter[0];
  yc = cellCenter[1];
  double egx[2] = { cos(angle), sin(angle) };
  double egy[2] = { -sin(angle), cos(angle) };
  double cx = (egx[0] * (xc - xcgrain) + egx[1] * (yc - ycgrain)) / length;
  double cy = (egy[0] * (xc - xcgrain) + egy[1] * (yc - ycgrain)) / length;
  if (fabs(cx) + fabs(cy) <= 1.0)
  {
    isCaptured = true;
  }
  return isCaptured;
}

/*--------------
   capture_cell
  --------------*/
void
CellularAutomataManager::capture_cell(double xcgrain,
  double ycgrain,
  double length,
  int orientationID,
  int iCell)
{
  Orientation * orientation = orientationVector_[orientationID];
  double angle = orientation->angle_;

  // Compute "grain coordinates" of cell center
  double xc, yc;
  double cellCenter[2];
  get_cell_center(iCell, cellCenter);
  xc = cellCenter[0];
  yc = cellCenter[1];
  double egx[2] = { cos(angle), sin(angle) };
  double egy[2] = { -sin(angle), cos(angle) };
  double cx = (egx[0] * (xc - xcgrain) + egx[1] * (yc - ycgrain)) / length;
  double cy = (egy[0] * (xc - xcgrain) + egy[1] * (yc - ycgrain)) / length;
  assert(fabs(cx) + fabs(cy) <= 1); // Cell center should be inside square

                    // Determine which quadrant
  int quad = 1;
  if (cx >= 0 && cy >= 0) quad = 1;
  else if (cx >= 0 && cy <  0) quad = 4;
  else if (cx <  0 && cy >= 0) quad = 2;
  else if (cx <  0 && cy <  0) quad = 3;

  // Find endpoints of "capture line"
  double angle1 = angle + (quad - 1) * M_PI / 2;
  double angle2 = angle + (quad)* M_PI / 2;
  double xp1 = xcgrain + length * cos(angle1);
  double yp1 = ycgrain + length * sin(angle1);
  double xp2 = xcgrain + length * cos(angle2);
  double yp2 = ycgrain + length * sin(angle2);

  // Project cell center onto capture line
  double xvp = xp2 - xp1;
  double yvp = yp2 - yp1;
  double xvc = xc - xp1;
  double yvc = yc - yp1;
  double normvp2 = xvp*xvp + yvp*yvp;
  double vpdotvc = xvp*xvc + yvp*yvc;
  double xc0 = xp1 + xvp*vpdotvc / normvp2;
  double yc0 = yp1 + yvp*vpdotvc / normvp2;

  // Compute line segments
  double L1 = sqrt((xc0 - xp1)*(xc0 - xp1) + (yc0 - yp1)*(yc0 - yp1));
  double L2 = sqrt((xc0 - xp2)*(xc0 - xp2) + (yc0 - yp2)*(yc0 - yp2));

  // Compute new grain size
  double newGrainLength = std::min(L1 / sqrt(2.), h_) + std::min(L2 / sqrt(2.), h_);

  // Compute closest corner
  double xcc = xp1;
  double ycc = yp1;
  if (L2 < L1)
  {
    xcc = xp2;
    ycc = yp2;
  }

  // Compute new grain center
  double norm_xcc_gc = sqrt((xcc - xcgrain)*(xcc - xcgrain) + (ycc - ycgrain)*(ycc - ycgrain));
  double xcnew = xcgrain + (length - newGrainLength) * (xcc - xcgrain) / norm_xcc_gc;
  double ycnew = ycgrain + (length - newGrainLength) * (ycc - ycgrain) / norm_xcc_gc;

  // Create new grain with same orientation and assign it to this cell
  double newGrainCenter[2];
  newGrainCenter[0] = xcnew;
  newGrainCenter[1] = ycnew;
  Grain * newGrain = create_new_grain(iCell, orientation, newGrainCenter, newGrainLength);
  cellGrainID_[iCell] = newGrain->id_;
  cellOrientationID_[iCell] = orientationID;
}

/*---------------------
   output_orientations
  ---------------------*/
void
CellularAutomataManager::output_orientations(const std::string & fileName)
{
  outputToFile(fileName, cellOrientationID_);
}

/*-------------------------------------------------------------
       mapped_grain_ID_array
  -------------------------------------------------------------*/
void
CellularAutomataManager::
global_mapped_grain_ID(double numberDensitySurface,
                          double numberDensityBulk)
{
  int numCellsSurface = (nDim_==2) ? (nx_ + ny_ - 2)*2
                                     :(nx_*ny_+nx_*(nz_-2)+(ny_-2)*(nz_-2))*2;
  
  int numCellsBulk = (nDim_ == 2) ? nx_*ny_
                                   : nx_*ny_*nz_;

  double areaPerCell = h_*h_;
  double volumePerCell = h_*h_*h_;
  int numSitesSurface = round(numberDensitySurface * numCellsSurface * areaPerCell) + 0.5;
  int numSitesBulk = round(numberDensityBulk * numCellsBulk * volumePerCell) + 0.5;
  int numGrains = numSitesSurface + numSitesBulk;
  int numSitesSurfaceAllProcesses;
  int numSitesBulkAllProcesses;
  MPI_Allreduce(&numSitesSurface_, &numSitesSurfaceAllProcesses, 1, MPI_INTEGER, MPI_SUM, cartComm_);
  MPI_Allreduce(&numSitesBulk_, &numSitesBulkAllProcesses, 1, MPI_INTEGER, MPI_SUM, cartComm_);
  numGrains = numSitesSurfaceAllProcesses + numSitesBulkAllProcesses;
  mappedGrainID= new int[numGrains];
  if (procID_ == 0)
  {
    // Random map from orientation ID to "shuffled" IDs
    for (int i = 0; i < numGrains; i++)
    {
      mappedGrainID[i] = i;
    }
    if(shafferVtrID_)
        std::random_shuffle(mappedGrainID, mappedGrainID + numGrains);
    
  }

  MPI_Bcast(mappedGrainID, numGrains, MPI_INTEGER, 0, cartComm_);
  if (!shafferVtrID_)
      mappedGrainID = NULL;

  // Brief Review
  CafeEnv::self().caOutputP0() << "\n" << "Nucleation Sites Review" << "\n";
  CafeEnv::self().caOutputP0() << "=============================" << "\n";
  CafeEnv::self().caOutputP0() << "Surface Nucleation Sites:  "<< "\n";
  CafeEnv::self().caOutputP0() << "                 should be = " << numSitesSurface << "\n";

  CafeEnv::self().caOutputP0() << "Bulk Nucleation Sites:  " << "\n";
  CafeEnv::self().caOutputP0() << "                 should be = " << numSitesBulk << "\n";

  
}

/*-------------------------------------------
   ~ setup_map_voxel_manager_and_accessory
  -------------------------------------------*/
void
CellularAutomataManager::
setup_map_voxel_manager_and_accessory(MapVoxelManager* p)
{
  mapVoxelManager_ = p;
}

/*----------------------------------------
    output_microstructure_information   
  ----------------------------------------*/
void
CellularAutomataManager::
output_microstructure_information()
{
  // do nothing here
}

/*--------------------------------------
    load_microstructure_information  
  --------------------------------------*/
void
CellularAutomataManager::
load_microstructure_information()
{
  // do nothing here
}

/*-----------------------------------
    output_last_step_result
  -----------------------------------*/
void 
CellularAutomataManager::output_last_step_result()
{
  // nothing here
}
