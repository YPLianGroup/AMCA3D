#include <yaml-cpp/yaml.h>
#include <fstream>
#include <mpi.h>
#include <random>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <algorithm>

#include "CafeParsing.h"
#include "CellularAutomataManager_3DRemelting.h"
#include "VtrManager.h"
#include "Grain.h"
#include "Orientation.h"
#include "MapVoxelManager.h"
#include "ProblemPhysics.h"
#include "Envelope.h"
#include "Octahedron.h"
#include "ParallelCommManager_3D.h"
#include "Timer.h"
#include "CafeMacro.h"
#define m_PI 3.14159265358979323846;
/*---------------------
   construction
  ---------------------*/
CellularAutomataManager_3DRemelting::CellularAutomataManager_3DRemelting()
  :CellularAutomataManager_3D()
{
  nullMatTemp_ = -1.0;
}

/*---------------------
   destruction
  ---------------------*/
CellularAutomataManager_3DRemelting::~CellularAutomataManager_3DRemelting()
{
  // nothing for now
}
/*--------------------------------------
   setup_several_data_member_pointer
  --------------------------------------*/
void
CellularAutomataManager_3DRemelting::setup_several_data_member_pointer()
{
  commManager_ = new ParallelCommManager_3D(this);
  vtrManager_ = new VtrManager(this);
}
/*---------------------
   load
  ---------------------*/
void
CellularAutomataManager_3DRemelting::load(const YAML::Node& node)
{
  CellularAutomataManager_3D::load(node);

  // solution option
  const YAML::Node * solutionOptions = node.FindValue("solution_options");
  if (solutionOptions) {
      const YAML::Node *options_node = expect_sequence(*solutionOptions, "options", true);
      if (options_node) {
          for (size_t ioption = 0; ioption < options_node->size(); ioption++) {
              const YAML::Node &option_node = (*options_node)[ioption];
              if (expect_map(option_node, "criterion_to_activate_grain", true)) {
                  const YAML::Node &grainOption = *option_node.FindValue("criterion_to_activate_grain");

                  get_required(grainOption, "using_nucleation_seed", setNucleationSite_);
                  CafeEnv::self().caOutputP0() << "Criteria for nucleation:" << "\n";
                  CafeEnv::self().caOutputP0() << "                         using nucleation seeds "
                                               << setNucleationSite_ << "\n";
                  CafeEnv::self().caOutputP0() << "                         activate existing grains "
                                               << activateExistingGrain_ << "\n";
              }
          }
      }
  }
}

/*---------------------
    initialization
  ---------------------*/
void
CellularAutomataManager_3DRemelting::initialize()
{
  nxy_ = nx_ * ny_;

  // set up CA netwrok
  setup_cells();
  // set up the length cutoff
  cellGrainLegnthCutoff_ *= h_;

  // allocate temperature_ and associations
  cellTemperature_ = new double[numLocGhostCells_];
  cooling_ = new bool[numLocGhostCells_];
  hasBeenMelted_.resize(numLocGhostCells_);
  cellLiquidNeigh.resize(numLocGhostCells_);
  for (int iCell = 0; iCell < numLocGhostCells_; iCell++) {
    cellTemperature_[iCell] = nullMatTemp_;
    cooling_[iCell] = false;
    cellLiquidNeigh[iCell] = false;

    if (setInitialStatedAsMelted_)
      hasBeenMelted_[iCell] = true;
    else
      hasBeenMelted_[iCell] = false;
  }

  // initialize vtr manager
  vtrManager_->initialize();

  // initialize microstructure information
  int numberExistingGrains = 0;
  if(randomMicrostructureInfo_) {
      random_microstructure_information();
  }
  else if (loadMicrostructureInfo_ && activateExistingGrain_) {
    load_microstructure_information(numberExistingGrains);
    // retrieve orientation ID for cells without material in previous simulation
    retrieve_grain_information_before_simulation();
  }

  if (setNucleationSite_){    
    // Create nucleation sites on surface
    create_nucleation_sites(ns_, deltaTs_max_, deltaTs_sigma_, true);

    // Create nucleation sites in volume
    create_nucleation_sites(nv_, deltaTv_max_, deltaTv_sigma_, false);

    // generate global mapped grain ID
    if (shafferVtrID_)
        global_mapped_grain_ID_plus(ns_, nv_, numberExistingGrains);
  }

  std::cout<<"initializing done!"<<std::endl<<"======================="<<std::endl;
  problemPhysics_->initialize();
}

/*-------------------------------
   Integrate up to a given time
  -------------------------------*/
void 
CellularAutomataManager_3DRemelting::
time_integration(double finalTime, int maxSteps)
{
  while (time_ < finalTime && iStep_ < maxSteps)
  {
    // update cell grain information
    checkup_current_temperature_to_update_phase_state();

    // calcute time step from active grains
    calculate_time_step();

    // grow grain and capture neighbors
    single_integration_step();
  }
 
  if (outputMicrostructureInfo_) {
    output_microstructure_information();
  }
}

/*----------------------------------
    load_microstructure_information
  ----------------------------------*/
void
CellularAutomataManager_3DRemelting::
load_microstructure_information(int &numGlobalExpectedGrains)
{
  numSitesBulk_ = 0;
  numSitesSurface_ = 0;

  int numCells, numOrientations, numLocalGrains;

  // Prepare the output file name
  std::size_t dotPos = microstructureInformationFileName_load_.find_last_of(".");
  std::string extension = std::to_string(procID_) + ".txt";
  std::string fileName = microstructureInformationFileName_load_.substr(0, dotPos);
  fileName += extension;

  std::fstream inputStream;
  inputStream.open(fileName.c_str(), std::ios::in);

  if (!inputStream) {
    std::cout << "fail to open the file" << fileName << std::endl;
    throw std::runtime_error("fail to open the specified file");
  }
  double X, Y, Z, nx, ny, nz, dcell;
  std::string line;
  // get the header
  {
    std::getline(inputStream, line); // skip the description line
    std::getline(inputStream, line); // get the value of the header
    std::istringstream iss(line);
    iss >> numCells >> numOrientations >> numLocalGrains >> numGlobalExpectedGrains >> X >> Y >> Z >> nx >> ny >> nz >> dcell;
    mappedGrainIDSizeFromFile_ = numGlobalExpectedGrains; // to save the value of numGloabalExpectedGrains
    // generate grain id mapped array
    if (!setNucleationSite_) {
      mappedGrainID = new int[numGlobalExpectedGrains];
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
      if (!shafferVtrID_)
          mappedGrainID = NULL;
    }    
  }
  // check the current cell number is equal to microstructure information file or not
  if((nxLocGhost_*nyLocGhost_*nzLocGhost_)!=numCells){
      std::cout << "Error : microstructure information file is not fit the current domain size" << fileName << std::endl;
      std::cout << "nxLocGhost: " << nxLocGhost_ << "  nyLocGhost: " << nyLocGhost_ << "  nzLocGhost: " <<  nzLocGhost_ << std::endl;
      std::cout << "numCells from microstructure information file: " << numCells << std::endl;
      std::cout << "nx: " << nx << "  ny: " << ny << "  nz: " << nz << std::endl;
      throw std::runtime_error("Error : microstructure information file is not fit the current domain size");
  }
  else
      std::cout << "Start loading cell information from files" << std::endl;


  // get cell information
  int count = 0;
  int cellID, cellGrainID, cellOrientationID;
  bool melted;
  double nucleationDt;
  std::getline(inputStream, line); // skip the description line

  while (count < numCells) {
    std::getline(inputStream, line);
    if (!line.size())
      continue; // this is an empty line

    std::istringstream iss(line);
    iss >> cellID >> cellGrainID >> cellOrientationID >> melted >> nucleationDt;
    //iss >> cellID >> cellGrainID >> cellOrientationID;
    cellGrainID_[cellID] = cellGrainID;
    cellOrientationID_[cellID] = cellOrientationID;
    hasBeenMelted_[cellID] = melted;
    if (load_nucleation_seeds_) {
      if (nucleationDt > -1) {
        cellNucleationDT_[cellID] = nucleationDt;
        nucleationCellIDList_.push_back(cellID);
      }
    }
    count++;
  }
  
  // get orientation information
  std::cout << "Start loading orientation information from files" << std::endl;
  count = 0;  
  std::getline(inputStream, line); // skip the description line
  while (count < numOrientations) {
    std::getline(inputStream, line);
    if (!line.size())
      continue; // this is an empty line

    std::istringstream iss(line);
    Orientation * orien = new Orientation;
    iss >> orien->id_
        >> orien->EulerAngle_[0]
        >> orien->EulerAngle_[1]
        >> orien->EulerAngle_[2];
    count++;
    orientationVector_.push_back(orien);
  }

  // get grain information
  std::cout << "Start loading grain information from files" << std::endl;
  count = 0;
  std::getline(inputStream, line); //skip the description line
  int orientationID;
  while (count < numLocalGrains) {
    std::getline(inputStream, line);
    if (!line.size())
      continue; // this is an empty line

    Grain *grain = new Grain(this, false);
    std::istringstream iss(line);
    iss >> grain->id_ >> grain->cell_
      >> grain->xc_ >> grain->yc_ >> grain->zc_
      >> orientationID;
    grain->length_ = 0.0;

    // Reset the grain center as the center of the cell
    if (true) { // FIEME: will add an boolean variable to control this
      double cellCenter[3];
      get_cell_center(grain->cell_, cellCenter);
      grain->xc_ = cellCenter[0];
      grain->yc_ = cellCenter[1];
      grain->zc_ = cellCenter[2];
    }

    grain->inactivate();
    grainVector_.push_back(grain);
    grainVector_[count]->orientation_ = orientationVector_[orientationID];
    count++;
  }

  inputStream.close();
  inputStream.clear();
  std::cout << "Loading microstructure files done!" << std::endl;
}

/*----------------------------------------------------
    retrieve_grain_information_before_simulation
  ----------------------------------------------------*/
void
CellularAutomataManager_3DRemelting::
retrieve_grain_information_before_simulation()
{
  for (int iCell = 0; iCell < numLocGhostCells_; iCell++)
  {
    if (cell_is_global_ghost((iCell)))
        cellOrientationID_[iCell] = -1;
    if (cell_is_ghost(iCell)) { // this is not the local process's business
      continue;
    }
    if (cellOrientationID_[iCell] > -1) {
      continue;
    }
    // get cell grain id
    int grainID = cellGrainID_[iCell];
    // retrieve the orientation id
    int orientationID = grainVector_[grainID]->orientation_->id_;
    cellOrientationID_[iCell] = orientationID;
  }
}

/*---------------------
   update_cell_grain
  ---------------------*/
void
CellularAutomataManager_3DRemelting::
checkup_current_temperature_to_update_phase_state()
{
  CAtimer_->start_timing_CA_temp();
  double meltingTemp = problemPhysics_->get_melting_temperature();
  double solidusTemp = problemPhysics_->get_solidus_temperature();
  int nx = nxLocGhost_;
  int nxy = nxyLocGhost_;
  for (int iCell = 0; iCell < numLocGhostCells_; iCell++)
      cellLiquidNeigh[iCell] = false; //reset cellLiquidNeigh vector
  for (int i = 0; i < nxLocGhost_; i++) {
      for (int j = 0; j < nyLocGhost_; j++) {
          for (int k = 0; k < nzLocGhost_ + 1; k++) {
              if (i < 0 || j < 0 || k < 0 || i >= nxLocGhost_ || j >= nyLocGhost_ || k >= nzLocGhost_)
                  continue;
              int iCell = k * nxy + j * nx + i;
              if (numLocGhostCells_ < iCell)
              {
                  throw std::runtime_error("error: numLocGhostCells_ < iCell!");
              }  
              double temperature = nullMatTemp_;
              bool withMat = false;
              withMat = mapVoxelManager_->evaluate_temperature(iCell, temperature); // get it back later

              if (iStep_ > 0) {
                  double tempRate = 410;
                  if (cellTemperature_[iCell] > 0)
                      tempRate = (temperature - cellTemperature_[iCell]) / dT_;
                  if (temperature <= 0)
                      tempRate = 0;
                  if (withMat) {
                      cooling_[iCell] = (tempRate < 0);
                  }
                  else
                      cooling_[iCell] = false;
              }

              cellTemperature_[iCell] = temperature;

              if (cell_is_global_ghost(iCell)) {
                  continue;
              }

              if (withMat && temperature < meltingTemp) { // no change
                  continue;
              }

              // only for cells with melted material or null material    
              cellOrientationID_[iCell] = -1;

              if (withMat) {
                  hasBeenMelted_[iCell] = true;
              }
              else { // since no material, it makes sense to reset hasBeenMelted_ be false.
                  hasBeenMelted_[iCell] = false;
              }

              int grainID = cellGrainID_[iCell];
              if (grainID > -1) { // Remark: keep the orientation pointer
                  grainVector_[grainID]->inactivate();
              }
          }
      }
  }
  // for (int iCell = 0; iCell < numLocGhostCells_; iCell++) 
  CAtimer_->stop_timing_CA_temp();
}

void
CellularAutomataManager_3DRemelting::checkup_current_temperature()
{
    // read the solidified track grain information, so it is no need to checkup cell phase
    for (int iCell = 0; iCell < numLocGhostCells_; iCell++)
    {
        if(cellOrientationID_[iCell] == -1){
            cellTemperature_[iCell] = nullMatTemp_;
            continue;
        }
        double temperature = nullMatTemp_;
        bool withMat = mapVoxelManager_->evaluate_temperature(iCell, temperature); // get it back later
        cellTemperature_[iCell] = temperature;
    }
}

/*---------------------
   check whether a cell has liquid cell Neighbor
  ---------------------*/
void
CellularAutomataManager_3DRemelting::
checkup_cell_neigh_to_growth(double time)
{
  CAtimer_->start_timing_CA_cell_active();
  double meltingTemp = problemPhysics_->get_melting_temperature();
  double solidusTemp = problemPhysics_->get_solidus_temperature();
  for (int iCell = 0; iCell < numLocGhostCells_; iCell++) 
  {    
    cellLiquidNeigh[iCell] = false;
    double temperature = nullMatTemp_;
    bool withMat = mapVoxelManager_->evaluate_temperature(iCell, temperature); // get it back later
    if(cellOrientationID_[iCell] == -1 && withMat && time >=maxTimestep_)//pass nonMat cell and liquid cell
    //time >= maxTimestep_ to avoid setting the cell in new added powder to melting 
    {
      hasBeenMelted_[iCell] = true;
      continue;
    }
    if(temperature>solidusTemp){
        if(cooling_[iCell])
            continue;
    }
    if(!withMat)
    {
        hasBeenMelted_[iCell] = false;
        continue;
    }
    std::vector<int> neighborVec;
    get_cell_neighbors(iCell, neighborVec);
    bool hasNonCapturedNeighbors = false;
    for (int i = 0; i < neighborVec.size(); ++i)
    {
      int iNeigh = neighborVec[i];
      double temperature = cellTemperature_[iNeigh];
      
      // there is a grain on the cell already
      if (cellOrientationID_[iNeigh] > -1) {
        continue;
      }

      // the cell does not contain material
      if (!(temperature > nullMatTemp_)) {
        continue;
      }
      //can't capture overhot cell
      if ((temperature > meltingTemp)) {
        continue;
      }
      hasNonCapturedNeighbors = true;
    } // end if (nonGrain nonVoid mushyTemperature)
    cellLiquidNeigh[iCell] = hasNonCapturedNeighbors;
  } // end neighbors loop
  CAtimer_->stop_timing_CA_cell_active();
}

/*---------------------
   check whether a cell has liquid cell Neighbor though a liquid cell
   loop every liquid cell's neighbor to set a growth label, this may speed up program because liquid cell is more less than solid
  ---------------------*/
void
CellularAutomataManager_3DRemelting::
checkup_liquid_cell_neigh_to_growth(double time)
{
  double meltingTemp = problemPhysics_->get_melting_temperature();
  double solidusTemp = problemPhysics_->get_solidus_temperature();
  int nx = nxLocGhost_;
  int nxy = nxyLocGhost_;
  for (int iCell = 0; iCell < numLocGhostCells_; iCell++)
      cellLiquidNeigh[iCell] = false; //reset cellLiquidNeigh vector
    for (int i = 1; i < nxLocGhost_-1; i++) {
        for (int j = 1; j < nyLocGhost_-1; j++) {
            for (int k = 1; k < nzLocGhost_ - 1; k++) {
              if (i<0 || j<0 || k<0 || i>=nxLocGhost_ || j>=nyLocGhost_|| k>= nzLocGhost_)
                  continue;
              int iCell = k * nxy + j * nx + i;
              if (iCell > numLocGhostCells_) {  
                  throw std::runtime_error("error!!iCell>=numLocGhostCells");
              }
                  
              double temperature = cellTemperature_[iCell];
              bool withMat = true;
              if(temperature<nullMatTemp_+0.01)
                  withMat = false;
              if(cellOrientationID_[iCell] == -1 && withMat && time >=maxTimestep_)//mat flow
                  //time >= maxTimestep_ to avoid setting the cell in new added powder to melting
              {
                  hasBeenMelted_[iCell] = true;
              }
              if(!withMat)
              {
                  hasBeenMelted_[iCell] = false;
                  continue;
              }
              if(cellOrientationID_[iCell] != -1 || cell_is_ghost(iCell))//pass global ghost cell, solid cell and non mat cell
                  // here, local ghost cell might be liquid cell, so it is important to loop the neighb of local ghost cell and then label as growth state.
                  continue;
              if(temperature >= meltingTemp) // cooling_[iCell]==false ||
                  continue;
              if(capture_mushy_only_ && temperature<solidusTemp)
                  continue;
              std::vector<int> neighborVec;
              get_cell_neighbors(iCell, neighborVec);
              for (int i = 0; i < neighborVec.size(); ++i)
              {
                  int iNeigh = neighborVec[i];
                  if(false)    // iNeigh might out of array's bounds, skipping this to avoid memery error
                      // Known bug, FixMe latter: iNeigh might not be a neighb of a ghost cell, e.g., iNeigh is localted in the other side of cubic face.      GhostCell|...CubicDomin...|iNeigh
                      continue;
                  //double temperature = cellTemperature_[iNeigh];
                  // if a liquid cell has a solid cell neighbor
                  if (cellOrientationID_[iNeigh] > -1 && !cell_is_global_ghost(iCell)) {
                      cellLiquidNeigh[iNeigh] = true;
                  }
              }
          }
      }
  }
}

void
CellularAutomataManager_3DRemelting::single_integration_step() {
    /*------------------------------------------------
       Grow grains if it has liquid neighbor
      ------------------------------------------------*/
    CAtimer_->start_timing_capture();
    int nx = nxLocGhost_;
    int nxy = nxyLocGhost_;
    for (int i = 1; i < nxLocGhost_-1; i++) {
        for (int j = 1; j < nyLocGhost_-1; j++) {
            for (int k = 1; k < nzLocGhost_ - 1; k++) {
                if (i < 0 || j < 0 || k < 0 || i >= nxLocGhost_ || j >= nyLocGhost_ || k >= nzLocGhost_)
                    continue;
                int iCell = k * nxy + j * nx + i;
                if (iCell > numLocGhostCells_) {
                    throw std::runtime_error("error!!iCell>=numLocGhostCells");
                }
                if (!cellLiquidNeigh[iCell] || cell_is_ghost(iCell))
                    continue;
                double solidusTemp = problemPhysics_->get_solidus_temperature();
                double meltingTemp = problemPhysics_->get_melting_temperature();
                if (cellGrainID_[iCell] < 0)
                    continue;
                Grain *grain = grainVector_[cellGrainID_[iCell]];

                double vel;
                double T = cellTemperature_[iCell];
                if (T < solidusTemp) {
                    vel = problemPhysics_->compute_growth_velocity_with_temperature_input(solidusTemp + 0.01);
                } else {
                    vel = problemPhysics_->compute_growth_velocity_with_temperature_input(T);
                }
                grain->grow(vel, dT_);
                std::vector<int> neighborVec;
                get_cell_neighbors(iCell, neighborVec);

                for (int i = 0; i < neighborVec.size(); ++i) {
                    int iNeigh = neighborVec[i];
                    double temperature = cellTemperature_[iNeigh];

                    // there is a grain on the cell already
                    if (cellOrientationID_[iNeigh] > -1) {
                        continue;
                    }
                    // the cell is overhot
                    if (temperature > meltingTemp) {
                        continue;
                    }
                    if (capture_mushy_only_ && temperature < solidusTemp)
                        continue;
                    // the cell does not contain material
                    if (!(temperature > nullMatTemp_)) {
                        continue;
                    }
                    // no need to capture global ghost
                    if (cell_is_global_ghost(iNeigh)) {
                        continue;
                    }
                    /*----------------------------------------
                       Capture neighbors of the active grains
                      ----------------------------------------*/
                    double cellCenter[3];
                    get_cell_center(iNeigh, cellCenter);
                    bool isCaptured = envelope_->cell_is_captured(grain, cellCenter);
                    if (isCaptured) {
                        bool cellIsGhost = cell_is_ghost(iNeigh);
                        double grainCenter[3];
                        grainCenter[0] = grain->xc_;
                        grainCenter[1] = grain->yc_;
                        grainCenter[2] = grain->zc_;
                        if (cellIsGhost) {
                            commManager_->pack_for_comm(grainCenter, grain->length_,
                                                        grain->orientation_->id_, iNeigh);
                        } else {
                            capture_cell(grainCenter, grain->length_, grain->orientation_->id_,
                                         iNeigh, false);
                        }
                    }
                } // end neighbors loop
            }
        }
    }

  CAtimer_->stop_timing_capture();
  // Parallel communication, and handle captures across proc boundaries
  CAtimer_->start_timing_crossProcCapture();
  commManager_->cross_proc_captures();
  CAtimer_->stop_timing_crossProcCapture();

  // Update time
  time_ += dT_;
  iStep_++;

  CAtimer_->stop_timing_CA_simulation();

  // call vtrManager
  vtrManager_->execute(iStep_, time_);

  /*-----------------------------------------------------------
              Output current time step log
    -----------------------------------------------------------*/
  CafeEnv::self().caOutputP0() << "*******************"
                               << "\n";
  CafeEnv::self().caOutputP0() << "Cellular Automaton Solver Part: " << "\n";
  CafeEnv::self().caOutputP0() << "time step count = " << iStep_
                               << " current time = " << time_ << "\n"
                               << " time step = " << dT_ << std::endl;

}
/*---------------------
   calculate_time_step
  ---------------------*/
double 
CellularAutomataManager_3DRemelting::calculate_time_step()
{ 
  CAtimer_->start_timing_nucleation();
  if (setNucleationSite_){
    // Nucleate new grains based on current undercooling
    nucleate_at_nucleation_sites(time_);
  }
  CAtimer_->stop_timing_nucleation();
  /*--------------------------
     calculate time step
     Get the Minimuim temp
    --------------------------*/
  int numGrains = grainVector_.size();

  double maxVelocity = 1e-10;
  double solidusTemp = problemPhysics_->get_solidus_temperature();
  double minT=problemPhysics_->get_melting_temperature();
  int maximumID = -1;
  int nx = nxLocGhost_;
  int nxy = nxyLocGhost_;
    for (int i = 1; i < nxLocGhost_-1; i++) {
        for (int j = 1; j < nyLocGhost_-1; j++) {
            for (int k = 1; k < nzLocGhost_ - 1; k++) {
              int iCell = k * nxy + j * nx + i;
              if (!cellLiquidNeigh[iCell] || cell_is_ghost(iCell))
                  continue;
              if (cellGrainID_[iCell] < 0)
                  continue;
              // only work on active cell
              double T = cellTemperature_[iCell];
              Grain* grain = grainVector_[cellGrainID_[iCell]];
              if (grain->length_ < cellGrainLegnthCutoff_) {
                  // Set the limit of vel to avoid a too small time step
                  // vel = problemPhysics_->compute_growth_velocity_with_temperature_input(T);
              }
              else { // Reset grain length and velocity since they do not make sense here!
                  grain->length_ = h_;
              }
              if (T <= solidusTemp) {
                  //skip low temperature cell to avoid extreme small time step
                  continue;
              }
              if (minT > T) {
                  minT = T;
                  maximumID = iCell;
              }
          }
      }
  }
  maxVelocity = (problemPhysics_->compute_growth_velocity_with_temperature_input(minT) > 0 ?
          problemPhysics_->compute_growth_velocity_with_temperature_input(minT) : 0);
  if (maximumID > -1) {
    std::vector<int> neighborVec;
    std::cout << "Time step is dertermined by Grain" << cellGrainID_[maximumID]
      << " Temp=" << cellTemperature_[maximumID] << " OrientationID="
      << cellOrientationID_[maximumID] << std::endl << "maxVelocity=" << maxVelocity << std::endl;

    get_cell_neighbors(maximumID, neighborVec);
    for (int i = 0; i < neighborVec.size(); i++)
    {
      int iNeigh = neighborVec[i];
      std::cout << "NeighborCell= " << iNeigh << "  Temp=" << cellTemperature_[iNeigh]
        << " OrientationID=" << cellOrientationID_[iNeigh]
        << " GrainID=" << cellGrainID_[iNeigh] << std::endl;
    }
  }

  double dtLoc = std::min(timeStepFactor_*h_ / maxVelocity, maxTimestep_);
  MPI_Allreduce(&dtLoc, &dT_, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  return dT_;
}

void
CellularAutomataManager_3DRemelting::
maintain_active_grains_on_interface_cells()
{
  // looping on liquid cells and active cells on the interface
  double meltingTemp = problemPhysics_->get_melting_temperature();
  double solidTemp = problemPhysics_->get_solidus_temperature();
  for (int iCell = 0; iCell < numLocGhostCells_; iCell++){
    // Inactivate all grain
    int grainID = cellGrainID_[iCell];
    if (grainID > -1){
       grainVector_[grainID]->inactivate();
    }
  }
  for (int iCell = 0; iCell < numLocGhostCells_; iCell++)
  {
    if (cell_is_ghost(iCell)) { // this is not the local process's business
      continue;
    }

    double iCellTemp = cellTemperature_[iCell];
    // updata the cell hasBeenMelted state
    if(nullMatTemp_==iCellTemp){
      hasBeenMelted_[iCell]=false;
      continue;
    }
    else{
      // cell has been melted if cell is liquid
      cellOrientationID_[iCell] < 0 ? hasBeenMelted_[iCell]=true:0;
    }
    // This following condition must be after the above one
    if (!hasBeenMelted_[iCell]) { // the cell has never been melted
      continue;
    }
    // skip overhot cell
    if(iCellTemp > meltingTemp || solidTemp > iCellTemp)
      continue;
    /*---------------------------------------------------
    activate all solid neighbourhood grain of liquid cell's
    ---------------------------------------------------*/

    if (cellOrientationID_[iCell] < 0)
    {// if cell is liquid, then activate all of its solid neighbourhood grain
      std::vector<int> neighborVec;
      get_cell_neighbors(iCell, neighborVec);
      int activeCellID = -1;
      for (int i = 0; i < neighborVec.size(); i++)
      {
        int iNeigh = neighborVec[i];
        int grainID = cellGrainID_[iNeigh];
        if (grainID > -1 && (!cell_is_ghost(iNeigh)) && grainVector_[grainID]->is_active()) {
          // Because there is an active grain, no need to activate one!
          activeCellID = -1;
          continue;
        }

        // active a solid grain with a liquid neighb
        if (cellOrientationID_[iNeigh] > -1) {
          double neighbCellTemp = cellTemperature_[iNeigh];
          activeCellID = iNeigh;
          if ((activeCellID > -1) && (!cell_is_global_ghost(activeCellID))) {// activate one
              int activeGrainID = cellGrainID_[activeCellID];
              if (cooling_[activeCellID] && activeGrainID>-1) {
                  grainVector_[activeGrainID]->activate();
              }
          }
          activeCellID = -1;
        }
      } // end for i
    } // end for if
  } // end for iCell
}

/*------------------------------------------------------------------
   Capture cell by creating a new grain with the appropriate
   orientation
   capture_cell
  -----------------------------------------------------------------*/
void
CellularAutomataManager_3DRemelting::
capture_cell(double* grainCenter, double length, int orientationID, int iCell, bool resetCoordMatrix)
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

  // Get information of the new active grain with same orientation and assign it to this cell
  double newGrainCenter[3];
  for (int i = 0; i < 3; i++)
    newGrainCenter[i] = newGrainInfo[i] + grainCenter[i];

  double newGrainLength = newGrainInfo[3];  
  if (cellGrainID_[iCell] < 0) { // create new grain if it does not exist
    Grain * newGrain = create_new_grain(iCell, orientation, newGrainCenter, newGrainLength);
    cellGrainID_[iCell] = newGrain->id_;
  }
  else { // update the old grain information
    Grain *oldGrain = grainVector_[cellGrainID_[iCell]];
    oldGrain->activate();
    oldGrain->length_ = newGrainLength;
    oldGrain->xc_ = newGrainCenter[0];
    oldGrain->yc_ = newGrainCenter[1];
    oldGrain->zc_ = newGrainCenter[2];
    oldGrain->orientation_ = orientation;
  }
  
  cellOrientationID_[iCell] = orientationID;
}

/*--------------------
   create_new_grain
  --------------------*/
Grain *
CellularAutomataManager_3DRemelting::create_new_grain(int cellID,
  Orientation * orientation,
  double * cellCenter, double length)
{
  Grain * grain = new Grain(this, false);
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

/*-------------------------------------------------------------
    global_mapped_grain_ID_plus
-------------------------------------------------------------*/
void
CellularAutomataManager_3DRemelting::
global_mapped_grain_ID_plus(double numberDensitySurface,
  double numberDensityBulk, int numberExistingGrains)
{
  int numCellsSurface = (nDim_ == 2) ? (nx_ + ny_ - 2) * 2
                        : (nx_*ny_ + nx_*(nz_ - 2) + (ny_ - 2)*(nz_ - 2)) * 2;

  int numCellsBulk = (nDim_ == 2) ? nx_*ny_ : nx_*ny_*nz_;

  double areaPerCell = h_*h_;
  double volumePerCell = h_*h_*h_;
  int numSitesSurface = round(numberDensitySurface * numCellsSurface * areaPerCell) + 0.5;
  int numSitesBulk = round(numberDensityBulk * numCellsBulk * volumePerCell) + 0.5;

  int numSitesSurfaceAllProcesses;
  int numSitesBulkAllProcesses;

  MPI_Allreduce(&numSitesSurface_, &numSitesSurfaceAllProcesses, 1, MPI_INTEGER, MPI_SUM, cartComm_);

  MPI_Allreduce(&numSitesBulk_, &numSitesBulkAllProcesses, 1, MPI_INTEGER, MPI_SUM, cartComm_);

  int numGrains = numSitesSurface + numSitesBulk + numberExistingGrains;
  if (numGrains < 0) {
      throw std::runtime_error("numGrains must be positive!\n");
  }
  mappedGrainID = new int[numGrains];
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
  CafeEnv::self().caOutputP0() << "Surface Nucleation Sites:  " << "\n";
  CafeEnv::self().caOutputP0() << "                 should be = " << numSitesSurface << "\n";

  CafeEnv::self().caOutputP0() << "Bulk Nucleation Sites:  " << "\n";
  CafeEnv::self().caOutputP0() << "                 should be = " << numSitesBulk << "\n";

  CafeEnv::self().caOutputP0() << "Existing Unique Grains:  " << "\n";
  CafeEnv::self().caOutputP0() << "                  may be " << numberExistingGrains << "\n";

  CafeEnv::self().caOutputP0() << "Final Unique Grains:  " << "\n";
  CafeEnv::self().caOutputP0() << "                  may be " << numGrains << "\n";

}

/*--------------------------------
nucleate_at_nucleation_sites
--------------------------------*/
void
CellularAutomataManager_3DRemelting::nucleate_at_nucleation_sites(double time)
{
  std::vector<Orientation *> newOrientationVector;
  std::vector<int> nucleatedCells;

  double meltingTemperature = problemPhysics_->get_melting_temperature();
  double minTempForNucleation = problemPhysics_->get_minimum_temperature_for_nucleation();
  // Loop over nucleation sites
  std::list<int>::iterator iter = nucleationCellIDList_.begin();
  while (iter != nucleationCellIDList_.end())
  {
    int iCell = *iter;

    // If cell has already been captured, don't nucleate here; just remove nucleation site from list
    bool isLiquid = (cellOrientationID_[iCell] == -1);
    bool beenMelted = hasBeenMelted_[iCell];
    double temperature = cellTemperature_[iCell];
    bool withMushyMat = (temperature > minTempForNucleation);
    if (cooling_[iCell] && isLiquid && beenMelted && withMushyMat)
    {
      double dTcrit = cellNucleationDT_[iCell];            
      double dT = meltingTemperature - cellTemperature_[iCell];

      if (dT > dTcrit)
      {
        // Nucleate with random orientation angle
        double angle[3];
        get_random_orientation_angle(angle);
        int iOrientation = newOrientationVector.size();
        Orientation * orientation = new Orientation(iOrientation, nDim_, angle);
        newOrientationVector.push_back(orientation);
        
        double cellCenter[3];
        get_cell_center(iCell, cellCenter);

        int grainID = cellGrainID_[iCell];
        if (grainID == -1) {// Create new grain at cell center          
          Grain * grain = create_new_grain(iCell, orientation, cellCenter, 0);
          int grainID = grain->id_;
          cellGrainID_[iCell] = grainID;
        }
        else { // reset the grain information          
          Grain * grain = grainVector_[grainID];
          grain->activate();
          grain->length_ = 0.0;
          grain->xc_ = cellCenter[0];
          grain->yc_ = cellCenter[1];
          grain->zc_ = cellCenter[2];
          grain->orientation_ = orientation;
        }
        
        int orientationID = orientation->id_;
        // Add to cell, and keep track of cells so we can correct the orientation IDs        
        cellOrientationID_[iCell] = orientationID;
        nucleatedCells.push_back(iCell);
      }
    }
    // Keep the cell in nucleation list
    ++iter;    
  }
  // Add local orientations to global list and update IDs stored on cells
  update_global_orientation_list(newOrientationVector, nucleatedCells);
}
