#include <iostream>
#include <fstream>
#include <random>

// General includes
#include "mpi.h"
#include "math.h"
#include "assert.h"
#include "time.h"
#include "stdlib.h"
#include <algorithm>
// Local includes
#include "CellularAutomataManager_2D.h"
#include "Grain.h"
#include "Orientation.h"
#include "ParallelCommManager_2D.h"
#include "ProblemPhysics.h"
#include "CafeParsing.h"

#define M_PI 3.14159265358979323846 
/// ------------------
//  Constructor
// ------------------
CellularAutomataManager_2D::CellularAutomataManager_2D()
  : CellularAutomataManager()
{
  // Parallel setup
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs_);
  MPI_Comm_rank(MPI_COMM_WORLD, &procID_);
  commManager_ = new ParallelCommManager_2D(this);

  // Random seed
  srand(procID_ * 10000);
}

// ------------------
//  Destructor
// ------------------
CellularAutomataManager_2D::~CellularAutomataManager_2D()
{
  // Delete cell arrays
  if (cellOrientationID_) delete cellOrientationID_;
  if (cellGrainID_) delete cellGrainID_;
  if (cellNucleationDT_) delete cellNucleationDT_;

  // Delete grain and orientation pointers
  for (int i = 0; i < grainVector_.size(); ++i)
  {
    delete grainVector_[i];
  }
  for (int i = 0; i < orientationVector_.size(); ++i)
  {
    delete orientationVector_[i];
  }

}
//-----------------
// load
// ----------------
void
CellularAutomataManager_2D::load(const YAML::Node& node)
{
  double original_point[2];
  double lateral_size[2];
  double dCell[2];
  int    numbCell[2];

  // get domain size
  const YAML::Node * domain_node = node.FindValue("domain");
  if (domain_node)
  {
    get_required(*domain_node, "original_point", original_point);
    get_required(*domain_node, "lateral_sizes", lateral_size);
  }
  else
    throw std::runtime_error("parser error: realm-domain");

  // get discretization info
  const YAML::Node * discretization_node = node.FindValue("discretization");
  if (discretization_node)
  {
    get_required(*discretization_node, "cell_sizes", dCell);
    get_required(*discretization_node, "number_cells", numbCell);
  }
  else
    throw std::runtime_error("parser error: realm-discretization");

  // original point coordinates
  x0_ = original_point[0];
  y0_ = original_point[1];

  // mesh information
  h_ = dCell[0];
  nx_ = numbCell[0];
  ny_ = numbCell[1];

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
        std::string type = "Gaussian";
        get_if_present(*surface_node, "type", type, type);
        get_required(*surface_node, "site_density", ns_);
        get_required(*surface_node, "mean", deltaTs_max_);
        get_required(*surface_node, "standard_deviation", deltaTs_sigma_);
      }
      const YAML::Node * bulk_node = nucleation_node.FindValue("bulk");
      if (bulk_node)
      {
        std::string type = "Gaussian";
        get_if_present(*bulk_node, "type", type, type);
        get_required(*bulk_node, "site_density", nv_);
        get_required(*bulk_node, "mean", deltaTv_max_);
        get_required(*bulk_node, "standard_deviation", deltaTv_sigma_);
      }
    }
  }
  else
    throw std::runtime_error("parser error: realm-nucleation_rules");
}

/*--------------------------------------
  setup_several_data_member_pointer
  --------------------------------------*/
void 
CellularAutomataManager_2D::setup_several_data_member_pointer()
{
  // nothing for now
}

// ----------------
// initialize
// ----------------
void
CellularAutomataManager_2D::initialize()
{
    // convert to 2D geometry (see R&G 1993, Table 1), 1/mm
    ns_ = sqrt(4 * ns_ / M_PI);
    nv_ = pow(nv_*sqrt(6. / M_PI), 2. / 3.);

    setup_cells();
      // Create nucleation sites on surface
      create_nucleation_sites(ns_, deltaTs_max_, deltaTs_sigma_, true);

      // Create nucleation sites in volume
      create_nucleation_sites(nv_, deltaTv_max_, deltaTv_sigma_, false);

      // generate global mapped grain ID
      global_mapped_grain_ID(ns_, nv_);
}
//-----------------
// Set up problem
//-----------------
void CellularAutomataManager_2D::setup_problem(int nDim, double* geometry, double h, int* meshInfo)
{
  // geometry info
  nDim_ = nDim;
  x0_ = geometry[0];
  y0_ = geometry[1];

  // mesh information
  h_ = h;
  nx_ = meshInfo[0];
  ny_ = meshInfo[1];
}
// ---------------
//  setup_cells
// ---------------
void
CellularAutomataManager_2D::setup_cells()
{
  // Set up Cartesian topology
  int dims[2];
  commManager_->setup_cartesian_comm();
  cartComm_ = commManager_->get_cart_comm();

  // Compute size of local grid
  get_local_size(cartComm_, 0, nx_, nxLocal_, xStart_);
  get_local_size(cartComm_, 1, ny_, nyLocal_, yStart_);
  numLocalCells_ = nxLocal_ * nyLocal_;
  nxLocGhost_ = nxLocal_ + 2;
  nyLocGhost_ = nyLocal_ + 2;
  numLocGhostCells_ = nxLocGhost_ * nyLocGhost_;

  // Allocate data
  allocate_arrays();

  // Initialize cells
  initialize_cells();
  initialize_boundaries();
}

//-------------------------
// calculate_time_step
//-------------------------
double
CellularAutomataManager_2D::calculate_time_step()
{
  // Nucleate new grains based on current undercooling
  nucleate_at_nucleation_sites(time_);

  // Loop over grains to compute timestep
  double maxVelocity = 1e-10;
  for (std::list<Grain *>::iterator iter = activeGrainList_.begin();
    iter != activeGrainList_.end();
    ++iter)
  {
    Grain * grain = *iter;
    double vel = problemPhysics_->compute_growth_velocity(grain, time_);
    if (vel > maxVelocity)
    {
      maxVelocity = vel;
    }
  }
  double dtLoc = std::min(timeStepFactor_*h_ / maxVelocity, maxTimestep_);
  MPI_Allreduce(&dtLoc, &dT_, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  
  return dT_;
}

// -------------------------
//  single_integration_step
// -------------------------
void
CellularAutomataManager_2D::single_integration_step()
{
  // Grow grains
  for (std::list<Grain *>::iterator iter = activeGrainList_.begin();
    iter != activeGrainList_.end();
    ++iter)
  {
    Grain * grain = *iter;
    double vel = problemPhysics_->compute_growth_velocity(grain, time_);
    grain->grow(vel, dT_);
  }

  // Check for capture of neighbors
  std::list<Grain *>::iterator iter = activeGrainList_.begin();

  while (iter != activeGrainList_.end())
  {
    Grain * grain = *iter;
    int iCell = grain->cell_;
    std::vector<int> neighborVec;
    get_cell_neighbors(iCell, neighborVec);
    bool hasLiquidNeighbors = false;
    for (int i = 0; i < neighborVec.size(); ++i)
    {
      int iNeigh = neighborVec[i];
      bool isLiquid = (cellGrainID_[iNeigh] == -1);
      if (isLiquid)
      {
        bool isCaptured = cell_is_captured(grain, iNeigh);

        // This is where we'll need something to handle parallelism, eventually
        if (isCaptured)
        {
          bool cellIsGhost = cell_is_ghost(iNeigh);
          if (cellIsGhost)
          {
            commManager_->pack_for_comm(grain->xc_,
              grain->yc_,
              grain->length_,
              grain->orientation_->id_,
              iNeigh);
          }
          else
          {
            capture_cell(grain->xc_,
              grain->yc_,
              grain->length_,
              grain->orientation_->id_,
              iNeigh);
          }
        }
        else
        {
          hasLiquidNeighbors = true;
        }

      } // end if (isLiquid)

    } // end neighbor loop

      // Deactivate grain if it has no more liquid neighbors; otherwise
      // update iterator
    if (!hasLiquidNeighbors)
    {
      grain->inactivate();
      iter = activeGrainList_.erase(iter);
    }
    else
    {
      ++iter;
    }

  }

  // Parallel communication, and handle captures across proc boundaries
  commManager_->cross_proc_captures();

  // Update time
  time_ += dT_;

  if (procID_ == 0)
  {
    std::cout << "time = " << time_ << std::endl;
  }

}


// ------------------
//  get_local_size
// ------------------
void
CellularAutomataManager_2D::get_local_size(MPI_Comm cartComm, int dim, int n, int & nLocal, int & offset)
{
  // Get my process coordinates
  int cartDims[2];
  int periods[2];
  int cartCoords[2];
  MPI_Cart_get(cartComm, 2, cartDims, periods, cartCoords);

  // Compute number of local interior nodes, and offset
  // global_I = local_i + offset
  nLocal = n / cartDims[dim];
  offset = cartCoords[dim] * nLocal;
  // Account for non-uniform distribution on last proc: offset + nlocal = n
  if (cartCoords[dim] == cartDims[dim] - 1)
  {
    nLocal = n - offset;
  }
}

// ------------------------
//  initialize_boundaries
// ------------------------
void
CellularAutomataManager_2D::initialize_boundaries()
{
  // Set all boundaries to have 0 grain and orientation id (so that
  // they are treated as non-liquid by their neighbors)

  // Get my process coordinates
  int cartDims[2];
  int periods[2];
  int cartCoords[2];
  MPI_Cart_get(cartComm_, 2, cartDims, periods, cartCoords);

  int nx = nxLocGhost_;
  int ny = nyLocGhost_;

  // Lower BC
  if (cartCoords[1] == 0)
  {
    int j = 0;
    for (int i = 0; i < nx; ++i)
    {
      cellOrientationID_[nx * j + i] = 0;
      cellGrainID_[nx * j + i] = 0;
    }
  }
  // Upper BC
  if (cartCoords[1] == cartDims[1] - 1)
  {
    int j = ny - 1;
    for (int i = 0; i < nx; ++i)
    {
      cellOrientationID_[nx * j + i] = 0;
      cellGrainID_[nx * j + i] = 0;
    }
  }
  // Left BC
  if (cartCoords[0] == 0)
  {
    int i = 0;
    for (int j = 0; j < ny; ++j)
    {
      cellOrientationID_[nx * j + i] = 0;
      cellGrainID_[nx * j + i] = 0;
    }
  }
  // Right BC
  if (cartCoords[0] == cartDims[0] - 1)
  {
    int i = nx - 1;
    for (int j = 0; j < ny; ++j)
    {
      cellOrientationID_[nx * j + i] = 0;
      cellGrainID_[nx * j + i] = 0;
    }
  }

}


// ------------------------------
//  get_random_orientation_angle
// ------------------------------
void
CellularAutomataManager_2D::get_random_orientation_angle(double * angle)
{
  // Compute random angle between 0 and pi/2
  angle[0] = static_cast <double> (rand()) / RAND_MAX * M_PI / 2.0;
}


// --------------------
//  create_new_grain
// --------------------
Grain *
CellularAutomataManager_2D::create_new_grain(int cellID,
  Orientation * orientation,
  double * newGrainCenter, double length)
{
  double xc, yc;
  xc = newGrainCenter[0];
  yc = newGrainCenter[1];
  Grain * grain = new Grain(this);
  grainVector_.push_back(grain);
  grain->id_ = grainVector_.size() - 1;
  grain->cell_ = cellID;
  grain->xc_ = xc;
  grain->yc_ = yc;
  grain->orientation_ = orientation;
  grain->length_ = length;

  return grain;
}

// ------------------
//  cell_is_captured
// ------------------
bool
CellularAutomataManager_2D::cell_is_captured(const Grain * grain, int iCell)
{
  bool isCaptured = false;
  double length = grain->length_;
  if (length == 0) return isCaptured;
  double xcgrain = grain->xc_;
  double ycgrain = grain->yc_;
  double angle = grain->orientation_->angle_;
  double xc, yc;
  double cellCenter[2];

  // Compute "grain coordinates" of cell center
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

// --------------
//  capture_cell
// --------------
void
CellularAutomataManager_2D::capture_cell(double xcgrain,
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

// -------------------
//  get_cell_center
// -------------------
void
CellularAutomataManager_2D::get_cell_center(int cellLocGhostID, double * cellCenter)
{
  int localI = cellLocGhostID % nxLocGhost_;
  int localJ = cellLocGhostID / nxLocGhost_;
  int globalI = xStart_ + localI;
  int globalJ = yStart_ + localJ;
  cellCenter[0] = x0_ + h_ * (globalI - 0.5);
  cellCenter[1] = y0_ + h_ * (globalJ - 0.5);
}

// ---------------------
//  get_cell_neighbors
// ---------------------
void
CellularAutomataManager_2D::get_cell_neighbors(int cellLocGhostID, std::vector<int> & neighborVec)
{
  int numNeighbors = 8;
  int i = cellLocGhostID;
  int nx = nxLocGhost_;
  int nbr[8] = { i - nx - 1,
    i - nx,
    i - nx + 1,
    i - 1,
    i + 1,
    i + nx - 1,
    i + nx,
    i + nx + 1 };
  neighborVec.clear();
  neighborVec.assign(&nbr[0], &nbr[0] + numNeighbors);
}

// ----------------
//  get_bulk_cells
// ----------------
void
CellularAutomataManager_2D::get_bulk_cells(std::vector<int> & cellIDs)
{
  // Return vector of all cells, excluding ghosts
  cellIDs.clear();
  cellIDs.resize(numLocalCells_);
  int count = 0;
  for (int j = 1; j < nyLocGhost_ - 1; ++j)
  {
    for (int i = 1; i < nxLocGhost_ - 1; ++i)
    {
      int id = j * nxLocGhost_ + i;
      cellIDs[count++] = id;
    }
  }
}

// -------------------
//  get_surface_cells
// -------------------
void
CellularAutomataManager_2D::get_surface_cells(std::vector<int> & cellIDs)
{
  // Loop over cells and, for surface cells, push back to vector.
  // Corner cells get included twice, and that's the intent.
  cellIDs.clear();
  for (int i = 0; i < numLocGhostCells_; ++i)
  {
    int locGhostI = i % nxLocGhost_;
    int locGhostJ = i / nxLocGhost_;
    int globalI = locGhostI + xStart_;
    int globalJ = locGhostJ + yStart_;
    if (globalI == 1 && globalJ != 0 && globalJ != ny_) cellIDs.push_back(i);//bug here
    else if (globalI == nx_ - 1 && globalJ != 0 && globalJ != ny_) cellIDs.push_back(i);
    if (globalJ == 1 && globalI != 0 && globalI != nx_) cellIDs.push_back(i);//bug here
    else if (globalJ == ny_ - 1 && globalI != 0 && globalI != nx_) cellIDs.push_back(i);
  }
}

// ----------------
//  cell_is_ghost
// ----------------
bool
CellularAutomataManager_2D::cell_is_ghost(int locGhostID)
{
  int locGhostI = locGhostID % nxLocGhost_;
  int locGhostJ = locGhostID / nxLocGhost_;
  bool isGhost = (locGhostI == 0 ||
    locGhostI == nxLocGhost_ - 1 ||
    locGhostJ == 0 ||
    locGhostJ == nyLocGhost_ - 1);
  return isGhost;
}

/*---------------------
  cell_is_global_ghost
  --------------------*/
bool 
CellularAutomataManager_2D::cell_is_global_ghost(int locGhostID)
{
  int locGhostI = locGhostID % nxLocGhost_;
  int locGhostJ = locGhostID / nxLocGhost_;

  int globalGhostI = locGhostI + xStart_;
  int globalGhostJ = locGhostJ + yStart_;

  bool isGlobalGhost = (globalGhostI == 0 ||
                  globalGhostI == nx_ + 1 ||
                        globalGhostJ == 0 ||
                     globalGhostJ == ny_ + 1);
  return isGlobalGhost;
}
// -------------------------------
//  convert_id_local_to_locghost
// -------------------------------
void
CellularAutomataManager_2D::convert_id_local_to_locghost(int localID, int & locGhostID)
{
  int localI = localID % nxLocal_;
  int localJ = localID / nxLocal_;
  locGhostID = (localJ + 1) * nxLocGhost_ + localI + 1;
}

// -------------------------------
//  convert_id_locghost_to_local
// -------------------------------
void
CellularAutomataManager_2D::convert_id_locghost_to_local(int locGhostID, int & localID)
{
  int locGhostI = locGhostID % nxLocGhost_;
  int locGhostJ = locGhostID / nxLocGhost_;
  int localI = locGhostI - 1;
  int localJ = locGhostJ - 1;
  // Make sure this cell is within local boundaries
  assert((localI >= 0) && (localI < nxLocal_));
  assert((localJ >= 0) && (localJ < nyLocal_));
  localID = localJ * nxLocal_ + localI;
}


// -------------------------------
//  convert_id_locghost_to_global
// -------------------------------
void
CellularAutomataManager_2D::convert_id_locghost_to_global(int locGhostID, int & globalID)
{
  int locGhostI = locGhostID % nxLocGhost_;
  int locGhostJ = locGhostID / nxLocGhost_;
  int globalI = locGhostI + xStart_;
  int globalJ = locGhostJ + yStart_;
  globalID = (nx_ + 2)*globalJ + globalI;
}

// -------------------------------
//  convert_id_global_to_locghost
// -------------------------------
void
CellularAutomataManager_2D::convert_id_global_to_locghost(int globalID, int & locGhostID)
{
  int globalI = globalID % (nx_ + 2);
  int globalJ = globalID / (nx_ + 2);
  int locGhostI = globalI - xStart_;
  int locGhostJ = globalJ - yStart_;
  locGhostID = nxLocGhost_ * locGhostJ + locGhostI;
}

/*---------------------------------
convert_id_global_to_locghost
meanwhile check up the position
---------------------------------*/
bool
CellularAutomataManager_2D::convert_id_global_to_locghost_and_check
(int globalID, int & locGhostID)
{
 
  int globalI = globalID % (nx_ + 2);
  int globalJ = globalID / (nx_ + 2);
  int locGhostI = globalI - xStart_;
  if (locGhostI >= nxLocGhost_)
    return false;
  int locGhostJ = globalJ - yStart_;
  if (locGhostJ >= nyLocGhost_)
    return false;
  locGhostID = nxLocGhost_ * locGhostJ + locGhostI;
  return true;
}

// ---------------------------------
//  update_global_orientation_list
// ---------------------------------
void
CellularAutomataManager_2D::update_global_orientation_list(const std::vector<Orientation *> newOrientations,
  const std::vector<int> updateCells)
{
  // Pack orientation angles into an array
  int numLocalOrientations = newOrientations.size();
  std::vector<double> localAngleArray(numLocalOrientations);
  int count = 0;
  for (std::vector<Orientation *>::const_iterator iter = newOrientations.begin();
    iter != newOrientations.end();
    ++iter)
  {
    localAngleArray[count++] = (*iter)->angle_;
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
  for (int i = 1; i < numProcs_; ++i)
  {
    displs[i] = displs[i - 1] + numOrientationsPerProc[i - 1];
  }
  int numGlobalOrientations = displs[numProcs_ - 1] + numOrientationsPerProc[numProcs_ - 1];
  std::vector<double> globalAngleArray(numGlobalOrientations);
  MPI_Allgatherv(&localAngleArray[0],
    numLocalOrientations,
    MPI_DOUBLE,
    &globalAngleArray[0],
    &numOrientationsPerProc[0],
    &displs[0],
    MPI_DOUBLE,
    MPI_COMM_WORLD);

  // Create new Orientations and add to global vector (make sure to
  // just copy the Orientations from the local proc, not create new
  // ones)
  int myDisp = displs[procID_];
  int numCurrentOrientations = orientationVector_.size();
  for (int i = 0; i < numGlobalOrientations; ++i)
  {
    // From current proc
    if ((myDisp <= i) && (i < myDisp + numLocalOrientations))
    {
      newOrientations[i - myDisp]->id_ = numCurrentOrientations + i;
      orientationVector_.push_back(newOrientations[i - myDisp]);
    }
    // From other procs
    else
    {
      orientationVector_.push_back(new Orientation(numCurrentOrientations + i, globalAngleArray[i]));
    }
  }

  // Update cells on local proc with correct global orientation ID
  for (std::vector<int>::const_iterator iter = updateCells.begin();
    iter != updateCells.end();
    ++iter)
  {
    int cellID = *iter;
    int oldOrientationID = cellOrientationID_[cellID];
    int newOrientationID = numCurrentOrientations + myDisp + oldOrientationID;
    cellOrientationID_[cellID] = newOrientationID;
  }

}

// ---------------
//  outputToFile
// ---------------
void
CellularAutomataManager_2D::outputToFile(const std::string & fileName, int * u)
{
  // Communicate all cells (including ghosted boundary cells) to proc 0 for output

  int myProcID, numProcs;
  MPI_Comm_rank(cartComm_, &myProcID);
  MPI_Comm_size(cartComm_, &numProcs);

  int gridSizeX = nxLocGhost_;
  int gridSizeY = nyLocGhost_;
  int xOffset = xStart_;
  int yOffset = yStart_;

  // On process 0, allocate data for entire u array, as well as arrays
  // to hold grid sizes and offsets gathered from individual procs
  int * uAll;
  int *gridSizeXArray, *gridSizeYArray, *xOffsetArray, *yOffsetArray;
  int fullSizeX = nx_ + 2; // Include boundary cells
  int fullSizeY = ny_ + 2; // Include boundary cells
  if (myProcID == 0)
  {
    uAll = new int[fullSizeX*fullSizeY];
    gridSizeXArray = new int[numProcs];
    gridSizeYArray = new int[numProcs];
    xOffsetArray = new int[numProcs];
    yOffsetArray = new int[numProcs];
  }

  // Gather grid sizes and offsets
  MPI_Gather(&gridSizeX, 1, MPI_INTEGER,
    gridSizeXArray, 1, MPI_INTEGER, 0, cartComm_);
  MPI_Gather(&gridSizeY, 1, MPI_INTEGER,
    gridSizeYArray, 1, MPI_INTEGER, 0, cartComm_);
  MPI_Gather(&xOffset, 1, MPI_INTEGER,
    xOffsetArray, 1, MPI_INTEGER, 0, cartComm_);
  MPI_Gather(&yOffset, 1, MPI_INTEGER,
    yOffsetArray, 1, MPI_INTEGER, 0, cartComm_);

  // On each processor, send grid data to process 0. Use non-blocking
  // communication to avoid deadlock.
  MPI_Request request;
  MPI_Isend(u, gridSizeX*gridSizeY, MPI_INTEGER, 0, 0, cartComm_, &request);

  // On process 0, loop over processes and receive the sub-block using
  // a derived data type
  if (myProcID == 0)
  {
    for (int proc = 0; proc < numProcs; ++proc)
    {
      MPI_Datatype subblockType;
      int count = gridSizeYArray[proc];
      int length = gridSizeXArray[proc];
      int stride = fullSizeX;
      MPI_Type_vector(count, length, stride, MPI_INTEGER, &subblockType);
      MPI_Type_commit(&subblockType);
      int * recvPointer = &uAll[yOffsetArray[proc] * fullSizeX + xOffsetArray[proc]];
      MPI_Status status;
      MPI_Recv(recvPointer, 1, subblockType, proc,
        0, cartComm_, &status);
      MPI_Type_free(&subblockType);
    }
  }

  MPI_Wait(&request, MPI_STATUS_IGNORE);
  int * meltTag;
  meltTag = new int[numLocGhostCells_];
  // Output to file from proc 0
  if (myProcID == 0)
  {
    //writeToFile(fileName, uAll, fullSizeX, fullSizeY);
    writeToFile(fileName, uAll, meltTag);
  }

  // Delete arrays on proc 0
  if (myProcID == 0)
  {
    delete uAll;
    delete gridSizeXArray;
    delete gridSizeYArray;
    delete xOffsetArray;
    delete yOffsetArray;
  }

}

// --------------
//  writeToFile
// --------------
void
CellularAutomataManager_2D::writeToFile(const std::string & fileName, 
                                        const int * u, const int * meltTag)
{
  int nx = nx_ + 2;
  int ny = ny_ + 2;
  std::ofstream vtkFile;
  vtkFile.open(fileName.c_str());

  vtkFile << "# vtk DataFile Version 2.0" << std::endl;
  vtkFile << "Test Data" << std::endl;
  vtkFile << "ASCII" << std::endl;
  vtkFile << "DATASET RECTILINEAR_GRID" << std::endl;

  // Dimensions
  vtkFile << "DIMENSIONS "
    << nx_ + 1 << " "
    << ny_ + 1 << " "
    << 1 << std::endl;

  // Coordinates
  vtkFile << "X_COORDINATES " << nx_ + 1 << " double" << std::endl;
  vtkFile << x0_;
  for (int i = 0; i < nx_; ++i)
  {
    vtkFile << " " << x0_ + (i + 1)*h_;
  }
  vtkFile << std::endl;
  vtkFile << "Y_COORDINATES " << ny_ + 1 << " double" << std::endl;
  vtkFile << y0_;
  for (int i = 0; i < ny_; ++i)
  {
    vtkFile << " " << y0_ + (i + 1)*h_;
  }
  vtkFile << std::endl;
  vtkFile << "Z_COORDINATES 1 double" << std::endl;
  vtkFile << "0.0" << std::endl << std::endl;

  // Cell data
  vtkFile << "CELL_DATA " << nx_*ny_ << std::endl;
  vtkFile << "SCALARS grain_id double 1" << std::endl;
  vtkFile << "LOOKUP_TABLE default" << std::endl;
  for (int j = 1; j < ny - 1; ++j)
  {
    for (int i = 1; i < nx - 1; ++i)
    {
      int val = u[j*nx + i];
      vtkFile << val << " ";
    }
  }
  vtkFile << std::endl;

  vtkFile.close();
}

/*--------------------------------------------------------------
   create_nucleation_sites_random_sampling
  --------------------------------------------------------------*/
void 
CellularAutomataManager_2D::
create_nucleation_sites_rejection_random_sampling(double numberDensity,
                                                        double dTcMean,
                                                       double dTcSigma,
                                                bool surfacePopulation)
{
  // Nothing for now
}
