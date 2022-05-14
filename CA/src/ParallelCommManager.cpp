// General includes
#include <mpi.h>

// Local includes
#include "CellularAutomataManager.h"
#include "ParallelCommManager.h"

/*=======================================================================
   Class Definition
    ParallelCommManager - manage parallel part
  =======================================================================*/

/*------------------
   Constructor
  ------------------*/
ParallelCommManager::ParallelCommManager()
  :controledByInputFile_(false)
{
  // Parallel setup
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs_);
  MPI_Comm_rank(MPI_COMM_WORLD, &procID_);
}  

/*-------------
   Destructor
  -------------*/
ParallelCommManager::~ParallelCommManager()
{
  // Nothing for now
}

/*-------------
    load
  -------------*/
void
ParallelCommManager::load(const YAML::Node& node)
{

}

