#ifndef PARALLEL_COMM_MANAGER_H
#define PARALLEL_COMM_MANAGER_H

#include <mpi.h>
#include "Simulation.h"

// Forward declarations
class CellularAutomataManager;

class ParallelCommManager
{
 public:
  // Constructor/destructor
  //ParallelCommManager(CellularAutomataManager * caManager);
  ParallelCommManager();
  virtual ~ParallelCommManager();

  virtual void load(const YAML::Node& node);
  //Control the direction of the space decomposition 
  bool controledByInputFile_;
  int specifiedDirection_;
  
  // MPI_Type for communicating comm info
  MPI_Datatype typeCommInfo_;
  
  // Parallel info
  int numProcs_;
  int procID_;
  MPI_Comm cartComm_;
  int cartRank_;
  //int procXP_, procXM_, procYP_, procYM_; // neighboring procs for communication in 4 directions

  // Vectors to hold comm info
  //std::vector<CommInfo> sendXP_, sendXM_, sendYP_, sendYM_;
  //std::vector<CommInfo> recvXP_, recvXM_, recvYP_, recvYM_;
  
  // Set up Cartesian communicator grid
  virtual void setup_cartesian_comm()=0;
  MPI_Comm & get_cart_comm() { return cartComm_; }
  virtual void setup_indexed_array() = 0;

  // Pack for communication
  virtual void pack_for_comm(double * coordinates, double length,
                     int orientationID, int iCell, double tMCS, double sumT) = 0;
  virtual void pack_for_comm(double * coordinates, double length,
                     int orientationID, int iCell) = 0;
  // Send/receive capture data
  virtual void communicate_capture_data()=0;

  // Handle captures across proc boundaries
  virtual void cross_proc_captures()=0;

  // Exchange ghost data across boundaries
  virtual void exchange_boundary_data(int * u)=0;
  
  // sending vector to all proc
  virtual void sending_int_vector_to_all_proc(const std::vector<int> &intVector, std::vector<int> &globalVector){throw std::runtime_error("Invalid function: sending_int_vector_to_all_proc in ParallelManager, this function should be rewrite in subclass!");};
  virtual void sending_double_vector_to_all_proc(const std::vector<double> &doubleVector, std::vector<double>& globalVector){throw std::runtime_error("Invalid function: sending_double_vector_to_all_proc in ParallelManager, this function should be rewrite in subclass!");};
 //protected:
  
  // Copy vector of CommInfo to double reals, and the reverse
  //void comminfo_to_reals(std::vector<CommInfo> & vecCommInfo, std::vector<double> & vecReals);
  //void reals_to_comminfo(std::vector<double> & vecReals, std::vector<CommInfo> & vecCommInfo);

  // Handle captures across a particular proc boundary
  //void cross_proc_captures(std::vector<CommInfo> & vecCommInfo);
  
};


#endif
