#ifndef PARALLEL_COMM_MANAGER_2D_H
#define PARALLEL_COMM_MANAGER_2D_H

#include <mpi.h>

// Forward declarations
class CellularAutomataManager_2D;

class CommInfo
{
public:
  CommInfo() {};
  CommInfo(double xg, double yg, double length,
    int orientationID, int iCellGlobal)
    : xg_(xg), yg_(yg), length_(length),
    orientationID_(orientationID), iCellGlobal_(iCellGlobal) {};
  ~CommInfo() {};

  double xg_;
  double yg_;
  double length_;
  int orientationID_;
  int iCellGlobal_;
  static int get_num_parameters() { return 5; }
};

class ParallelCommManager_2D
{
public:
  // Constructor/destructor
  ParallelCommManager_2D(CellularAutomataManager_2D * caManager);
  ~ParallelCommManager_2D();

  CellularAutomataManager_2D * caManager_;

  // MPI_Type for communicating comm info
  MPI_Datatype typeCommInfo_;

  // Parallel info
  int numProcs_;
  int procID_;
  MPI_Comm cartComm_;
  int cartRank_;
  int procXP_, procXM_, procYP_, procYM_; // neighboring procs for communication in 4 directions

                      // Vectors to hold comm info
  std::vector<CommInfo> sendXP_, sendXM_, sendYP_, sendYM_;
  std::vector<CommInfo> recvXP_, recvXM_, recvYP_, recvYM_;

  // Set up Cartesian communicator grid
  void setup_cartesian_comm();
  MPI_Comm & get_cart_comm() { return cartComm_; }

  // Pack for communication
  void pack_for_comm(double xg, double yg, double length,
    int orientationID, int iCell);
  void pack_for_comm(double xg, double yg, double length,
    int orientationID, int iCell, double tMCS, double sumT){};

  // Send/receive capture data
  void communicate_capture_data();

  // Handle captures across proc boundaries
  void cross_proc_captures();

  // Exchange ghost data across boundaries
  void exchange_boundary_data(int * u);

protected:

  // Copy vector of CommInfo to double reals, and the reverse
  void comminfo_to_reals(std::vector<CommInfo> & vecCommInfo, std::vector<double> & vecReals);
  void reals_to_comminfo(std::vector<double> & vecReals, std::vector<CommInfo> & vecCommInfo);

  // Handle captures across a particular proc boundary
  void cross_proc_captures(std::vector<CommInfo> & vecCommInfo);
};


#endif
