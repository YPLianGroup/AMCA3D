#ifndef PARALLEL_COMM_MANAGER_3D_H
#define PARALLEL_COMM_MANAGER_3D_H

#include <mpi.h>
// Base class
#include "ParallelCommManager.h"
// Forward delcarations
class CellularAutomataManager_3D;

class CommInfo_3D
{
public:
  CommInfo_3D() {};
  CommInfo_3D(double xg, double yg, double zg, double length,
    int orientationID, int iCellGlobal)
    : xg_(xg), yg_(yg), zg_(zg), length_(length),
    orientationID_(orientationID), iCellGlobal_(iCellGlobal), sumT_(0) ,tMCS_(0) {};
  CommInfo_3D(double xg, double yg, double zg, double length,
    int orientationID, int iCellGlobal,double sumT,double tMCS)
    : xg_(xg), yg_(yg), zg_(zg), length_(length),
    orientationID_(orientationID), iCellGlobal_(iCellGlobal), sumT_(sumT) ,tMCS_(tMCS) {};
  ~CommInfo_3D() {};

  double xg_;
  double yg_;
  double zg_;
  double length_;
  int orientationID_;
  int iCellGlobal_;
  double sumT_;
  double tMCS_;
  //static int get_num_parameters() { return 6; }
  static int get_num_parameters() { return 8;}
};

class ParallelCommManager_3D : public ParallelCommManager
{
public:
  // Construction/destructor
  ParallelCommManager_3D(CellularAutomataManager_3D * caManager);
  ~ParallelCommManager_3D();

  CellularAutomataManager_3D * caManager_;

  // Parallel info: neighboring procs for communication in 6 directions
  int procXP_, procXM_, procYP_, procYM_, procZP_, procZM_;
  int *xLayer_blockLengths;
  int *yLayer_blockLengths;
  int *xLayer_displacements;
  int *yLayer_displacements;
  
  // Vectors to hold comm info
  std::vector<CommInfo_3D> sendXP_, sendXM_, sendYP_, sendYM_, sendZP_, sendZM_;
  std::vector<CommInfo_3D> recvXP_, recvXM_, recvYP_, recvYM_, recvZP_, recvZM_;
  std::vector<int> sendCellIDXP_, sendCellIDXM_, sendCellIDYP_, sendCellIDYM_, sendCellIDZP_, sendCellIDZM_;
  std::vector<CommInfo_3D> recvCellIDXP_, recvCellIDXM_, recvCellIDYP_, recvCellIDYM_, recvCellIDZP_, recvCellIDZM_;
  // Set up Cartesian communicator grid
  virtual void setup_cartesian_comm();
  virtual void setup_indexed_array();
  // Pack for communication
  virtual void pack_for_comm(double * positions, double length,
    int orientationID, int iCell);
  virtual void pack_for_comm(double * positions, double length,
    int orientationID, int iCell, double tMCS, double sumT);
  // Send/receive capture data
  virtual void communicate_capture_data();

  // Handle captures across proc boundaries
  virtual void cross_proc_captures();

  // Exchange ghost data across boundaries
  virtual void exchange_boundary_data(int * u);
  virtual void exchange_boundary_data(double * u);

  // sending vector to all proc
  virtual void sending_int_vector_to_all_proc(const std::vector<int> &intVector, std::vector<int> &globalVector);
  virtual void sending_double_vector_to_all_proc(const std::vector<double> &doubleVector, std::vector<double>& globalVector);

protected:

  // Copy vector of CommInfo to double reals, and the reverse
  void comminfo_to_reals(std::vector<CommInfo_3D> & vecCommInfo, std::vector<double> & vecReals);
  void reals_to_comminfo(std::vector<double> & vecReals, std::vector<CommInfo_3D> & vecCommInfo);

  // Handle captures across a particular proc boundary
  void cross_proc_captures(std::vector<CommInfo_3D> & vecCommInfo);

};

#endif
