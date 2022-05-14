// General includes
#include <mpi.h>
#include <assert.h>

// Local includes
#include "CellularAutomataManager_2D.h"
#include "ParallelCommManager_2D.h"

/*=======================================================================
   Class Definition
     ParallelCommManager_2D - manage parallel part for 2D problem
  =======================================================================*/

/*------------------
   Constructor
  ------------------*/
ParallelCommManager_2D::ParallelCommManager_2D(CellularAutomataManager_2D * caManager)
  : caManager_(caManager)
{
  // Parallel setup
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs_);
  MPI_Comm_rank(MPI_COMM_WORLD, &procID_);
}


/*-------------
  Destructor
  -------------*/
ParallelCommManager_2D::~ParallelCommManager_2D()
{
}

/*---------------------------
  setup_cartesian_comm
  ---------------------------*/
void
ParallelCommManager_2D::setup_cartesian_comm()
{
  int ndims = 2;
  int dims[2] = { 0, 0 };
  MPI_Dims_create(numProcs_, ndims, dims);
  int periods[2] = { 0, 0 };
  int reorder = 1;
  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims,
    periods, reorder, &cartComm_);
  MPI_Comm_rank(cartComm_, &cartRank_);
  MPI_Cart_shift(cartComm_, 0, 1, &procXM_, &procXP_);
  MPI_Cart_shift(cartComm_, 1, 1, &procYM_, &procYP_);

}

/*---------------------------
   pack_for_comm
  ---------------------------*/
void
ParallelCommManager_2D::pack_for_comm(double xg,
  double yg,
  double length,
  int orientationID,
  int iCell)
{
  // Determine which boundary (if any) this cell is on
  int nx = caManager_->nxLocGhost_;
  int ny = caManager_->nyLocGhost_;
  int I = iCell % nx;
  int J = iCell / nx;
  std::vector<CommInfo> * commVector = NULL;
  if (I == 0) commVector = &sendXM_;
  else if (I == nx - 1) commVector = &sendXP_;
  else if (J == 0) commVector = &sendYM_;
  else if (J == ny - 1) commVector = &sendYP_;
  assert(commVector != NULL);

  // Convert cell ID to global
  int iCellGlobal;
  caManager_->convert_id_locghost_to_global(iCell, iCellGlobal);

  // Push the CommInfo onto the correct vector
  CommInfo commInfo(xg, yg, length, orientationID, iCellGlobal);
  commVector->push_back(commInfo);

}

/*---------------------------
   communicate_capture_data
  ---------------------------*/
void
ParallelCommManager_2D::communicate_capture_data()
{
  // Convert vectors of CommInfo to reals
  std::vector<double> sendXPf, sendXMf, sendYPf, sendYMf;
  comminfo_to_reals(sendXP_, sendXPf);
  comminfo_to_reals(sendXM_, sendXMf);
  comminfo_to_reals(sendYP_, sendYPf);
  comminfo_to_reals(sendYM_, sendYMf);

  // Communicate sizes for send/recv
  int nSendXP = sendXPf.size();
  int nSendXM = sendXMf.size();
  int nSendYP = sendYPf.size();
  int nSendYM = sendYMf.size();
  int nRecvXP = 0, nRecvXM = 0, nRecvYP = 0, nRecvYM = 0;
  MPI_Status status;
  MPI_Sendrecv(&nSendXP, 1, MPI_INTEGER, procXP_, 0,
    &nRecvXM, 1, MPI_INTEGER, procXM_, 0, cartComm_, &status);
  MPI_Sendrecv(&nSendXM, 1, MPI_INTEGER, procXM_, 0,
    &nRecvXP, 1, MPI_INTEGER, procXP_, 0, cartComm_, &status);
  MPI_Sendrecv(&nSendYP, 1, MPI_INTEGER, procYP_, 0,
    &nRecvYM, 1, MPI_INTEGER, procYM_, 0, cartComm_, &status);
  MPI_Sendrecv(&nSendYM, 1, MPI_INTEGER, procYM_, 0,
    &nRecvYP, 1, MPI_INTEGER, procYP_, 0, cartComm_, &status);

  // Communicate comminfo buffers
  std::vector<double> recvXPf, recvXMf, recvYPf, recvYMf;
  recvXPf.resize(nRecvXP);
  recvXMf.resize(nRecvXM);
  recvYPf.resize(nRecvYP);
  recvYMf.resize(nRecvYM);
  MPI_Sendrecv(&sendXPf[0], nSendXP, MPI_DOUBLE, procXP_, 0,
    &recvXMf[0], nRecvXM, MPI_DOUBLE, procXM_, 0, cartComm_, &status);
  MPI_Sendrecv(&sendXMf[0], nSendXM, MPI_DOUBLE, procXM_, 0,
    &recvXPf[0], nRecvXP, MPI_DOUBLE, procXP_, 0, cartComm_, &status);
  MPI_Sendrecv(&sendYPf[0], nSendYP, MPI_DOUBLE, procYP_, 0,
    &recvYMf[0], nRecvYM, MPI_DOUBLE, procYM_, 0, cartComm_, &status);
  MPI_Sendrecv(&sendYMf[0], nSendYM, MPI_DOUBLE, procYM_, 0,
    &recvYPf[0], nRecvYP, MPI_DOUBLE, procYP_, 0, cartComm_, &status);

  // Convert back to CommInfo
  reals_to_comminfo(recvXPf, recvXP_);
  reals_to_comminfo(recvXMf, recvXM_);
  reals_to_comminfo(recvYPf, recvYP_);
  reals_to_comminfo(recvYMf, recvYM_);

  // Clear send vectors
  sendXP_.clear();
  sendXM_.clear();
  sendYP_.clear();
  sendYM_.clear();

}

/*----------------------
   cross_proc_captures
  ----------------------*/
void
ParallelCommManager_2D::cross_proc_captures()
{
  // Handle captures across each proc boundary
  communicate_capture_data();
  cross_proc_captures(recvXP_);
  cross_proc_captures(recvXM_);
  cross_proc_captures(recvYP_);
  cross_proc_captures(recvYM_);
  exchange_boundary_data(caManager_->cellOrientationID_);

  // Corner cases require that we do all this one more time
  communicate_capture_data();
  cross_proc_captures(recvXP_);
  cross_proc_captures(recvXM_);
  cross_proc_captures(recvYP_);
  cross_proc_captures(recvYM_);
  exchange_boundary_data(caManager_->cellOrientationID_);

  // At this point communication lists should be clear (in 3D, require one more iteration)
  assert(sendXP_.size() == 0);
  assert(sendXM_.size() == 0);
  assert(sendYP_.size() == 0);
  assert(sendYM_.size() == 0);

}

/*----------------------------------------
   cross_proc_captures
  ----------------------------------------*/
void
ParallelCommManager_2D::cross_proc_captures(std::vector<CommInfo> & vecCommInfo)
{
  for (int i = 0; i < vecCommInfo.size(); ++i)
  {
    CommInfo commInfo = vecCommInfo[i];
    double xg = commInfo.xg_;
    double yg = commInfo.yg_;
    double length = commInfo.length_;
    int orientationID = commInfo.orientationID_;
    int iCell;
    caManager_->convert_id_global_to_locghost(vecCommInfo[i].iCellGlobal_, iCell);
    bool cellIsGhost = caManager_->cell_is_ghost(iCell);
    bool isLiquid = (caManager_->cellGrainID_[iCell] == -1);
    if (isLiquid)
    {
      if (cellIsGhost)
      {
        pack_for_comm(xg, yg, length, orientationID, iCell);
      }
      else
      {
        caManager_->capture_cell(xg, yg, length, orientationID, iCell);
      }
    }
  }

}

/*-------------------------
   exchange_boundary_data
  -------------------------*/
void
ParallelCommManager_2D::exchange_boundary_data(int * u)
{
  int xSize = caManager_->nxLocGhost_;
  int ySize = caManager_->nyLocGhost_;

  // Define row and column data types
  MPI_Datatype rowType, colType;
  MPI_Type_contiguous(xSize, MPI_INTEGER, &rowType);
  MPI_Type_vector(ySize, 1, xSize, MPI_INTEGER, &colType);
  MPI_Type_commit(&rowType);
  MPI_Type_commit(&colType);

  // Get the processes at right, left, up and down
  int procR, procL, procU, procD;
  MPI_Cart_shift(cartComm_, 0, 1, &procL, &procR);
  MPI_Cart_shift(cartComm_, 1, 1, &procD, &procU);

  MPI_Status status;

  // Communicate up (send my top row up, recv my bottom boundary from down)
  MPI_Sendrecv(&u[xSize * (ySize - 2)], 1, rowType, procU, 0,
    &u[0], 1, rowType, procD, 0,
    cartComm_, &status);

  // Communicate down (send my bottom row down, recv my top boundary from up)
  MPI_Sendrecv(&u[xSize], 1, rowType, procD, 0,
    &u[xSize * (ySize - 1)], 1, rowType, procU, 0,
    cartComm_, &status);

  // Communicate right (send my right col right, recv my left boundary from left)
  MPI_Sendrecv(&u[xSize - 2], 1, colType, procR, 0,
    &u[0], 1, colType, procL, 0,
    cartComm_, &status);

  // Communicate left (send my left col left, recv my right boundary from right)
  MPI_Sendrecv(&u[1], 1, colType, procL, 0,
    &u[xSize - 1], 1, colType, procR, 0,
    cartComm_, &status);

  MPI_Type_free(&rowType);
  MPI_Type_free(&colType);

}

/*-------------------
   comminfo_to_reals
  -------------------*/
void
ParallelCommManager_2D::comminfo_to_reals(std::vector<CommInfo> & vecCommInfo,
  std::vector<double> & vecReals)
{
  int numParams = CommInfo::get_num_parameters();
  int numInfo = vecCommInfo.size();
  vecReals.resize(numParams*numInfo);
  int j = 0;
  for (int i = 0; i < numInfo; ++i)
  {
    vecReals[j++] = vecCommInfo[i].xg_;
    vecReals[j++] = vecCommInfo[i].yg_;
    vecReals[j++] = vecCommInfo[i].length_;
    vecReals[j++] = (double)(vecCommInfo[i].orientationID_ + 0.5); // add 0.5 to avoid truncation error
    vecReals[j++] = (double)(vecCommInfo[i].iCellGlobal_ + 0.5);
  }
}


/*-------------------
   reals_to_comminfo
  -------------------*/
void
ParallelCommManager_2D::reals_to_comminfo(std::vector<double> & vecReals,
  std::vector<CommInfo> & vecCommInfo)
{
  int numParams = CommInfo::get_num_parameters();
  int numInfo = vecReals.size() / numParams;
  vecCommInfo.resize(numInfo);
  int j = 0;
  for (int i = 0; i < numInfo; ++i)
  {
    vecCommInfo[i].xg_ = vecReals[j++];
    vecCommInfo[i].yg_ = vecReals[j++];
    vecCommInfo[i].length_ = vecReals[j++];
    vecCommInfo[i].orientationID_ = (int)vecReals[j++];
    vecCommInfo[i].iCellGlobal_ = (int)vecReals[j++];
  }

}
