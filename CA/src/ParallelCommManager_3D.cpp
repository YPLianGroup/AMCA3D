// General includes
#include <mpi.h>
#include <assert.h>
#include <iostream>
#include <sstream>

// Local includes
#include "CellularAutomataManager_3D.h"
#include "ParallelCommManager_3D.h"
#include "ProblemPhysics.h"

/*=======================================================================
   Class Definition
     ParallelCommManager_3D - manage parallel part for 3D problem
  =======================================================================*/

/*------------------
  Constructor
  ------------------*/
ParallelCommManager_3D::ParallelCommManager_3D(CellularAutomataManager_3D * caManager)
  : caManager_(caManager)
{
  // Parallel setup
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs_);
  MPI_Comm_rank(MPI_COMM_WORLD, &procID_);
}

/*-------------
  Destructor
  -------------*/
ParallelCommManager_3D::~ParallelCommManager_3D()
{
  if(xLayer_blockLengths)
     delete[] xLayer_blockLengths;
  if(yLayer_blockLengths)
     delete[] yLayer_blockLengths;
  if(xLayer_displacements)
     delete[] xLayer_displacements;
  if(yLayer_displacements)
     delete[] yLayer_displacements;
}


/*------------------------
   setup_cartesian_comm
  ------------------------*/
void
ParallelCommManager_3D::setup_cartesian_comm()
{
  int ndims = 3;
  int dims[3] = { 0, 0, 0 };
  
  if (controledByInputFile_) {
    for (int i = 0; i < 3; i++) {
      dims[i] = 1;
    }
    dims[specifiedDirection_] = numProcs_;
  }
  else {
    MPI_Dims_create(numProcs_, ndims, dims);
  }
  
  int periods[3] = { 0, 0, 0 };

  int reorder = 1;
  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims,
    periods, reorder, &cartComm_);
  MPI_Comm_rank(cartComm_, &cartRank_);
  MPI_Cart_shift(cartComm_, 0, 1, &procXM_, &procXP_);
  MPI_Cart_shift(cartComm_, 1, 1, &procYM_, &procYP_);
  MPI_Cart_shift(cartComm_, 2, 1, &procZM_, &procZP_);

}
/*-----------------------------
   setup_indexed_array
  -----------------------------*/
void
ParallelCommManager_3D::setup_indexed_array()
{
  int xSize = caManager_->nxLocGhost_;
  int ySize = caManager_->nyLocGhost_;
  int zSize = caManager_->nzLocGhost_;
  int xySize = xSize * ySize;

  xLayer_blockLengths = new int[ySize*zSize];
  yLayer_blockLengths = new int[xSize*zSize];
  xLayer_displacements = new int[ySize*zSize];
  yLayer_displacements = new int[xSize*zSize];

  // Define the block lengths and displacements array for X-layer
  int count = 0;
  for (int k = 0; k < zSize; k++)
  {
    for (int j = 0; j < ySize; j++)
    {
      xLayer_blockLengths[count] = 1;
      xLayer_displacements[count] = j * xSize + k * xySize;
      count++;
    }
  }
  // Define the block lengths and displacements array for Y-layer
  count = 0;
  for (int k = 0; k < zSize; k++)
  {
    for (int i = 0; i < xSize; i++)
    {
      yLayer_blockLengths[count] = 1;
      yLayer_displacements[count] = i + k * xySize;
      count++;
    }
  }
}

/*-----------------
   pack_for_comm
  -----------------*/
void
ParallelCommManager_3D::pack_for_comm(double* coordinates,
  double length,
  int orientationID,
  int iCell, double tMCS, double sumT)
{
  // Determine which boundary (if any) this cell is on
  int nx = caManager_->nxLocGhost_;
  int ny = caManager_->nyLocGhost_;
  int nz = caManager_->nzLocGhost_;
  int nxy = caManager_->nxyLocGhost_;

  int K = iCell / nxy;
  int I = (iCell % nxy) % nx;
  int J = (iCell % nxy) / nx;
  std::vector<CommInfo_3D> * commVector = NULL;

  if (I == 0) commVector = &sendXM_;
  else if (I == nx - 1) commVector = &sendXP_;
  else if (J == 0) commVector = &sendYM_;
  else if (J == ny - 1) commVector = &sendYP_;
  else if (K == 0) commVector = &sendZM_;
  else if (K == nz - 1) commVector = &sendZP_;
  assert(commVector != NULL);

  // Convert cell ID to global
  int iCellGlobal;
  caManager_->convert_id_locghost_to_global(iCell, iCellGlobal);

  // Push the CommInfo onto the correct vector
  double xg = coordinates[0];
  double yg = coordinates[1];
  double zg = coordinates[2];
  CommInfo_3D commInfo(xg, yg, zg, length, orientationID, iCellGlobal, sumT, tMCS);
  commVector->push_back(commInfo);
}

/*-----------------
   pack_for_comm
  -----------------*/
void
ParallelCommManager_3D::pack_for_comm(double* coordinates,
                                      double length,
                                      int orientationID,
                                      int iCell)
{
    // Determine which boundary (if any) this cell is on
    int nx = caManager_->nxLocGhost_;
    int ny = caManager_->nyLocGhost_;
    int nz = caManager_->nzLocGhost_;
    int nxy = caManager_->nxyLocGhost_;

    int K = iCell / nxy;
    int I = (iCell % nxy) % nx;
    int J = (iCell % nxy) / nx;
    std::vector<CommInfo_3D> * commVector = NULL;

    if (I == 0) commVector = &sendXM_;
    else if (I == nx - 1) commVector = &sendXP_;
    else if (J == 0) commVector = &sendYM_;
    else if (J == ny - 1) commVector = &sendYP_;
    else if (K == 0) commVector = &sendZM_;
    else if (K == nz - 1) commVector = &sendZP_;
    assert(commVector != NULL);

    // Convert cell ID to global
    int iCellGlobal;
    caManager_->convert_id_locghost_to_global(iCell, iCellGlobal);

    // Push the CommInfo onto the correct vector
    double xg = coordinates[0];
    double yg = coordinates[1];
    double zg = coordinates[2];
    CommInfo_3D commInfo(xg, yg, zg, length, orientationID, iCellGlobal);
    commVector->push_back(commInfo);
}

/*---------------------------
  communicate_capture_data
  ---------------------------*/
void
ParallelCommManager_3D::communicate_capture_data()
{
  // Convert vectors of CommInfo to reals
  std::vector<double> sendXPf, sendXMf, sendYPf, sendYMf, sendZMf, sendZPf;
  comminfo_to_reals(sendXP_, sendXPf);
  comminfo_to_reals(sendXM_, sendXMf);
  comminfo_to_reals(sendYP_, sendYPf);
  comminfo_to_reals(sendYM_, sendYMf);
  comminfo_to_reals(sendZP_, sendZPf);
  comminfo_to_reals(sendZM_, sendZMf);

  // Communicate sizes for send/recv
  int nSendXP = sendXPf.size();
  int nSendXM = sendXMf.size();
  int nSendYP = sendYPf.size();
  int nSendYM = sendYMf.size();
  int nSendZP = sendZPf.size();
  int nSendZM = sendZMf.size();
  int nRecvXP = 0, nRecvXM = 0, nRecvYP = 0, nRecvYM = 0, nRecvZP = 0, nRecvZM = 0;
  MPI_Status status;
  MPI_Sendrecv(&nSendXP, 1, MPI_INTEGER, procXP_, 0,
    &nRecvXM, 1, MPI_INTEGER, procXM_, 0, cartComm_, &status);
  MPI_Sendrecv(&nSendXM, 1, MPI_INTEGER, procXM_, 0,
    &nRecvXP, 1, MPI_INTEGER, procXP_, 0, cartComm_, &status);
  MPI_Sendrecv(&nSendYP, 1, MPI_INTEGER, procYP_, 0,
    &nRecvYM, 1, MPI_INTEGER, procYM_, 0, cartComm_, &status);
  MPI_Sendrecv(&nSendYM, 1, MPI_INTEGER, procYM_, 0,
    &nRecvYP, 1, MPI_INTEGER, procYP_, 0, cartComm_, &status);
  MPI_Sendrecv(&nSendZP, 1, MPI_INTEGER, procZP_, 0,
    &nRecvZM, 1, MPI_INTEGER, procZM_, 0, cartComm_, &status);
  MPI_Sendrecv(&nSendZM, 1, MPI_INTEGER, procZM_, 0,
    &nRecvZP, 1, MPI_INTEGER, procZP_, 0, cartComm_, &status);

  // Communicate comminfo buffers
  std::vector<double> recvXPf, recvXMf, recvYPf, recvYMf, recvZPf, recvZMf;
  recvXPf.resize(nRecvXP);
  recvXMf.resize(nRecvXM);
  recvYPf.resize(nRecvYP);
  recvYMf.resize(nRecvYM);
  recvZPf.resize(nRecvZP);
  recvZMf.resize(nRecvZM);

  MPI_Sendrecv(&sendXPf[0], nSendXP, MPI_DOUBLE, procXP_, 0,
    &recvXMf[0], nRecvXM, MPI_DOUBLE, procXM_, 0, cartComm_, &status);
  MPI_Sendrecv(&sendXMf[0], nSendXM, MPI_DOUBLE, procXM_, 0,
    &recvXPf[0], nRecvXP, MPI_DOUBLE, procXP_, 0, cartComm_, &status);
  MPI_Sendrecv(&sendYPf[0], nSendYP, MPI_DOUBLE, procYP_, 0,
    &recvYMf[0], nRecvYM, MPI_DOUBLE, procYM_, 0, cartComm_, &status);
  MPI_Sendrecv(&sendYMf[0], nSendYM, MPI_DOUBLE, procYM_, 0,
    &recvYPf[0], nRecvYP, MPI_DOUBLE, procYP_, 0, cartComm_, &status);
  MPI_Sendrecv(&sendZPf[0], nSendZP, MPI_DOUBLE, procZP_, 0,
    &recvZMf[0], nRecvZM, MPI_DOUBLE, procZM_, 0, cartComm_, &status);
  MPI_Sendrecv(&sendZMf[0], nSendZM, MPI_DOUBLE, procZM_, 0,
    &recvZPf[0], nRecvZP, MPI_DOUBLE, procZP_, 0, cartComm_, &status);

  // Convert back to CommInfo
  reals_to_comminfo(recvXPf, recvXP_);
  reals_to_comminfo(recvXMf, recvXM_);
  reals_to_comminfo(recvYPf, recvYP_);
  reals_to_comminfo(recvYMf, recvYM_);
  reals_to_comminfo(recvZPf, recvZP_);
  reals_to_comminfo(recvZMf, recvZM_);

  // Clear send vectors
  sendXP_.clear();
  sendXM_.clear();
  sendYP_.clear();
  sendYM_.clear();
  sendZP_.clear();
  sendZM_.clear();
}

/*----------------------
   cross_proc_captures
  ----------------------*/
void
ParallelCommManager_3D::cross_proc_captures()
{
  // change data for grain coarsening
  // bool mc = caManager_->type_;
  // double *MCS = caManager_->problemPhysics_->get_sumT();

  // Handle captures across each proc boundary
  communicate_capture_data();
  cross_proc_captures(recvXP_);
  cross_proc_captures(recvXM_);
  cross_proc_captures(recvYP_);
  cross_proc_captures(recvYM_);
  cross_proc_captures(recvZP_);
  cross_proc_captures(recvZM_);
  exchange_boundary_data(caManager_->cellOrientationID_);
  //if(mc==1 && NULL!=MCS)
  //  exchange_boundary_data(MCS);

  if (!controledByInputFile_) {
    // Edge cases require that we do all this one more time
    communicate_capture_data();
    cross_proc_captures(recvXP_);
    cross_proc_captures(recvXM_);
    cross_proc_captures(recvYP_);
    cross_proc_captures(recvYM_);
    cross_proc_captures(recvZP_);
    cross_proc_captures(recvZM_);
    exchange_boundary_data(caManager_->cellOrientationID_);
    //if(mc==1 && NULL!=MCS)
    //  exchange_boundary_data(MCS);

    // Corner cases require that we do all this one more time
    communicate_capture_data();
    cross_proc_captures(recvXP_);
    cross_proc_captures(recvXM_);
    cross_proc_captures(recvYP_);
    cross_proc_captures(recvYM_);
    cross_proc_captures(recvZP_);
    cross_proc_captures(recvZM_);
    exchange_boundary_data(caManager_->cellOrientationID_);
    //if(mc==1 && NULL!=MCS)
    //  exchange_boundary_data(MCS);
  }
  
  // At this point communication lists should be clear 
  assert(sendXP_.size() == 0);
  assert(sendXM_.size() == 0);
  assert(sendYP_.size() == 0);
  assert(sendYM_.size() == 0);
  assert(sendZP_.size() == 0);
  assert(sendZM_.size() == 0);

}

/*-----------------------------------
   cross_proc_captures
  -----------------------------------*/
void
ParallelCommManager_3D::cross_proc_captures(std::vector<CommInfo_3D> & vecCommInfo)
{
  for (int i = 0; i < vecCommInfo.size(); ++i)
  {
    CommInfo_3D commInfo = vecCommInfo[i];
    double grainCenter[3];
    grainCenter[0] = commInfo.xg_;
    grainCenter[1] = commInfo.yg_;
    grainCenter[2] = commInfo.zg_;
    double length = commInfo.length_;
    int orientationID = commInfo.orientationID_;
    int iCell;
    caManager_->convert_id_global_to_locghost(vecCommInfo[i].iCellGlobal_, iCell);
    bool cellIsGhost = caManager_->cell_is_ghost(iCell);
    //bool isLiquid = (caManager_->cellGrainID_[iCell] == -1);
    bool isLiquid = (caManager_->cellOrientationID_[iCell] == -1);
    if (isLiquid)
    {
      if (cellIsGhost)
          pack_for_comm(grainCenter, length, orientationID, iCell,0,0);
      else
      {
        caManager_->capture_cell(grainCenter, length, orientationID, iCell, true);
      }
    }
  }
}

/*-------------------------
   exchange_boundary_data
  -------------------------*/
void
ParallelCommManager_3D::exchange_boundary_data(int * u)
{
  int xSize = caManager_->nxLocGhost_;
  int ySize = caManager_->nyLocGhost_;
  int zSize = caManager_->nzLocGhost_;
  int xySize = caManager_->nxyLocGhost_;
  int yzSize = ySize*zSize;
  int xzSize = xSize*zSize;

  // Define x-, y-, z-layer data types
  MPI_Datatype xLayerType, yLayerType, zLayerType;
  MPI_Type_contiguous(xySize, MPI_INTEGER, &zLayerType);
  MPI_Type_indexed(yzSize, xLayer_blockLengths, xLayer_displacements, MPI_INTEGER, &xLayerType);
  MPI_Type_indexed(xzSize, yLayer_blockLengths, yLayer_displacements, MPI_INTEGER, &yLayerType);

  MPI_Type_commit(&xLayerType);
  MPI_Type_commit(&yLayerType);
  MPI_Type_commit(&zLayerType);

  // Get the processes at up and down in each direction
  int procXD, procXU, procYD, procYU, procZD, procZU;
  MPI_Cart_shift(cartComm_, 0, 1, &procXD, &procXU);
  MPI_Cart_shift(cartComm_, 1, 1, &procYD, &procYU);
  MPI_Cart_shift(cartComm_, 2, 1, &procZD, &procZU);

  MPI_Status status;

  // Communicate up (send my top z-layer up, recv my bottom z-layer from down)
  MPI_Sendrecv(&u[xySize * (zSize - 2)], 1, zLayerType, procZU, 0,
    &u[0], 1, zLayerType, procZD, 0,
    cartComm_, &status);

  // Communicate down (send my bottom z-layer down, recv my top z-layer from up)
  MPI_Sendrecv(&u[xySize], 1, zLayerType, procZD, 0,
    &u[xySize * (zSize - 1)], 1, zLayerType, procZU, 0,
    cartComm_, &status);

  // Communicate right (send my top x-layer up (right to right),
  //                    recv my bottom x-layer from down (left to left))
  MPI_Sendrecv(&u[xSize - 2], 1, xLayerType, procXU, 0,
    &u[0], 1, xLayerType, procXD, 0,
    cartComm_, &status);

  // Communicate left (send my bottom x-layer down (left to left),
  //                   recv my top x-layer from up (right to right))
  MPI_Sendrecv(&u[1], 1, xLayerType, procXD, 0,
    &u[xSize - 1], 1, xLayerType, procXU, 0,
    cartComm_, &status);

  // Communicate rear (send my top y-layer up (rear to rear),
  //                    recv my bottom y-layer from down (front to front))
  MPI_Sendrecv(&u[xSize*(ySize-2)], 1, yLayerType, procYU, 0,
    &u[0], 1, yLayerType, procYD, 0,
    cartComm_, &status);

  // Communicate front (send my bottom y-layer down (front to front),
  //                   recv my top z-layer from up (rear to rear))
  MPI_Sendrecv(&u[xSize], 1, yLayerType, procYD, 0,
    &u[xSize*(ySize-1)], 1, yLayerType, procYU, 0,
    cartComm_, &status);

  MPI_Type_free(&xLayerType);
  MPI_Type_free(&yLayerType);
  MPI_Type_free(&zLayerType);
}

/*-------------------------
   exchange_boundary_data
  -------------------------*/
void
ParallelCommManager_3D::exchange_boundary_data(double * u)
{
    int xSize = caManager_->nxLocGhost_;
    int ySize = caManager_->nyLocGhost_;
    int zSize = caManager_->nzLocGhost_;
    int xySize = caManager_->nxyLocGhost_;
    int yzSize = ySize*zSize;
    int xzSize = xSize*zSize;

    // Define x-, y-, z-layer data types
    MPI_Datatype xLayerType, yLayerType, zLayerType;
    MPI_Type_contiguous(xySize, MPI_INTEGER, &zLayerType);
    MPI_Type_indexed(yzSize, xLayer_blockLengths, xLayer_displacements, MPI_INTEGER, &xLayerType);
    MPI_Type_indexed(xzSize, yLayer_blockLengths, yLayer_displacements, MPI_INTEGER, &yLayerType);

    MPI_Type_commit(&xLayerType);
    MPI_Type_commit(&yLayerType);
    MPI_Type_commit(&zLayerType);

    // Get the processes at up and down in each direction
    int procXD, procXU, procYD, procYU, procZD, procZU;
    MPI_Cart_shift(cartComm_, 0, 1, &procXD, &procXU);
    MPI_Cart_shift(cartComm_, 1, 1, &procYD, &procYU);
    MPI_Cart_shift(cartComm_, 2, 1, &procZD, &procZU);

    MPI_Status status;

    // Communicate up (send my top z-layer up, recv my bottom z-layer from down)
    MPI_Sendrecv(&u[xySize * (zSize - 2)], 1, zLayerType, procZU, 0,
                 &u[0], 1, zLayerType, procZD, 0,
                 cartComm_, &status);

    // Communicate down (send my bottom z-layer down, recv my top z-layer from up)
    MPI_Sendrecv(&u[xySize], 1, zLayerType, procZD, 0,
                 &u[xySize * (zSize - 1)], 1, zLayerType, procZU, 0,
                 cartComm_, &status);

    // Communicate right (send my top x-layer up (right to right),
    //                    recv my bottom x-layer from down (left to left))
    MPI_Sendrecv(&u[xSize - 2], 1, xLayerType, procXU, 0,
                 &u[0], 1, xLayerType, procXD, 0,
                 cartComm_, &status);

    // Communicate left (send my bottom x-layer down (left to left),
    //                   recv my top x-layer from up (right to right))
    MPI_Sendrecv(&u[1], 1, xLayerType, procXD, 0,
                 &u[xSize - 1], 1, xLayerType, procXU, 0,
                 cartComm_, &status);

    // Communicate rear (send my top y-layer up (rear to rear),
    //                    recv my bottom y-layer from down (front to front))
    MPI_Sendrecv(&u[xSize*(ySize-2)], 1, yLayerType, procYU, 0,
                 &u[0], 1, yLayerType, procYD, 0,
                 cartComm_, &status);

    // Communicate front (send my bottom y-layer down (front to front),
    //                   recv my top z-layer from up (rear to rear))
    MPI_Sendrecv(&u[xSize], 1, yLayerType, procYD, 0,
                 &u[xSize*(ySize-1)], 1, yLayerType, procYU, 0,
                 cartComm_, &status);

    MPI_Type_free(&xLayerType);
    MPI_Type_free(&yLayerType);
    MPI_Type_free(&zLayerType);
}

/*-------------------
   comminfo_to_reals
  -------------------*/
void
ParallelCommManager_3D::comminfo_to_reals(std::vector<CommInfo_3D> & vecCommInfo,
  std::vector<double> & vecReals)
{
  int numParams = CommInfo_3D::get_num_parameters();
  int numInfo = vecCommInfo.size();
  vecReals.resize(numParams*numInfo);
  int j = 0;
  for (int i = 0; i < numInfo; ++i)
  {
    vecReals[j++] = vecCommInfo[i].xg_;
    vecReals[j++] = vecCommInfo[i].yg_;
    vecReals[j++] = vecCommInfo[i].zg_;
    vecReals[j++] = vecCommInfo[i].length_;
    vecReals[j++] = (double)(vecCommInfo[i].orientationID_ + 0.5); // add 0.5 to avoid truncation error
    vecReals[j++] = (double)(vecCommInfo[i].iCellGlobal_ + 0.5);
    vecReals[j++] = vecCommInfo[i].sumT_;
    vecReals[j++] = vecCommInfo[i].tMCS_;
  }
}

/* -------------------
   reals_to_comminfo
  -------------------*/
void
ParallelCommManager_3D::reals_to_comminfo(std::vector<double> & vecReals,
  std::vector<CommInfo_3D> & vecCommInfo)
{
  int numParams = CommInfo_3D::get_num_parameters();
  int numInfo = vecReals.size() / numParams;
  vecCommInfo.resize(numInfo);
  int j = 0;
  for (int i = 0; i < numInfo; ++i)
  {
    vecCommInfo[i].xg_ = vecReals[j++];
    vecCommInfo[i].yg_ = vecReals[j++];
    vecCommInfo[i].zg_ = vecReals[j++];
    vecCommInfo[i].length_ = vecReals[j++];
    vecCommInfo[i].orientationID_ = (int)vecReals[j++];
    vecCommInfo[i].iCellGlobal_ = (int)vecReals[j++];
    vecCommInfo[i].sumT_ = vecReals[j++];
    vecCommInfo[i].tMCS_ = vecReals[j++];
  }

}

/*----------------------------------------
     MPI communication
     sending vector to all proc 
  ----------------------------------------*/
void
ParallelCommManager_3D::
sending_int_vector_to_all_proc(const std::vector<int> &intVector, std::vector<int>& globalVector)
// para : intVector : local vector, this vector should be communicated and stortred in globalVector
// para : globalVector : golabe vector, this vector stortred local vector
{
  // Pack vector into an array
  int lenVector = intVector.size();

  // Communicate with other procs through Allgatherv
  // FIXME: maybe MPI_Gather is better
  std::vector<int> lenVectorPerProc(numProcs_);
  MPI_Allgather(&lenVector,
                1,
                MPI_INTEGER,
                &lenVectorPerProc[0],
                1,
                MPI_INTEGER,
                MPI_COMM_WORLD);

  std::vector<int> displs(numProcs_);
  displs[0] = 0;
  for (int i = 1; i < numProcs_; ++i)
  {
    displs[i] = displs[i-1] + lenVectorPerProc[i-1];
  }
  globalVector.resize(displs[numProcs_-1] + lenVectorPerProc[numProcs_-1]);
  if (true)
  {
    // MPI_Allgatherv means recived number of data is different
    // FIXME: maybe MPI_Gatherv is better
	  MPI_Allgatherv(&intVector[0],
		  lenVector,
		  MPI_INTEGER,
		  &globalVector[0],
		  &lenVectorPerProc[0],
		  &displs[0],
		  MPI_INTEGER,
		  MPI_COMM_WORLD);
  }
}

void
ParallelCommManager_3D::
sending_double_vector_to_all_proc(const std::vector<double> &doubleVector, std::vector<double>& globalVector)
// para : doubleVector : local vector, this vector should be communicated and stortred in globalVector
// para : globe Vector : golabe vector, this vector stortred local vector
{
  // Pack vector into an array
  int lenVector = doubleVector.size();

  // Communicate with other procs through Allgatherv
  // FIXME: maybe MPI_Gather is better
  std::vector<int> lenVectorPerProc(numProcs_);
  MPI_Allgather(&lenVector,
                1,
                MPI_INTEGER,
                &lenVectorPerProc[0],
                1,
                MPI_INTEGER,
                MPI_COMM_WORLD);

  std::vector<int> displs(numProcs_);
  displs[0] = 0;
  for (int i = 1; i < numProcs_; ++i)
  {
    displs[i] = displs[i-1] + lenVectorPerProc[i-1];
  }
  //std::vector<int> globalVector(displs[numProcs_-1]);
  globalVector.resize(displs[numProcs_-1]);
  if (true)
  {
    // MPI_Allgatherv means recived number of data is different
    // FIXME: maybe MPI_Gatherv is better
	  MPI_Allgatherv(&doubleVector[0],
		  lenVector,
		  MPI_DOUBLE,
		  &globalVector[0],
		  &lenVectorPerProc[0],
		  &displs[0],
		  MPI_DOUBLE,
		  MPI_COMM_WORLD);
  }
}
