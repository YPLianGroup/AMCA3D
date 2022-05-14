/*==========================================================
     Define methods used in the FiniteElementManager class
     to load data from fiel including:
     - generate mesh
     - intialize element birth curve
  ==========================================================*/
#include <sstream>
#include <string>
#include <iostream>
#include <mpi.h>
#include "Node.h"
#include "Element.h"
#include "FiniteElementManager.h"
#include "VtuManager.h"
#include "CafeParsing.h"
#include "CafeMacro.h"

/*------------------------------------------
   ~ load_mesh_from_txt_file
  ------------------------------------------*/
void
FiniteElementManager::
load_mesh_and_first_two_steps_results_from_txt_file()
{
  // link the file to the ifstream
  inStreamTracking_.open(inFileName_, std::ios::in);
  if (!inStreamTracking_)
  {
    std::cerr << "fail to open the file" << inFileName_ << std::endl;
    throw std::runtime_error("fail to open the specified file");
  }
  std::string line;

  int dummy[2];
  double time[2];
  int nodeRange[6];

  /* Title */
  int count = 0;
  int skipLines = inNumTitleLines_ + inNumSubtitleLines_; // maybe 11
  while (count < skipLines)
  {
    std::getline(inStreamTracking_, line);
    if (count == inNumTitleLines_+3)
    {
      std::istringstream iss(line);
      
      iss >> dummy[0] >> dummy[1] >> time[0] >> time[1];
      for (int i = 0; i < 6; i++)
        iss >> nodeRange[i];
    }
    count++;
  }

  // store the initial time from the input file
  inTimeOffset_ = time[0];
  inTimeSeries_[0] = 0.0;

  // calculate number of nodes and cells in each direction
  int numNodes[3];
  int numCells[3];

  for (int i = 0; i < 3; i++)
  {    
    numCells[i] = nodeRange[2 * i + 1] - nodeRange[2 * i] ;
    numNodes[i] = numCells[i] +1;
  }

  numNodes_ = numNodes[0] * numNodes[1] * numNodes[2];
  if (numNodes_ == 1)
      throw std::runtime_error("fail to read temperature file, please check the file type and the lines of title!!");
  // Initialize nodeVector_
  while (count < numNodes_ + skipLines)
  {
    std::getline(inStreamTracking_, line);
    std::istringstream iss(line);

    Node node(2);
   
    // get the nodal coordinates
    for (int iCom = 0; iCom < 3; iCom++)
    {
      iss >> node.coordinates_[iCom];
      node.coordinates_[iCom] += positionOffset_[iCom];
      // length unit conversion
      node.coordinates_[iCom] *= lengthScale_;
    }
    // get the nodal temperature of first step
    iss >> node.theta_[0];

    // push the current node into node vector
    nodeVector_.push_back(node);
        
    count++;  // keep on counting the lines number
  }

  // Initialize elementVector_  
  numElements_ = numCells[0] * numCells[1] * numCells[2];
  int nID[8];
  int numNodeXY = numNodes[0] * numNodes[1];
  double coordinates[8][3], iJac[3][3];

  int iEle = 0;
  for (int k = 0; k < numCells[2]; k++)
  {
    for (int j = 0; j < numCells[1]; j++)
    {
      for (int i = 0; i < numCells[0]; i++)
      {
        nID[0] = numNodeXY*k + numNodes[0] * j + i;
        nID[1] = nID[0] + 1;
        nID[2] = nID[1] + numNodes[0];
        nID[3] = nID[2] - 1;
        for (int L = 4; L < 8; L++)
          nID[L] = nID[L - 4] + numNodeXY;

        HexahedralElement* ele = new HexahedralElement;
        double dCell[3];
        for (int M = 0; M < 8; M++)
        {
          ele->nID_[M] = nID[M];
          for (int N = 0; N < 3; N++)
          {
            coordinates[M][N] = nodeVector_[nID[M]].coordinates_[N];
          }
        }
      
        elementVector_.push_back(ele);
        iEle++;                
      } // end i
    } // end j
  } // end k

  inStepRecord_ = 0;
  // load on the second step's result
  load_next_step_nodal_data();
}

/*-------------------------------------------
   ~ load_nodal_data_step_by_step
  -------------------------------------------*/
bool
FiniteElementManager::load_next_step_nodal_data()
{
  std::string line;
  int count = 0;
  double nodalCoord[3];
  double nodalResult;
  while (count < inNumSubtitleLines_ + numNodes_)
  {
    if (! inStreamTracking_.eof()){
      std::getline(inStreamTracking_, line);      
    }
    else {
      return false;
    }
    if (count == 0 && "\r" != line)
        if(""!=line)
            throw std::runtime_error("fail to read next step temperature file, please check title of line!!");
    if (count == 3)
    {
      int dummy[2];
      double time[2];

      std::istringstream iss(line);
      iss >> dummy[0] >> dummy[1] >> time[0] >> time[1];

      // Move current time stamp backward
      inTimeSeries_[1] = inTimeSeries_[0];
      // Update current time stamp 
      inTimeSeries_[0] = (time[0] - inTimeOffset_)*timeScale_;
    }

    // Read in nodal temperature
    if (count >= inNumSubtitleLines_)
    {
      std::istringstream iss(line);
      iss >> nodalCoord[0] >> nodalCoord[1] >> nodalCoord[2] >> nodalResult;
      int nodeID = count - inNumSubtitleLines_;
      Node * node = &nodeVector_[nodeID];
      // Move current nodal result backward
      node->theta_[1] = node->theta_[0];
      // Update current nodal result
      node->theta_[0] = nodalResult;
    }
    count++;
  }

  inStepRecord_++;
  return true;
}


/*-------------------------------------------
   ~ load_nodal_data_step_by_step
  -------------------------------------------*/
double
FiniteElementManager::get_max_time_step()
{

    if (!solveProblemMyself_)
        return inTimeSeries_[0] - time_;
    else
        return get_time_step();
}

bool
FiniteElementManager::pseudo_single_integration_step()
{
  // initial FEM time
  if (time_ == 0)
    time_ = inTimeSeries_[1];
  // Update time
  time_ += dT_;
  // check up time and load next step result if necessary
  while (time_ >= inTimeSeries_[0])
  {
      if (!load_next_step_nodal_data())
      {
        CafeEnv::self().caOutputP0() << "*********************************************************************" << "\n";
        CafeEnv::self().caOutputP0() << "\n" << "Simulation terminated as reading to end of the FLOW3D file" << "\n";
        CafeEnv::self().caOutputP0() << "*********************************************************************" << std::endl;
        return false;
      }
  }
  // inTimeSereis_: 0 for the new one, 1 for the old one
  double timeRatio = (time_ - inTimeSeries_[1]) / (inTimeSeries_[0] - inTimeSeries_[1]);
  if (timeRatio < 0.0)
  {
      throw std::runtime_error("time ratio error: timeRatio is less then 0\ntime 0: " + std::to_string(inTimeSeries_[0]) + "\ntime 1: " + std::to_string(inTimeSeries_[1]));
  }
  
  // update the step count first
  iStep_++;
  bool output = (iStep_ % vtuManager_->outputFreq_ == 0);
  
  if (output)
  {
    FEMtimer_->start_timing_FEM_nodal_time_integration();

    // initialize array used in MPI_Allreduce    
    for (int i = 0; i < numNodes_; i++)
    {
      arrayForAllReduce_[i] = 0.0;
      thetaArray_[i] = 0.0;
    }

    // only loop over local nodes
    for (int iNode = nodeStart_; iNode < nodeEnd_; iNode++) {
      Node* node = &nodeVector_[iNode];
	  // edit nullMatTempToLiquid, also exit in MapVoxelManager.cpp
      bool nullMatTempToLiquid = false;
      if (node->theta_[0] <= nullMatTemperature_) {
        arrayForAllReduce_[iNode] = nullMatTemperature_;
      }
      else if(node->theta_[1] <= nullMatTemperature_){
        arrayForAllReduce_[iNode] = node->theta_[0];
      }
      else {
        arrayForAllReduce_[iNode] = timeRatio*(node->theta_[0] - node->theta_[1]) + node->theta_[1];
      }
    }
    
    FEMtimer_->stop_timing_FEM_nodal_time_integration();
    FEMtimer_->stop_timing_FEM_simulation();

    // equal to broadcast
    MPI_Allreduce(arrayForAllReduce_, thetaArray_, numNodes_,
      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  else 
  {
    FEMtimer_->start_timing_FEM_nodal_time_integration();

    // only loop over nodes with grain around
    int numbGANodes = grainAroundNodeVector_.size();
    for (int iCount = 0; iCount < numbGANodes; iCount++)
    {
      int iNode = grainAroundNodeVector_[iCount];
      Node* node = &nodeVector_[iNode];
      if (node->theta_[0] <= nullMatTemperature_) {
        thetaArray_[iNode] = nullMatTemperature_;
        // avoid the temperature interpolation between a normal cell and nullmat cell
      }
      else if(node->theta_[1] <= nullMatTemperature_){
        thetaArray_[iNode] = node->theta_[0];
      }
      else {
        thetaArray_[iNode] = timeRatio*(node->theta_[0] - node->theta_[1]) + node->theta_[1];
      }
    }

    FEMtimer_->stop_timing_FEM_nodal_time_integration();
    FEMtimer_->stop_timing_FEM_simulation();
  }

  // call vtuManager
  if (output) {
    vtuManager_->execute(iStep_, time_);
  }

  /*-----------------------------------------------------------
     Output current time step log
    -----------------------------------------------------------*/
  CafeEnv::self().caOutputP0() << "************************************************"
                               << "\n";
  CafeEnv::self().caOutputP0() << "Finite Element Method Solver Part: " << "\n";
  CafeEnv::self().caOutputP0() << "time step count = " << iStep_
                               << " current time = " << time_ << "\n"
                               << " time step = " << dT_ << std::endl;

  return true;
}


/*--------------------------------------------
    ~ get_output_time_step_from_txt_file
  --------------------------------------------*/
double
FiniteElementManager::get_output_time_step_from_txt_file()
{
  double timeStep = inTimeSeries_[0] - inTimeSeries_[1];
  return timeStep;
}

