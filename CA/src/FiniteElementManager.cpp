/*=======================================================================
                   FiniteElementManager
           Definition of Finite Element Manager class.
  =======================================================================*/
#include <mpi.h>

#include "Node.h"
#include "Element.h"
#include "CafeParsing.h"
#include "FiniteElementManager.h"
#include "VtuManager.h"
#include "CafeMacro.h"
#define __max(a,b)  (((a) > (b)) ? (a) : (b))
#define __min(a,b)  (((a) < (b)) ? (a) : (b))
/*---------------lo
   constructor
  ---------------*/
FiniteElementManager::FiniteElementManager()
  :timeStepFactor_(0.9),
  initialTheta_(0.0),
  time_(0.0),
  iStep_(0),
  nodalForceArray_(NULL),
  arrayForAllReduce_(NULL),
  thetaArray_(NULL),
  solveProblemMyself_(false),
  lengthScale_(1.0),
  timeScale_(1.0),
  skippingHead_(true),
  numNodes_(0),
  numElements_(0),
  microTimeStepStrategy_(false),
  nullMatTemperature_(-1.0e8)
{
    // Parallel setup
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs_);
    MPI_Comm_rank(MPI_COMM_WORLD, &procID_);
    vtuManager_ = new VtuManager(this);
}

/*---------------
    destructor
  ---------------*/
FiniteElementManager::~FiniteElementManager()
{
  if (nodalForceArray_)
    delete[] nodalForceArray_;

  if (arrayForAllReduce_)
    delete[] arrayForAllReduce_;

  if (thetaArray_) {
    delete[] thetaArray_;
  }

  if (vtuManager_)
    delete vtuManager_;

  if (inStreamTracking_)
    inStreamTracking_.close();
}

/*------------
     load
  ------------*/
void
FiniteElementManager::load(const YAML::Node& node)
{
  // Second: get solution options
  const YAML::Node *solutionOptions = node.FindValue("solution_options");
  if (solutionOptions)
  {
    // Make some things pretty in the log file
    CafeEnv::self().caOutputP0() << "\n" << "FEM: Solution Options Review" << std::endl;
    CafeEnv::self().caOutputP0() << "=============================" << std::endl;
    std::string name_;
    get_required(*solutionOptions, "name", name_);

    const YAML::Node * options_node = expect_sequence(*solutionOptions, "options", true);
    if (options_node)
    {
      for (size_t ioption = 0; ioption < options_node->size(); ioption++)
      {
        const YAML::Node & option_node = (*options_node)[ioption];
        if (expect_map(option_node, "load_data_from_file", true))
        {
          bool loadWholeTimeData = false;
          const YAML::Node& loadData = *option_node.FindValue("load_data_from_file");
          get_if_present(loadData, "for_whole_time", loadWholeTimeData, loadWholeTimeData);
          get_if_present(loadData, "lines_for_title", inNumTitleLines_, inNumTitleLines_);
          get_if_present(loadData, "lines_for_subtitle", inNumSubtitleLines_, inNumSubtitleLines_);
          get_if_present(loadData, "length_scale", lengthScale_, lengthScale_);
          get_if_present(loadData, "time_scale", timeScale_, timeScale_);
          get_if_present(loadData, "skipping_head", skippingHead_, skippingHead_);
          std::string keyWord[3] = {"x_offset", "y_offset", "z_offset"};
          for (int iOffset = 0; iOffset < 3; iOffset++) {
            positionOffset_[iOffset] = 0.0;                        
            get_if_present(loadData, keyWord[iOffset], positionOffset_[iOffset], positionOffset_[iOffset]);
          }
          if (loadWholeTimeData)
          {
            solveProblemMyself_ = false;
            CafeEnv::self().caOutputP0() << "FEM solver will not solve the problem" << std::endl;
          }
        }
      } // end for ioption
    } // end for option_nodes
  }

  // Third: get mesh
  const YAML::Node * mesh = node.FindValue("mesh");
  // get parameterized mesh
  const YAML::Node * domain_node = node.FindValue("domain");
  if (mesh)
  {
    CafeEnv::self().caOutputP0() << "\n" << "FEM Domain Size Review" << "\n";
    CafeEnv::self().caOutputP0() << "=============================" << "\n";

    get_required(node, "mesh", inFileName_);
    int matID = 0;
    get_if_present(node, "material_id", matID, matID);
    CafeEnv::self().caOutputP0() << "Load Mesh From File: " << inFileName_ << std::endl;
    
    // check the file name extension    
    std::size_t dotPos = inFileName_.find_last_of(".");
    std::string extension = inFileName_.substr(dotPos+1);

    load_mesh_and_first_two_steps_results_from_txt_file();
    CafeEnv::self().caOutputP0() << "Total node number: " << numNodes_ << "\n";
    CafeEnv::self().caOutputP0() << "Total element number:" << numElements_ << "\n";
    CafeEnv::self().caOutputP0() << "Original Point: " << "\n";
    CafeEnv::self().caOutputP0() << "    x0=" << nodeVector_[0].coordinates_[0]
        << "    y0=" << nodeVector_[0].coordinates_[1]
        << "    z0=" << nodeVector_[0].coordinates_[2] << "\n";
    CafeEnv::self().caOutputP0() << "End Point: " << "\n";
    CafeEnv::self().caOutputP0() << "    x0=" << nodeVector_[nodeVector_.size()-1].coordinates_[0]
        << "    y0=" << nodeVector_[nodeVector_.size()-1].coordinates_[1]
        << "    z0=" << nodeVector_[nodeVector_.size()-1].coordinates_[2] << "\n";
  }
  // get output info
  vtuManager_->load(node);
}

/*-------------
   initialize
  -------------*/
void 
FiniteElementManager::initialize()
{
  // allocate and initialize the arrays
  if (solveProblemMyself_){
    nodalForceArray_ = new double[numNodes_];
    for (int i = 0; i < numNodes_; i++) {
      nodalForceArray_[i] = 0.0;
    }
  }
  arrayForAllReduce_ = new double[numNodes_];
  thetaArray_ = new double[numNodes_];
  for (int i = 0; i < numNodes_; i++)
  {    
    arrayForAllReduce_[i] = 0.0;
    
    if (!solveProblemMyself_) {
      thetaArray_[i] = nodeVector_[i].theta_[1];
    }
    else {
      thetaArray_[i] = initialTheta_;
      if (microTimeStepStrategy_) {
        nodeVector_[i].theta_[0] = initialTheta_;
        nodeVector_[i].theta_[1] = initialTheta_;
      }
    }
      nodeVector_[i].birth_ = true;
  }
 

    for (int iEle = 0; iEle < numElements_; iEle++){
      elementVector_[iEle]->birth_ = true;
    }
  
  // for parallel part
  node_element_list_decomposition();

  // for output
  vtuManager_->initialize();

  // temperature for null material

}

/*----------------------------------
   node_element_list_decomposition 
  ----------------------------------*/
void
FiniteElementManager::node_element_list_decomposition()
{
    // Compute number of local nodes, and offset
    numNodesLocal_ = numNodes_ / numProcs_;
    nodeStart_ = procID_ * numNodesLocal_;
    // Account for non-uniform distribution on last proc: offset + nlocal = n
    if (procID_ == numProcs_ - 1)
    {
        numNodesLocal_ = numNodes_ - nodeStart_;
    }
    nodeEnd_ = nodeStart_ + numNodesLocal_;

    // Compute number of local elements, and offset
    numElementsLocal_ = numElements_ / numProcs_;
    elementStart_ = procID_ * numElementsLocal_;
    // Account for non-uniform distribution on last proc: offset + nlocal = n
    if (procID_ == numProcs_ - 1)
    {
        numElementsLocal_ = numElements_ - elementStart_;
    }
    elementEnd_ = elementStart_ + numElementsLocal_;
}

/*-----------------------
   setup_integration
  -----------------------*/
void 
FiniteElementManager::setup_integration(double timeStepFactor)
{
  if (solveProblemMyself_)
  {
    timeStepFactor_ = timeStepFactor;
  }
  else
  {
    dT_ = get_output_time_step_from_txt_file();
  }
}

/*--------------------------
    set_coupling_time_step
  --------------------------*/
void
FiniteElementManager::setup_time_step_for_coupling(double dT)
{
  dT_ = dT;
}

/*----------------------
     time_integration
  ----------------------*/
void
FiniteElementManager::time_integration(double finalTime)
{
  while (time_ < finalTime) {
    if (iStep_ == 0){
      vtuManager_->execute(iStep_, time_);
    }

    FEMtimer_->start_timing_FEM_simulation();
    dT_ = inTimeSeries_[0] - inTimeSeries_[1];
      if (!pseudo_single_integration_step()) {
        CafeEnv::self().caOutputP0() << "*********************************************************************"<< "\n";
        CafeEnv::self().caOutputP0() << "\n" << "Simulation is terminated because reading to end of the file" << "\n";
        CafeEnv::self().caOutputP0() << "*********************************************************************" << std::endl;
        break;
      }
  } // end while
}

/*--------------------------------
   set_grain_around_node_vector()
  --------------------------------*/
void
FiniteElementManager::set_grain_around_node_vector()
{
  for (int iNode = 0; iNode < numNodes_; iNode++) {
    if (nodeVector_[iNode].grainAround_) {
      grainAroundNodeVector_.push_back(iNode);
    }
  }
}