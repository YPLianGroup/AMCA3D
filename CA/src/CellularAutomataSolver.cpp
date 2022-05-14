#include "CafeParsing.h"
#include "Simulation.h"
#include "CellularAutomataSolver.h"
#include "CellularAutomataManager.h"
#include "CellularAutomataManager_2D.h"
#include "CellularAutomataManager_3D.h"
#include "CellularAutomataManager_3DRemelting.h"
#include "ProblemPhysics.h"
#include "CafeEnv.h"
#include "Timer.h"

#define M_PI 3.14159265358979323846 

/*--------------------------------------------------------------------------
    constructor
  --------------------------------------------------------------------------*/
CellularAutomataSolver::CellularAutomataSolver():
  maxSteps_(1000),
  caManager_(NULL),
  problemPhysics_(NULL)
{
  nDim_ = 3; 
}

/*--------------------------------------------------------------------------
    destructor
  --------------------------------------------------------------------------*/
CellularAutomataSolver::~CellularAutomataSolver()
{
  if (caManager_)
    delete caManager_;
  if (problemPhysics_)
    delete problemPhysics_;
}

/*---------------------------------------------------------------------------
   load
  ---------------------------------------------------------------------------*/
void
CellularAutomataSolver::load(const YAML::Node &node)
{  
  const YAML::Node * realms = node.FindValue("realms");
  if (realms)
  {
    for (size_t irealm = 0; irealm < realms->size(); irealm++)
    {
      const YAML::Node & realm_node = (*realms)[irealm];
      // check for realm type
      std::string realmType; //= "cellular_automata";
      get_required(realm_node, "type", realmType);
      get_if_present(realm_node, "dimension", nDim_, nDim_);
      
      if (realmType == "cellular_automata")
      {
        if (nDim_ == 2)
        {
          caManager_ = new CellularAutomataManager_2D();
        }
        else
        {
          caManager_ = new CellularAutomataManager_3D();
        }
        caManager_->nDim_ = nDim_;
        problemPhysics_ = new ProblemPhysics_RappazGandin1993(caManager_);
        
        // setup pointers about output and parallel managers
        caManager_->setup_several_data_member_pointer();
        
        // Load CellularAutomataManager
        caManager_->load(realm_node);
       
        // Load ProblemPhysics
        problemPhysics_->load(realm_node);
      }
      else if(realmType == "cellular_automata_remelting"){
          caManager_ = new CellularAutomataManager_3DRemelting();

          caManager_->nDim_ = nDim_;
          problemPhysics_ = new ProblemPhysics_RappazGandin1993(caManager_);

          // setup pointers about output and parallel managers
          caManager_->setup_several_data_member_pointer();

          // Load CellularAutomataManager
          caManager_->load(realm_node);

          // Load ProblemPhysics
          problemPhysics_->load(realm_node);
      }
      else
        throw std::runtime_error("parser error: realm-type");
      
    }
  }
  else
    throw std::runtime_error("parser error: realms");

  // FIEME: will move this to some class, say time integration class
  const YAML::Node * timeIntegration = node.FindValue("time_integrators");
  if (timeIntegration)
  {
    // Make some things pretty in the log file
    CafeEnv::self().caOutputP0() << "\n" << "Time Integration Review" << "\n";
    CafeEnv::self().caOutputP0() << "=============================" << "\n";

    for (size_t iTimeInte = 0; iTimeInte < timeIntegration->size(); iTimeInte++)
    {
      const YAML::Node &timeInte_node = (*timeIntegration)[iTimeInte];
      const YAML::Node *standardTimeInte_node = timeInte_node.FindValue("standard_time_integrator");
      if (standardTimeInte_node)
      {
        //(*standardTimeInte_node)["name"] >> name_;   // will be used in the future
        (*standardTimeInte_node)["termination_time"] >> finalTime_;
        (*standardTimeInte_node)["time_step"] >> maxDt_;
        get_if_present(*standardTimeInte_node, "termination_step_count", maxSteps_, maxSteps_);
        get_if_present(*standardTimeInte_node, "time_step_factor_ca", caManager_->timeStepFactor_, caManager_->timeStepFactor_);
        // yaml-cpp Receipt
        CafeEnv::self().caOutputP0() << "Time Integration details gathered from input file and/or defaults:" << "\n";
        CafeEnv::self().caOutputP0() << "Termination Time: " << finalTime_ << "\n";
        CafeEnv::self().caOutputP0() << "Termination Step Count: " << maxSteps_ << "\n";
        CafeEnv::self().caOutputP0() << "Initial Time Step: " << maxDt_ << "\n";
        CafeEnv::self().caOutputP0() << "Time Step Factor for CA: " << caManager_->timeStepFactor_ << "\n";
      }
    }

  }
  else
    throw std::runtime_error("parser error: tiem_integrations");

}

/*---------------------------------------------------------------------------
   initialize 
  ---------------------------------------------------------------------------*/
void
CellularAutomataSolver::initialize()
{
  timer_->set_CA_timer();
  // Set up problem and output files
  caManager_->initialize();
  caManager_->set_timer_pointer(timer_);
}

/*---------------------------------------------------------------------------
   execute
  ---------------------------------------------------------------------------*/
void
CellularAutomataSolver::execute()
{
  timer_->start_timing_entire_modeling();

  CafeEnv::self().caOutputP0() << "*********************************************************************";
  CafeEnv::self().caOutputP0() << "\n" << "Simulation shall commence: number of processors = " 
                                       << CafeEnv::self().parallel_size()<< "\n";
  CafeEnv::self().caOutputP0() << "*********************************************************************" << "\n";;
  // Integration
  caManager_->setup_integration(maxDt_);
  caManager_->time_integration(finalTime_, maxSteps_);
}

/*---------------------------------------------------------------------------
   finalize - will be moved to main()
   FIXME: will move this member function to other place
  ---------------------------------------------------------------------------*/
void
CellularAutomataSolver::finalize()
{
  // Stop timing
  timer_->stop_timing_entire_modeling();

  // Output final number of grains
  int numberOrientation = 0;
  caManager_->get_amount_of_final_orientations(numberOrientation);
  CafeEnv::self().caOutputP0() << "=====================================" << "\n";
  CafeEnv::self().caOutputP0() << " Final Number of Grains:" << "\n";
  CafeEnv::self().caOutputP0() << "             " << numberOrientation << "\n";
  CafeEnv::self().caOutputP0() << "=====================================" << std::endl;
  
  // Output time log
  double maxTimingArray[10];
  double minTimingArray[10];
  for (int i = 0; i < 10; i++) {
    maxTimingArray[i] = 0.0;
    minTimingArray[i] = 0.0;
  }
  MPI_Allreduce(timer_->timingArray_, maxTimingArray, 10,
    MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(timer_->timingArray_, minTimingArray, 10,
    MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  if (caManager_->procID_ == 0)
  {
    timer_->output_timing_info(caManager_->numProcs_,
                                maxTimingArray,minTimingArray);
  }

  // Output final results
  caManager_->output_last_step_result();
  std::string vtkFileName = "finalGrains.vtk";
  caManager_->output_orientations(vtkFileName);

}