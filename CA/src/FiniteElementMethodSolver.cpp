#include <yaml-cpp/yaml.h>
#include "CafeParsing.h"
#include "Simulation.h"
#include "FiniteElementMethodSolver.h"
#include "FiniteElementManager.h"
#include "ProblemPhysics.h"
#include "CafeEnv.h"
#include "Timer.h"

/*--------------
   constructor
  --------------*/
FiniteElementMethodSolver::FiniteElementMethodSolver()
{
  // nothing for now
  nDim_ = 3;
}

/*--------------
   destructor
  --------------*/
FiniteElementMethodSolver::~FiniteElementMethodSolver()
{
  if (feManager_)
    delete feManager_;
}

/*-----------
     load 
  -----------*/
void
FiniteElementMethodSolver::
load(const YAML::Node& node)
{
  const YAML::Node * realms = node.FindValue("realms");
  if (realms)
  {
    for (size_t irealm = 0; irealm < realms->size(); irealm++)
    {
      const YAML::Node & realm_node = (*realms)[irealm];
      // check for realm type
      std::string realmType; // = "finite_element";
      get_required(realm_node, "type", realmType);
      get_if_present(realm_node, "dimension", nDim_, nDim_);

      if (realmType == "finite_element")
      {
        if (nDim_ == 2)
        {
          //nothing for now
        }
        else
        {
          feManager_ = new FiniteElementManager();
        }
        feManager_->nDim_ = nDim_;
        //matHeatConduction_ = new Material_HeatConduction(feManager_);
        
        // Load ProblemPhysics first because material properties are needed in feManager_
        //matHeatConduction_->load(realm_node);
        // Load FiniteElementManager
        feManager_->load(realm_node);
      }
    }
  }
  else
    throw std::runtime_error("parser error: realms");

  const YAML::Node * timeIntegration = node.FindValue("time_integrators");
  if (timeIntegration)
  {
    // Make some things pretty in the log file
    CafeEnv::self().caOutputP0() << "\n" << "Time Integration Review For FEM" << "\n";
    CafeEnv::self().caOutputP0() << "=============================" << "\n";

    for (size_t iTimeInte = 0; iTimeInte < timeIntegration->size(); iTimeInte++)
    {
      const YAML::Node &timeInte_node = (*timeIntegration)[iTimeInte];
      const YAML::Node *standardTimeInte_node = timeInte_node.FindValue("standard_time_integrator");
      if (standardTimeInte_node)
      {
        //(*standardTimeInte_node)["name"] >> name_;   // will be used in the future
        (*standardTimeInte_node)["termination_time"] >> finalTime_;
        (*standardTimeInte_node)["time_step_factor"] >> timeStepFactor_;
        (*standardTimeInte_node)["time_step"] >> feManager_->dT_;

        // yaml-cpp Receipt
        CafeEnv::self().caOutputP0() << "Time Integration details gathered from input file and/or defaults:" << "\n";
        CafeEnv::self().caOutputP0() << "Termination Time: " << finalTime_ << "\n";
        CafeEnv::self().caOutputP0() << "Time Step: " << feManager_->dT_ << "\n";
        CafeEnv::self().caOutputP0() << "Time Step Factor: " << timeStepFactor_ << "\n";
      }
    }

  }
  else
    throw std::runtime_error("parser error: tiem_integrations");
}

/*-------------
    initialize
  -------------*/
void
FiniteElementMethodSolver::initialize()
{
  feManager_->initialize();
  timer_->set_FEM_timer();
  feManager_->set_timer_pointer(timer_);
}

/*------------
     execute
  ------------*/
void
FiniteElementMethodSolver::execute()
{
  timer_->start_timing_entire_modeling();
  CafeEnv::self().caOutputP0() << "*********************************************************************";
  CafeEnv::self().caOutputP0() << "\n" << "Simulation shall commence: number of processors = "
    << CafeEnv::self().parallel_size() << "\n";
  CafeEnv::self().caOutputP0() << "*********************************************************************" << std::endl;
  // Integration
  feManager_->setup_integration(timeStepFactor_);
  feManager_->time_integration(finalTime_);
}

/*----------------
    finalize
  ----------------*/
void 
FiniteElementMethodSolver::finalize()
{
  // Output timing information  
  timer_->stop_timing_entire_modeling();

  double maxTimingArray[10];
  double minTimingArray[10];

    timer_->output_timing_info(1,
      maxTimingArray, minTimingArray);
}