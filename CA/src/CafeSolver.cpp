#include <sstream>
#include <string>
#include <cmath>
#include <iostream>
#include "CafeParsing.h"
#include "CafeSolver.h"
#include "CellularAutomataManager.h"
#include "CellularAutomataManager_2D.h"
#include "CellularAutomataManager_3D.h"
#include "CellularAutomataManager_3DRemelting.h"
#include "FiniteElementManager.h"
#include "MapVoxelManager.h"
#include "CafeEnv.h"
#include "Timer.h"
#include <mpi.h>
#define __max(a,b)  (((a) > (b)) ? (a) : (b))
#define __min(a,b)  (((a) < (b)) ? (a) : (b))
/*---------------*/
/* ~ constructor  
/*---------------*/
CafeSolver::CafeSolver()
{
//nothing for now
}

/*  ~  destructor   */
CafeSolver::~CafeSolver()
{
  if (feManager_)
    delete feManager_;
  if (caManager_)
    delete caManager_;
  if (mapVoxelManager_)
    delete mapVoxelManager_;
}

/*---------------*/
/*  ~ load  
/*---------------*/
void
CafeSolver::load(const YAML::Node& node)
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

      if (realmType == "cellular_automata" || realmType == "cellular_automata_remelting")
      {
        if (nDim_ == 2)
        {
          caManager_ = new CellularAutomataManager_2D();
        }
        else
        {
          if (realmType == "cellular_automata")
            caManager_ = new CellularAutomataManager_3D();
          else if (realmType == "cellular_automata_remelting")
            caManager_ = new CellularAutomataManager_3DRemelting();
        }
        caManager_->nDim_ = nDim_;
        problemPhysics = new ProblemPhysics_RappazGandin1993(caManager_);

        // setup pointers about output and parallel managers
        caManager_->setup_several_data_member_pointer();

        // Load CellularAutomataManager
        caManager_->load(realm_node);

        // Load ProblemPhysics
        problemPhysics->load(realm_node);
      }
      
      else if (realmType == "finite_element")
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

        // Load FiniteElementManager
        feManager_->load(realm_node);
      }
      else
        throw std::runtime_error("parser error: realm-type");
    }
  }
  else
    throw std::runtime_error("parser error: realms");

  const YAML::Node * timeIntegration = node.FindValue("time_integrators");
  if (timeIntegration)
  {
    // Make some things pretty in the log file
    CafeEnv::self().caOutputP0() << "\n" << "Time Integration Review" << "\n";
    CafeEnv::self().caOutputP0() << "=============================" << std::endl;

    for (size_t iTimeInte = 0; iTimeInte < timeIntegration->size(); iTimeInte++)
    {
      const YAML::Node &timeInte_node = (*timeIntegration)[iTimeInte];
      const YAML::Node *standardTimeInte_node = timeInte_node.FindValue("standard_time_integrator");
      if (standardTimeInte_node)
      {
        //(*standardTimeInte_node)["name"] >> name_;   // will be used in the future
        (*standardTimeInte_node)["termination_time"] >> finalTime_;
        (*standardTimeInte_node)["time_step"] >> maxDt_;
        (*standardTimeInte_node)["time_step_factor"] >> timeStepFactorFE_;
        get_if_present(*standardTimeInte_node, "termination_step_count", maxSteps_, maxSteps_);
        get_if_present(*standardTimeInte_node, "time_step_factor_ca", caManager_->timeStepFactor_, caManager_->timeStepFactor_);
        std::string timeWindowFileName = "none";
        // yaml-cpp Receipt
        CafeEnv::self().caOutputP0() << "Time Integration details gathered from input file and/or defaults:" << "\n";
        CafeEnv::self().caOutputP0() << "Termination Time: " << finalTime_ << "\n";
        CafeEnv::self().caOutputP0() << "Termination Step Count (from CA): " << maxSteps_ << "\n";
        CafeEnv::self().caOutputP0() << "Initial Time Step (from CA): " << maxDt_ << std::endl;
        CafeEnv::self().caOutputP0() << "Time Step Factor for CA: " << caManager_->timeStepFactor_ << "\n";
      }
    }

  }
  else
    throw std::runtime_error("parser error: tiem_integrations");

  mapVoxelManager_ = new MapVoxelManager(feManager_, caManager_);
  mapVoxelManager_->load(node);
}

/*---------------*/
/*  ~ initialize 
/*---------------*/
void
CafeSolver::initialize()
{
  feManager_->initialize();
  caManager_->initialize();
  mapVoxelManager_->initialize();  
  caManager_->setup_map_voxel_manager_and_accessory(mapVoxelManager_);

  timer_->set_FEM_timer();
  timer_->set_CA_timer();
  feManager_->set_timer_pointer(timer_);
  caManager_->set_timer_pointer(timer_);
}

/*---------------*/
/*   ~ execute   
/*---------------*/
void
CafeSolver::execute()
{
  mapVoxelManager_->execute_with_focus_on_CA_region();
  feManager_->set_grain_around_node_vector();

  timer_->start_timing_entire_modeling();

  CafeEnv::self().caOutputP0() << "*********************************************************************";
  CafeEnv::self().caOutputP0() << "\n" << "CAFE Simulation shall commence: number of processors = "
                               << CafeEnv::self().parallel_size() << "\n";
  CafeEnv::self().caOutputP0() << "*********************************************************************" << std::endl;
  // Integration
  time_integration_weak_coupling();
}

/*--------------------------------------------------------------*/
/*  whether the temperature is higher than solidius or not, solid cell growth*/
/*--------------------------------------------------------------*/
void
CafeSolver::time_integration_weak_coupling()
{
  // Get original time step of FEM and use it as the maximum time step of CA
  feManager_->setup_integration(timeStepFactorFE_);  
  double dT_FE = feManager_->get_time_step();
  caManager_->setup_integration(dT_FE);
  // Time integration
  double time = 0.0;
  int iStep = 0;
  while (time < finalTime_)
  {
    // Detect phase state of each cell in CA based on the current temperature field 
    timer_->start_timing_CA_simulation();
    caManager_->checkup_current_temperature_to_update_phase_state();
    caManager_->checkup_liquid_cell_neigh_to_growth(time);
    // Determine the minimum time step
    double dT_CA = caManager_->calculate_time_step();
    double dT = __min(dT_FE, dT_CA);
    double dT_min;
    MPI_Allreduce(&dT, &dT_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    caManager_->timeStep_ = dT_min;
    timer_->stop_timing_CA_simulation();
    // Time integration in FEM
    timer_->start_timing_FEM_simulation();
    feManager_->setup_time_step_for_coupling(dT_min);
    // Time integration in CA
    timer_->start_timing_CA_simulation();
    caManager_->setup_time_step_for_coupling(dT_min);
    caManager_->single_integration_step();
    // Update information
    time += dT_min;
    iStep++;
      if (!feManager_->pseudo_single_integration_step()) {
        break; // due to reading to the end of the file
      }
      if (!feManager_->solveProblemMyself_) {
          dT_FE = feManager_->get_max_time_step();
          caManager_->setup_integration(dT_FE);
      }
  }
}

/*----------------------------------
    finalize
  ----------------------------------*/
void
CafeSolver::finalize()
{
  // Stop timing
  timer_->stop_timing_entire_modeling();

  // Statistics final grains
  int numberOrientation = 0;
  caManager_->get_amount_of_final_orientations(numberOrientation);
  CafeEnv::self().caOutputP0() << "=====================================" << "\n";
  CafeEnv::self().caOutputP0() << " Final Number of Grains:" << "\n";
  CafeEnv::self().caOutputP0() << "             " << numberOrientation << "\n";
  CafeEnv::self().caOutputP0() << "=====================================" << std::endl; 

  // Output time log
  double maxTimingArray[10];
  double minTimingArray[10];
  MPI_Allreduce(timer_->timingArray_, maxTimingArray, 10,
    MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(timer_->timingArray_, minTimingArray, 10,
    MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  if (feManager_->procID_ == 0)
  {
    timer_->output_timing_info(feManager_->numProcs_,
      maxTimingArray, minTimingArray);
  }

  // output final microstructure information
  if (caManager_->outputMicrostructureInfo_) {
    caManager_->output_microstructure_information();
  }

  // Output final results
  std::string vtkFileName = "finalGrains.vtk";
  caManager_->output_last_step_result();
  caManager_->output_orientations(vtkFileName);
}
