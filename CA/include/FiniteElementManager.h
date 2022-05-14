#ifndef FINITEELEMENTMANAGER_H
#define FINITEELEMENTMANAGER_H

/*=======================================================================
                 FiniteElementManager
   Define the Manager to manage finite element list and finite node list
  =======================================================================*/
#include <vector>
#include <mpi.h>
#include <fstream>
#include "Simulation.h"
#include "ProblemPhysics.h"
#include "Timer.h"
class Node;
class Element;
class HexahedralElement;
class VtuManager;
class ThermalBoundaryCondition;

class FiniteElementManager
{
public:
  // Constructor
  FiniteElementManager();

  // Destructor
  virtual ~FiniteElementManager();

  // time cost record
  Timer * FEMtimer_;
  void set_timer_pointer(Timer *timer) { FEMtimer_ = timer; };

  /*------------------------
     Nodes and Elements
   -------------------------*/
  int nDim_;           // dimension of the problem
  int numNodes_;       // global nodes number
  int numElements_;    // global elements number
  std::vector<Node> nodeVector_;          // global nodes list
  std::vector<Element*> elementVector_;   // global element list
  double* nodalForceArray_;    // global, external and internal nodal force after all reduce  
  double* arrayForAllReduce_;  // global, tricky array using for all reduce operation
  double* thetaArray_;         // global, nodal value (temperature) array

  std::vector<int> grainAroundNodeVector_;
  void set_grain_around_node_vector();

  // local nodes with global ID
  int nodeStart_;            // starting node ID
  int nodeEnd_;              // ending node ID
  int numNodesLocal_;        // local nodes number
    
  // local elements with global ID
  int elementStart_;         // starting element ID
  int elementEnd_;           // ending element ID
  int numElementsLocal_;     // local elements number
      
  /*-------------------------
      time integration info
    -------------------------*/
  double initialTheta_;   // initial condition
  double dT_;             // time step, dT_ : timestep for CA, dT_max: timestep for FEM
  double time_;           // current time
  int iStep_;             // current steps
  double timeStepFactor_; // time step factor

  bool microTimeStepStrategy_;

  /*-----------------
     Parallel info
    -----------------*/
  int numProcs_;
  int procID_;
  // MPI_Comm cartComm_;

  /*---------------
     Output 
    ---------------*/
  friend class VtuManager;
  VtuManager* vtuManager_;
  /*---------------
      Material  
    ---------------*/
  double nullMatTemperature_;

  /*----------------------------------------
     Methods associated with initialization
    ----------------------------------------*/
  void load(const YAML::Node& node);

  void initialize();

  void node_element_list_decomposition();

  /*-----------------------------------------
      Mehtods associated with integration
    -----------------------------------------*/
  void setup_integration(double timeStepFactor);

  void setup_time_step_for_coupling(double dT);

  void time_integration(double finalTime);

  double get_time_step() { return dT_; };
  double get_max_time_step() ;

//  void set_problem_physics(Material_HeatConduction * p) { matHeatConduction_ = p; };

  /*--------------------------------------------
      Methods about loading simulation result 
      from others Software, say Flow-3D
      all the associated member dataa are
      with a prefix of in-
    --------------------------------------------*/
  bool solveProblemMyself_;
  std::string inFileName_;  
  std::ifstream inStreamTracking_;
  double lengthScale_;  // used for length unit conversion
  double timeScale_;    // used for time unit conversion
  double positionOffset_[3]; // location offset
  bool skippingHead_;   // skip several time step at begining if they are the same(heat source is outside the domain)
  
  int inNumTitleLines_;
    // the number of lines as the file-level title description in .txt file
  int inNumSubtitleLines_;
    // the number of lines as the step-level title description in .txt file
  int inStepRecord_;
    // record the number of steps in the .txt file
  double inTimeOffset_;
    // the initial time in the .txt file and serves as the offset in the current simulation
  double inTimeSeries_[2];
    // the time stamps of two adjacent steps in .txt file
  bool pseudo_single_integration_step();
    // get the current result by interpolating from a sequence of two steps results
  void load_mesh_and_first_two_steps_results_from_txt_file();
    // load mesh from the .txt file and first two step's result as well
  bool load_next_step_nodal_data();
    // load next step results from .txt file 
  double get_output_time_step_from_txt_file();
    // get output time step from txt file
};
#endif
