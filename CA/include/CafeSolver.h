#ifndef CafeSolver_h
#define CafeSolver_h

// TPL
#include <mpi.h>
#include <math.h>
#include <boost/program_options.hpp>
#include <yaml-cpp/yaml.h>
// Local includes
#include "Simulation.h"
#include "Solver.h"

class FiniteElementManager;
class CellularAutomataManager;
class MapVoxelManager;
class ProblemPhysics;

class CafeSolver : public Solver
{
public:
  CafeSolver();
  virtual ~CafeSolver();
  /*  time integration info  */
  double finalTime_;
  double timeStepFactorFE_;
  double maxDt_;
  int maxSteps_;

  /*  Managers  */
  CellularAutomataManager * caManager_;
  FiniteElementManager * feManager_;
  MapVoxelManager * mapVoxelManager_;

  /*  Material Models  */
  ProblemPhysics * problemPhysics;

  /*------------------------
       Interface Methods
  -----------------------*/
  void load(const YAML::Node& node);

  void initialize();

  void execute();

  void finalize();
  /*---------------------------
      Coupling Methods via 
      Time Integration Method
    ---------------------------*/
  // one-way coupling method without subcycle of time step for CA   
  void time_integration_weak_coupling();
};

#endif
