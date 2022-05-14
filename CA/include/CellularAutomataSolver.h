#ifndef CaSolver_H
#define CaSolver_H

// TPL
#include <mpi.h>
#include <math.h>
#include <boost/program_options.hpp>
#include <yaml-cpp/yaml.h>
// Local includes
#include "Simulation.h"
#include "Solver.h"

class CellularAutomataManager;
class ProblemPhysics;

class CellularAutomataSolver: public Solver
{
public:
  CellularAutomataSolver();
  virtual ~CellularAutomataSolver();

  /*-------------------------
     high-level data member
    -------------------------*/
  CellularAutomataManager * caManager_;
  ProblemPhysics * problemPhysics_;

  /*------------------------------
    time integration info
    (will move to other class)
    -------------------------------*/
  double maxDt_;
  double finalTime_;
  int maxSteps_;

  /*---------------------
      Interface methods
    ---------------------*/
  void load(const YAML::Node &node);
  void initialize();
  void execute();

  /*------------------------
      output final result
    ------------------------*/
  void finalize();

};
#endif
