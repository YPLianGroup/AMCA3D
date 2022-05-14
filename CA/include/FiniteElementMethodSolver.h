#ifndef FiniteElementMethodSolver_h
#define FiniteElementMethodSolver_h

// TPL
#include <mpi.h>
#include <math.h>
#include <boost/program_options.hpp>
#include <yaml-cpp/yaml.h>
// Local includes
#include "Simulation.h"
#include "Solver.h"
class FiniteElementManager;
class Material_HeatConduction;

class FiniteElementMethodSolver: public Solver
{
public:
  FiniteElementMethodSolver();
  ~FiniteElementMethodSolver();

  /*-------------------------
     time integration info
    -------------------------*/
  double finalTime_;
  double timeStepFactor_;
  /*-------------------------
     Manage finite elements
    -------------------------*/
  FiniteElementManager * feManager_;

  /*-------------------------
      Material law
    -------------------------*/
  //Material_HeatConduction * matHeatConduction_;
  /*------------------------
       Methods
    -----------------------*/
  void load(const YAML::Node& node);
  
  void initialize();

  void execute();

  void finalize();
};
#endif
