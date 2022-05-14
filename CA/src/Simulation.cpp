//yaml for parsing
#include <yaml-cpp/yaml.h>
#include "CafeParsing.h"
#include "Solver.h"
#include "CellularAutomataSolver.h"
#include "FiniteElementMethodSolver.h"
#include "CafeSolver.h"
#include "Simulation.h"
#include "CafeEnv.h"
#include "Timer.h"

/*-------------------------------------------------------------
    constructor
  -------------------------------------------------------------*/
Simulation::Simulation(const YAML::Node& root_node)
  :m_root_node(root_node),
   cellularAutomata_(false),
   finiteElementMethod_(false)
{  
  //nothing for now
}

/*-------------------------------------------------------------
   destructor
  -------------------------------------------------------------*/
Simulation::~Simulation()
{
  if (solver_)
    delete solver_;
}

/*-------------------------------------------------------------
   load
  -------------------------------------------------------------*/
void
Simulation::load(const YAML::Node& node)
{
  high_level_banner();
  const YAML::Node * mySolver = node.FindValue("solvers");
  if (mySolver)
  {   
    for (size_t iSolver = 0; iSolver < mySolver->size(); iSolver++)
    {
      std::string name;
      (*mySolver)[iSolver] >> name;
      if (name == "cellular_automata")
      {
        cellularAutomata_ = true;
      }
      if (name == "finite_element_method")
      {
        finiteElementMethod_ = true;
      }
    }
  }
  else
    throw std::runtime_error("parser error: at least one solver should be specified!");
  
if (cellularAutomata_ && finiteElementMethod_)
  {
    solver_ = new CafeSolver();
  }
  else if(finiteElementMethod_)
  {
    solver_ = new FiniteElementMethodSolver();    
  }
  else if (cellularAutomata_)
  {
    solver_ = new CellularAutomataSolver();
  }
 
  solver_->load(node);
  
}

/*-----------------------------------------------------------
    initialize
  -----------------------------------------------------------*/
void
Simulation::initialize()
{
  solver_->initialize();
}

/*-----------------------------------------------------------
    run
  -----------------------------------------------------------*/
void 
Simulation::run()
{
  solver_->execute();
  solver_->finalize();
}

/*------------------------------------------------------------
    high_level_banner
  ------------------------------------------------------------*/
void
Simulation::high_level_banner()
{
  CafeEnv::self().caOutputP0() << std::endl;
  CafeEnv::self().caOutputP0() << "=====================================================================" << std::endl;
  CafeEnv::self().caOutputP0() << "  ______   __       __   ______    ______    ______   _______  \n"
                                  " /      \\ |  \\     /  \\ /      \\  /      \\  /      \\ |       \\ \n"
                                  "|  $$$$$$\\| $$\\   /  $$|  $$$$$$\\|  $$$$$$\\|  $$$$$$\\| $$$$$$$\\\n"
                                  "| $$__| $$| $$$\\ /  $$$| $$   \\$$| $$__| $$ \\$$__| $$| $$  | $$\n"
                                  "| $$    $$| $$$$\\  $$$$| $$      | $$    $$  |     $$| $$  | $$\n"
                                  "| $$$$$$$$| $$\\$$ $$ $$| $$   __ | $$$$$$$$ __\\$$$$$\\| $$  | $$\n"
                                  "| $$  | $$| $$ \\$$$| $$| $$__/  \\| $$  | $$|  \\__| $$| $$__/ $$\n"
                                  "| $$  | $$| $$  \\$ | $$ \\$$    $$| $$  | $$ \\$$    $$| $$    $$\n"
                                  " \\$$   \\$$ \\$$      \\$$  \\$$$$$$  \\$$   \\$$  \\$$$$$$  \\$$$$$$$ \n";
  CafeEnv::self().caOutputP0() << "=====================================================================" << std::endl;
  CafeEnv::self().caOutputP0() << " A grain growing code targeting AM ......                " << std::endl;
  CafeEnv::self().caOutputP0() << std::endl;
  CafeEnv::self().caOutputP0() << "TPLS: Boost, yaml_cpp                              " << std::endl;
  CafeEnv::self().caOutputP0() << std::endl;
  CafeEnv::self().caOutputP0() << "--------------------------------------------------------------------" << std::endl;
  CafeEnv::self().caOutputP0() << "             Prof.  Yanping Lian Group             " << std::endl;
  CafeEnv::self().caOutputP0() << std::endl;
  CafeEnv::self().caOutputP0() << "Major developers:  Dr.   Yanping Lian              " << std::endl;
  CafeEnv::self().caOutputP0() << "                   PhD.  Feiyu Xiong               " << std::endl;
  CafeEnv::self().caOutputP0() << std::endl;
  CafeEnv::self().caOutputP0() << "--------------------------------------------------------------------" << std::endl;
  CafeEnv::self().caOutputP0() << std::endl;
}