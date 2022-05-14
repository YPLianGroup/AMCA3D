#ifndef Solver_h
#define Solver_h

//TPL
#include <boost/program_options.hpp>
#include <yaml-cpp/yaml.h>

// Local includes 
#include "Simulation.h"
#include "Timer.h"

class Solver
{
public:
  Solver();
  virtual ~Solver();

  unsigned int nDim_;
  Timer * timer_;
  virtual void load(const YAML::Node& node);
  virtual void initialize()=0;
  virtual void execute()=0;
  virtual void finalize();
};

#endif
