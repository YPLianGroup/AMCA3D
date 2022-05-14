#include "Solver.h"

/*  ~  constructor   */
Solver::Solver()
{
  timer_ = new Timer;
}

/*  ~  destructor   */
Solver::~Solver()
{
  if (timer_) {
    delete timer_;
  }
}

/*  ~ finalize   */
void
Solver::finalize()
{
  //nothing for now
}
void
Solver::load(const YAML::Node &node)
{
  //nothing for now
}
