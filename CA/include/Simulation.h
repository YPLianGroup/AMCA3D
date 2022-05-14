#ifndef SIMULATION_H
#define SIMULATION_H

namespace YAML
{
  class Node;
}

class Solver;
class CellularAutomataSolver;

class Simulation
{
public:
  Simulation(const YAML::Node &root_node);
  ~Simulation();

  void load(const YAML::Node &node);
  void initialize();
  void run();
  void high_level_banner();

  const YAML::Node &m_root_node;

  Solver * solver_;
  bool cellularAutomata_;
  bool finiteElementMethod_;
};

#endif

