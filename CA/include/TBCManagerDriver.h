#ifndef TBCManagerDriver_H
#define TBCManagerDriver_H

#include <vector>
#include "TBCManager.h"
#include "yaml-cpp/yaml.h"

class FiniteElementManager;

class TBCManagerDriver
{
public:
  // constructor
  TBCManagerDriver();

  // destructor
  virtual ~TBCManagerDriver();

  // vector of TBCManager
  std::vector<TBCManager> tbcManagerVector_;
  // number of boundary condition managers
  int numBcM_;

  /*------------------------------
      Interface methods
    ------------------------------*/
  // load parameters from input file
  void load(const YAML::Node& node, FiniteElementManager * feManager);

  // initialize
  void initialize();

  // execute boundary condtions calculation
  void execute(double time);
};

#endif
