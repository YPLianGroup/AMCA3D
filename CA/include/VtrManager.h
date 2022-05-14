#ifndef VTRMANAGER_H
#define VTRMANAGER_H

// General include
#include <iostream>
#include <sstream>
#include <iomanip>
#include <yaml-cpp/yaml.h>
#include "Simulation.h"
#include "VtkManager.h"
//Forward declaration
class CellularAutomataManager_3D;
class MonteCarloManager_3D;

class VtrManager: public VtkManager
{
public:
  // Constructor
  VtrManager(CellularAutomataManager_3D *caManager);

  // Destructor
  virtual ~VtrManager();

  CellularAutomataManager_3D * caManager_;

  /*------------------------
    Vtr Methods
    --------------------------*/
  void vtr_initialization();
   //vtr_initialization -  Open file stream to .vtr file and write headers
  void vtr_write_coordinates();
   // vtr_write_coordinates - write point coordinates to .vtr file
  void vtr_write_cell_data();
   // vtr_write_cell_data - write cell data to .vtr file
  void vtr_finalization();
   // vtr_finalization - close vtr headers to .vtr file and clos vtr filestream

  /*-----------------------
    Pvtr Methods
    -------------------------*/
  void pvtr_generator();  
   // pvtr_generator - create, write and close a .pvtr file

  /*-----------------------
     Additional Methods
    -------------------------*/
  void initialize();
   // initialize several member data
  void load(const YAML::Node& node);
   // load - load information relevant to vtu output from input file
  bool execute(int numStep, double time);
   // execute - execute necessary commands to write .vtr file and time stamp .pvd file  

private:
  // Streams
  std::ofstream vtrOut_;
  std::ofstream pvtrOut_;

  // extent
  int x1ID;
  int x2ID;
  int y1ID;
  int y2ID;
  int z1ID;
  int z2ID;
  int * x1IDArray;
  int * x2IDArray;
  int * y1IDArray;
  int * y2IDArray;
  int * z1IDArray;
  int * z2IDArray;
  // coordinates offset
  double coordOffset[3];
  // cell infor
  int nxLocGhost_;
  int nyLocGhost_;
  int nzLocGhost_;
  int nxyLocGhost_;
  double h_;

};

#endif
