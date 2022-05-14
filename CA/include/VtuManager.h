#ifndef VtuManager_h
#define VtuManager_h
// General include
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <yaml-cpp/yaml.h>
#include "Simulation.h"
#include "VtkManager.h"
//Forward declaration
class FiniteElementManager;

class VtuManager: public VtkManager
{
public:
  // constructor
  VtuManager(FiniteElementManager* feManager);

  // destructor
  virtual ~VtuManager();

  // vtu manager target
  FiniteElementManager* feManager_;

  /*--------------------------
     Vtu Methods
    --------------------------*/
  void vtu_initialization();
  //vtu_initialization -  Open file stream to .vtu file and write headers
  void vtu_write_coordinates();
  // vtu_write_coordinates - write point coordinates to .vtu file
  void vtu_write_cell_data_scalar();
  // vtu_write_cell_data_scalar - write scalar values on cell to .vtu file
  void vtu_finalization();
  // vtu_finalization - close vtu headers to .vtu file and clos vtu filestream
  void vtu_write_point_data();
  // vtu_write_scalar -  write scalar values on node to .vtu file
  void vtu_set_connect();

  /*-----------------------
    Pvtr Methods
    -------------------------*/ 
  void pvtu_generator();
  // pvtu_generator - create, write and close a .pvtu file

  /*-----------------------
  Additional Methods
  -------------------------*/
  void initialize();
  // initialize several member data
  void load(const YAML::Node& node);
  // load - load information relevant to vtu output from input file
  void execute(int numStep, double time);
  // execute - execute necessary commands to write .vtr file and time stamp .pvd file
  
private:
  // Streams
  std::ofstream vtuOut_;
  std::ofstream pvtuOut_;

  
  // point/cell infor
  int numPoints_;
  int numCells_;
  int numCellsOutput_;
  
};

#endif
