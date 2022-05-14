#ifndef MapVoxelManager_h
#define MapVoxelManager_h
#include <vector>
#include "Simulation.h"

class Bin {
public:
  Bin();
  ~Bin();

  std::vector <int> elements_;
  std::vector< double > birthTime_;
  std::vector< std::vector<double> > eleBounds_;
  std::vector< std::vector< std::vector<double> > > nodes_;
  std::vector< std::vector< int>  > nodesID_;
  bool hasElements_;

};



/*------------------------------------------------------
                   MapVoxelManager
  Used to map the CA network to corresponding FEM mesh 
  and interporlate the temperature of CA cell fromd FEM
  solution
 -------------------------------------------------------*/
class Element;
class FiniteElementManager;
class CellularAutomataManager;
class CellularAutomataManager_3D;

class MapVoxelManager
{
public:
  // constructor
  MapVoxelManager(FiniteElementManager* feManager, CellularAutomataManager* caManager);
  // destructor
  virtual ~MapVoxelManager();

  double nullMatTemperature_;    // temperature indicates void
  /*-------------------------
     Two managers
    -------------------------*/
  FiniteElementManager* feManager_;
  CellularAutomataManager* caManager_;

  /*---------------------------
      domain size info
    ---------------------------*/
  int nDim_;
  int numVoxels_;
    // number of voxels from CA netwrok
  double xyzMinMax_[6];
    // size of the FE mesh domain
  int sizeOfBin[3];
    // the size of bin in 3 directions in terms of CM cell size
  
  /*------------------------------
     Lists for linking 
    ------------------------------*/
  int* elementTagList_;  
    // element tag for the CA cell: plus means ownning an active element
    //                              minus means with inactive or without element
  double* parCoordsList_;
    // parametric coordinates of the CA cell
  
  /*----------------------
      Interface Methods  
    ----------------------*/
  void load(const YAML::Node& node);
    // load the parameters from input file
  void initialize();
    // initialize all the lists
  void execute();
  void execute_with_focus_on_CA_region();
    // fill up the lists
  void update();
    // update some part of the lists if necessary
  bool evaluate_temperature(int voxelID, double &interpolatinValue);
    // get the temperature from heat conduction solver

  /*--------------------------
     Methods used for mapping
    --------------------------*/    
  void get_domain_xyz_min_max();
  void get_CA_domain_xyz_min_max();  
    // get the cubioc domain that cover the mesh domain
  void get_cell_xyz_min_max(double* range, Element* ele);
    // get the cubioc domain that cover the cell 
  void attach_elements(double* dBin, int* numBin, std::vector<Bin> &binList);
    // attach elements to each bin
  void tag_voxels(double* dBin, int* numBin, std::vector<Bin> &binList);
    // tag voxel to each bin

  /*--------------------------------
     Methods about Newton-Raphson iteration
    --------------------------------*/
  int maxNonlinearIter_;
  double isInElemTolerance_;
  bool is_in_element(Element* ele, double *pointCoords, double* paraCoords);
    // check if the given poing is located in the element
  bool within_tolerance(const double &val, const double &tol);
    // check up error
  double vector_norm_sq(const double * vect, int len);
    // get L2 norm
  double parametric_distance(std::vector<double> x);
};

#endif
