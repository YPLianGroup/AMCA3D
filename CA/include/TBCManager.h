/*==============================================
       Thermal Boundary Condition Manager
  ==============================================*/
#ifndef ThermalBoundaryConditionManager_H
#define ThermalBoundaryConditionManager_H

#include <vector>
#include <list>
#include "ThermalBoundaryCondition.h"
#include "yaml-cpp/yaml.h"
struct Surface {
public:
  int elementID_;       // the element with the surface
  int faceID_;          // the face of the element
  
  double birthTime_;     
  double deathTime_;

  bool convectionBc_;  // convection B.C.
  bool radiationBc_;   // radiation B.C.
  bool movingFluxBc_;  // moving Flux 
  bool fluxBc_;
};

class Element;
class FiniteElementManager;

class TBCManager
{
public:
  // constructor
  TBCManager(FiniteElementManager * feManager);

  // destructor
  virtual ~TBCManager();

  // element birth and death technology
  bool birthDeathOption_;

  // target
  FiniteElementManager * feManager_;

  // current time
  double currentTime_;

  // flux boundary conditions vector
  std::vector<ThermalBoundaryCondition*> tbcVector_;
  
  // moving flux boundary condition
  std::vector<MovingFluxBoundaryCondition *> movingFluxVector_;

  // surfaces that will have b.c.s for a moment or forever
  std::list<Surface *> surfaceList_;
  
  // manximum element ID in the surfaceList_
  int maximumEleID_;

  // number of Gaussian points along each direction
  int numIp_;

  /*----------------------------------
      Interface methods
    ----------------------------------*/
  // load parameters from input file
  void load(const YAML::Node& node);

  // initialize
  void initialize();

  // execute boundary condtions calculation
  void execute(double time);

  /*----------------------------------------
     Additional methods
    ----------------------------------------*/
  // impose a bunch of flux boundary conditions
  void impose_boundary_condition(Surface* surf, Element* ele);

  // calculate the co-rotation coordinates systme
  void map_3D_coordinates_to_2D_one(double xN3D[4][3], double xN3Dlocal[4][3], bool additionalWork,
                                    double* gCoords, double *lCoords);
  // calculate the coordinate transformation matrix 
  void calculate_coordinate_transofrmation_matrix(double xN3D[4][3], double g2LTransfMatrix[3][3]);
  
  // read in surface data from file and take the ones belong to the current process
  void load_surface_info_to_get_my_portion(std::string fileName);
};

#endif