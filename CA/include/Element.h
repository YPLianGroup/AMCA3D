#ifndef ELEMENT_H
#define ELEMENT_H

/*=======================================================================
                      Element
     Define the element class base for FEM
  =======================================================================*/
#include <vector>

class Element
{
public:
  // Constructor
  Element();

  // Destructor
  virtual ~Element();

  /*-------------------------
      Basic Attributes
    -------------------------*/    
  int* nID_;            // nodal ID    
  bool birth_;          // active state  
  double birthTime_;      // birth time
  bool liquid_;         // phase state    
  double volume_;       // element volume    
  double* stiffMatrix_; // element-level stiffness matrix    
  int matID_;           // material id
  virtual int get_nodes_per_element() = 0;
  /*---------------------------------------
     Basic methods associated with element
     type and geometry
    ---------------------------------------*/
  // Jacobian determinant and inverse Jacobian matrix calculation
  virtual void 
    Jacobian(double* deriv, double* coords, double* iJac, double &detJac) = 0;

  // shape functin of each node at a given point  
  virtual void shape_fcn(double* parCoord, double* shapeFcn) = 0;
    
  // calculate charateristic length of element
  virtual double charateristic_length(double* coordinates) = 0;
    
  // derivative of shape function with respect to physical coordinates
  virtual void 
    derivative_of_shape_fuction_about_real_coords
    (double* deriv, double* iJac, double* gradN) = 0;
  
  // derivative of shape function with respect to natural coordinates
  virtual void 
    derivative_of_shape_function_about_parametic_coords
    (double* parCoord, double* deriv) = 0;

  /*---------------------------------------
    Basic methods associated with FEM
    ---------------------------------------*/
  // stiffness matrix of element level
  virtual void 
    element_level_stiffness_matrix
    (int numIp, double* nodalCoords, double* eleStiffMatrix) = 0;

  // time step of element level  
  virtual double element_level_time_step()=0;    

  // Gaussian quadrature to a given function
  virtual void 
    integration_over_the_element_domain
    (int numIp, double* nodalCoords, double* ipValue, double* nodalValue) = 0;

  // Element volume
  virtual double element_volume(int numIp, double *nodalCoords) = 0;
};

/*=======================================================================
  Hexahedral Element
  Define the 3D element class for FEM
  =======================================================================*/
class HexahedralElement: public Element
{
public:
  HexahedralElement();
  ~HexahedralElement();

  int get_nodes_per_element() { return 8; };
  /*---------------------------------------
     Basic methods associated with element
     type and geometry
    ---------------------------------------*/
  // Jacobian determinant and inverse Jacobian matrix calculation
  void Jacobian(double* deriv, double* coords, 
                double* iJac, double &detJac);
  
  // shape functin of each node at a given point
  void shape_fcn(double* parCoord, double* shapeFcn);
  
  // derivative of shape function with respect to physical coordinates
  void derivative_of_shape_fuction_about_real_coords
       (double* deriv, double* iJac, double* gradN);

  // derivative of shape function with respect to natural coordinates
  void derivative_of_shape_function_about_parametic_coords
       (double* parCoord, double* deriv);

  // calculate charateristic length of element
  double charateristic_length(double* coordinates);
  

  /*---------------------------------------
  Basic methods associated with FEM
  ---------------------------------------*/
  // stiffness matrix of element level
  void element_level_stiffness_matrix
       (int numIp, double* nodalCoords,double* eleStiffMatrix);

  // time step of element level
  double element_level_time_step();

  // Gaussian quadrature to a given function  
  void integration_over_the_element_domain
  (int numIp, double* nodalCoords, double* ipValue, double* nodalValue);
  
  // Element volume
  double element_volume(int numIp, double *nodalCoords);
};

/*=======================================================================
   Quadrilateral Element
   Define the 2D element class for FEM
  =======================================================================*/
class QuadrilateralElement: public Element
{
public:
  QuadrilateralElement();
  ~QuadrilateralElement();

  int get_nodes_per_element() { return 4; };
  /*---------------------------------------
    Basic methods associated with element
    type and geometry
  ---------------------------------------*/
  // Jacobian determinant and inverse Jacobian matrix calculation
  void Jacobian(double* deriv, double* coords,
                double* iJac, double &detJac);

  // shape functin of each node at a given point
  void shape_fcn(double* parCoord, double* shapeFcn);

  // derivative of shape function with respect to physical coordinates
  void derivative_of_shape_fuction_about_real_coords
       (double* deriv, double* iJac, double* gradN);

  // derivative of shape function with respect to natural coordinates
  void derivative_of_shape_function_about_parametic_coords
       (double* parCoord, double* deriv);

  // calculate charateristic length of element
  double charateristic_length(double* coordinates);

  /*---------------------------------------
    Basic methods associated with FEM
    ---------------------------------------*/
  // stiffness matrix of element level
  void element_level_stiffness_matrix
       (int numIp, double* nodalCoords, double* eleStiffMatrix);

  // time step of element level
  double element_level_time_step();

  // Gaussian quadrature to a given function  
  void integration_over_the_element_domain
       (int numIp, double* nodalCoords, double* ipValue, double* nodalValue);
  
  // Element volume
  double element_volume(int numIp, double *nodalCoords);
};
#endif
