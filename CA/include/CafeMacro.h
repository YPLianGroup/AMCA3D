/*-------------------------------------
   Define common macros and paramters
  -------------------------------------*/
#ifndef CafeMacro_H
#define CafeMacro_H

#include <string>
namespace CAFE {
  // several common variables
  //static const double nullMatTemperature = -1.0e8;
  //static const double M_PI = 3.14159265358979323846;

  /*----------------------------------------------
      For output controls
    ----------------------------------------------*/
  const int nbNames = 13;
  static const std::string variableName[nbNames] = {
    "temperature",   "temperature_rate",
    "grain_ID",      "grain_orientation",
    "grain_velocity", "element_birth",
    "melted_flag",    "phase_ID",
    "alpha_phase",   "phase_orientation",
    "velocity",      "debug"            ,        "debugMC"
  };

  static const std::string 
   variableComponents3D[nbNames] = {"1", "1", "1", "3","3","1","1","1","1","3","3" };

  /*------------------------------------------------
      For Gaussian quadrature
    ------------------------------------------------*/
  // Gaussian quadrature
  const double GaussianWeight[5][5] =
  {
    { 2.0,               0.0,               0.0,               0.0,               0.0 },
    { 1.0,               1.0,               0.0,               0.0,               0.0 },
    { 5.0 / 9.0,         8.0 / 9.0,         5.0 / 9.0,         0.0,               0.0 },
    { 0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454, 0.0 },
    { 0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189 }
  };
  const double GaussianPoint[5][5] =
  {
    { 0.0,                 0.0,               0.0,               0.0,               0.0 },
    { -0.577350269189626,  0.577350269189626, 0.0,               0.0,               0.0 },
    { -0.774596669241483,  0.0,               0.774596669241483, 0.0,               0.0 },
    { -0.86113631159405,  -0.33998104358486,  0.33998104358486,  0.86113631159405,  0.0 },
    { -0.906179845938664, -0.538469310105683, 0.0,               0.538469310105683, 0.906179845938664 }
  };

  /*---------------------------------------------
      For VTK format
    ---------------------------------------------*/
  // Vtk indentation rule
  const std::string SecondEleWidth_ = "  ";
  const std::string ThirdEleWidth_  = "    ";
  const std::string FourthEleWidth_ = "      ";
  const std::string FifthEleWidth_  = "        ";
  const std::string SixthEleWidth_  = "          ";

  /*----------------------------------------------------
                  Hexahedral Element 
      Node numbering convention.
    ----------------------------------------------------*/
  /*
             7 *-------------* 6
              /:            /:
             / :           / : 
            /  :          /  :
           /   :         /   :
         4*-------------* 5  :
          |    *--------|----*        
          |   / 3       |   / 2
          |  /          |  /
          | /           | /
        0 *-------------* 1
  */ 
  const int mappedNode[6][4] = { // following counter-clock-wise order
    {4,5,6,7},{0,3,2,1},         // top and bottom
    {0,1,5,4},{3,7,6,2},         // front and rear
    {0,4,7,3},{1,2,6,5}          // left and right
  };
  // used for the element-level stiffness matrix
  const int mapIndex[8][8] = { { 0,  1,  2,  3,  4,  5,  6,  7 },
  { 1,  8,  9, 10, 11, 12, 13, 14 },
  { 2,  9, 15, 16, 17, 18, 19, 20 },
  { 3, 10, 16, 21, 22, 23, 24, 25 },
  { 4, 11, 17, 22, 26, 27, 28, 29 },
  { 5, 12, 18, 23, 27, 30, 31, 32 },
  { 6, 13, 19, 24, 28, 31, 33, 34 },
  { 7, 14, 20, 25, 29, 32, 34, 35 } };

  /*--------------------------------------
      Methods used for math calculation
    --------------------------------------*/
  // cross product of vector: Vector_c = Vector_a x Vector_b
  void vector_cross_product(double* a, double* b, double* c);

  // calculate vector norm
  double vector_2nd_norm(double *vector, unsigned int n);

  // inner product of vector
  double vector_inner_product(double *vector1, double *vector2, unsigned int n);
  
  // vector normalization
  double vector_normalization(double* v);  

  // calculate the inverse of a third order matrix
  void matrix_inverse(double a[][3], double inverseA[][3], double &determinant);

  // get the transposition of matrix 
  void matrix_transposition(double matrix[][3], double matrixT[][3]);

  // transform coordinates between two cooridnates system
  void coordinate_transformation(double transformMatrix[3][3], double inputVector[3], double outputVector[3]);

}

#endif
