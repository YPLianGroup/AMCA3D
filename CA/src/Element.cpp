/*=======================================================================
                        Element
              Definition of Element base class.
  =======================================================================*/
#include "Element.h"
#include <assert.h>
#include <algorithm>
#include <cmath>
#include "CafeMacro.h"
#define __max(a,b)  (((a) > (b)) ? (a) : (b))
#define __min(a,b)  (((a) < (b)) ? (a) : (b))

/*-----------------
    Constructor
  -----------------*/
Element::Element()
  : birth_(true),     // this has to be true
    birthTime_(-1.0),
    liquid_(false),
    matID_(0),
    nID_(NULL),
    stiffMatrix_(NULL)
{
  //nothing for now
}

/*-----------------
    Destructor
  -----------------*/
Element::~Element()
{
  // nothing for now
}


/*=======================================================================
                           Hexahedral Element
                   Declaration of Element base class.
  =======================================================================*/

/*---------------
   - Constructor
  ---------------*/
HexahedralElement::HexahedralElement()
  :Element()
{
  nID_ = new int[8];
  stiffMatrix_ = NULL;
}

/*--------------
   - Destructor
  --------------*/
HexahedralElement::~HexahedralElement()
{
  if (nID_)
    delete[] nID_;
  if (stiffMatrix_)
    delete[] stiffMatrix_;
}

/*-------------------------
   - Jacobian
     coords[24] = [8][3]
     deriv[24] = [3][8];
  -------------------------*/
void
HexahedralElement::
Jacobian(double* deriv, double* coords, double* iJac, double &detJac)
{
  int count = 0;
  double Jac[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      Jac[i][j] = 0.0;

  for (int k = 0; k < 3; k++) {
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 8; i++) {
        Jac[k][j] += deriv[k * 8 + i] * coords[i * 3 + j];
      }
    }
  }

  // calculate determinant
  detJac = Jac[0][0] * Jac[1][1] * Jac[2][2] + Jac[0][1] * Jac[1][2] * Jac[2][0] +
           Jac[1][0] * Jac[2][1] * Jac[0][2] - Jac[0][2] * Jac[1][1] * Jac[2][0] -
           Jac[0][0] * Jac[1][2] * Jac[2][1] - Jac[0][1] * Jac[1][0] * Jac[2][2];
  assert(detJac >= 0.0);

  // cofactor matrix of Jacobian matrix
  iJac[0] = Jac[1][1] * Jac[2][2] - Jac[1][2] * Jac[2][1];
  iJac[1] = Jac[1][0] * Jac[2][2] - Jac[1][2] * Jac[2][0];
  iJac[2] = Jac[1][0] * Jac[2][1] - Jac[1][1] * Jac[2][0];

  iJac[3] = Jac[0][1] * Jac[2][2] - Jac[0][2] * Jac[2][1];
  iJac[4] = Jac[0][0] * Jac[2][2] - Jac[0][2] * Jac[2][0];
  iJac[5] = Jac[0][0] * Jac[2][1] - Jac[0][1] * Jac[2][0];

  iJac[6] = Jac[0][1] * Jac[1][2] - Jac[0][2] * Jac[1][1];
  iJac[7] = Jac[0][0] * Jac[1][2] - Jac[0][2] * Jac[1][0];
  iJac[8] = Jac[0][0] * Jac[1][1] - Jac[0][1] * Jac[1][0];

  for (unsigned char i = 0; i < 3; i++) {
    for (unsigned char j = 0; j < 3; j++) {
      iJac[i * 3 + j] *= std::pow(-1.0, i + j) / detJac;
    }
  }
  // inverse matrix of Jacobian matrix
  double swap;
  swap = iJac[1];
  iJac[1] = iJac[3];
  iJac[3] = swap;

  swap = iJac[2];
  iJac[2] = iJac[6];
  iJac[6] = swap;

  swap = iJac[5];
  iJac[5] = iJac[7];
  iJac[7] = swap;
}
/*----------------------------------------------------------
   - derivative_of_shape_function_about_real_coords
   gradN[24]=[3][8]
  ----------------------------------------------------------*/
void
HexahedralElement::derivative_of_shape_fuction_about_real_coords
(double* deriv, double* iJac, double* gradN)
{
  for (int k = 0; k < 8; k++)
    for (int j = 0; j < 3; j++)
    {
      for (int i = 0; i < 3; i++)
        gradN[j*8 + k] += iJac[j*3 + i] * deriv[i*8 + k];
    }
}

/*--------------------------------------------------------
  - gradient_of_shape_function_about_parametic_coordinates
  -------------------------------------------------------*/
void
HexahedralElement::derivative_of_shape_function_about_parametic_coords
(double* parCoord, double* deriv)
{
  double oneMinusChsi = 1.0 - parCoord[0];
  double onePlusChsi  = 1.0 + parCoord[0];
  double oneMinusEta  = 1.0 - parCoord[1];
  double onePlusEta   = 1.0 + parCoord[1];
  double oneMinusZeta = 1.0 - parCoord[2];
  double onePlusZeta  = 1.0 + parCoord[2];
  
  // with respect to chsi
  deriv[0] = -0.1250 * oneMinusEta * oneMinusZeta;
  deriv[2] =  0.1250 * onePlusEta * oneMinusZeta;
  deriv[4] = -0.1250 * oneMinusEta * onePlusZeta;
  deriv[6] =  0.1250 * onePlusEta * onePlusZeta;
  for (int i = 0; i < 4; i++)
    deriv[i * 2 + 1] = -deriv[i * 2];

  // with respect to eta
  deriv[0 + 8] = -0.1250 * oneMinusChsi * oneMinusZeta;
  deriv[3 + 8] = -deriv[8];
  deriv[1 + 8] = -0.1250 * onePlusChsi * oneMinusZeta;
  deriv[2 + 8] = -deriv[9];
  deriv[4 + 8] = -0.1250 * oneMinusChsi * onePlusZeta;
  deriv[7 + 8] = -deriv[12];
  deriv[5 + 8] = -0.1250 * onePlusChsi * onePlusZeta;
  deriv[6 + 8] = -deriv[13];

  // with respect to zeta
  deriv[4 + 16] = 0.1250 * oneMinusChsi * oneMinusEta;
  deriv[5 + 16] = 0.1250 * onePlusChsi * oneMinusEta;
  deriv[6 + 16] = 0.1250 * onePlusChsi * onePlusEta;
  deriv[7 + 16] = 0.1250 * oneMinusChsi * onePlusEta;
  for (int i = 0; i < 4; i++)
    deriv[i + 16] = -deriv[i + 20];
  
}

/*---------------
    shape_fcn
  ---------------*/
void
HexahedralElement::shape_fcn(double* parCoord, double* N)
{
  double chsi = parCoord[0];
  double eta = parCoord[1];
  double zeta = parCoord[2];

  N[0] = 0.125*(1.0 - chsi)*(1.0 - eta)*(1.0 - zeta);
  N[3] = 0.125*(1.0 - chsi)*(1.0 + eta)*(1.0 - zeta);
  N[2] = 0.125*(1.0 + chsi)*(1.0 + eta)*(1.0 - zeta);
  N[1] = 0.125*(1.0 + chsi)*(1.0 - eta)*(1.0 - zeta);
  N[4] = 0.125*(1.0 - chsi)*(1.0 - eta)*(1.0 + zeta);
  N[7] = 0.125*(1.0 - chsi)*(1.0 + eta)*(1.0 + zeta);
  N[6] = 0.125*(1.0 + chsi)*(1.0 + eta)*(1.0 + zeta);
  N[5] = 0.125*(1.0 + chsi)*(1.0 - eta)*(1.0 + zeta);
}

/*------------------------
   characteristic_length
  ------------------------*/
double
HexahedralElement::charateristic_length(double* coordinates)
{
  double characterLength;
  // FIXME: replace xyz by coordinates
  double xyz[3][8];
  for (int i = 0; i < 8; i++)
    for (int j = 0; j < 3; j++)
    {
      xyz[j][i] = coordinates[i * 3 + j];
    }
  unsigned char fac[4][6] = { { 0,4,0,1,2,3 },{ 1,5,1,2,3,0 },
                              { 2,6,5,6,7,4 },{ 3,7,4,5,6,7 } };
  double dt, at;
  double e, f, g, atest, x13[3], x24[3], fs[3], ft[3];
  double areal = 1.0e20, aream = 0.0;
  unsigned char k1, k2, k3, k4;

  // Calculate area of each surface
  for (unsigned char j = 0; j<6; j++)
  {
    for (unsigned char i = 0; i<3; i++)
    {
      k1 = fac[0][j];
      k2 = fac[1][j];
      k3 = fac[2][j];
      k4 = fac[3][j];

      x13[i] = xyz[i][k3] - xyz[i][k1];
      x24[i] = xyz[i][k4] - xyz[i][k2];

      fs[i] = x13[i] - x24[i];
      ft[i] = x13[i] + x24[i];
    }// end loop i

    e = fs[0] * fs[0] + fs[1] * fs[1] + fs[2] * fs[2];
    f = fs[0] * ft[0] + fs[1] * ft[1] + fs[2] * ft[2];
    g = ft[0] * ft[0] + ft[1] * ft[1] + ft[2] * ft[2];

    atest = e*g - f*f;     // area/4  (4*area)^2

    aream = __max(atest, aream);
  }// end loop j

  characterLength = 4 * (volume_) / std::sqrt(aream);

  return characterLength;
}

/*----------------------------------
   element_level_stiffness_matrix
  ----------------------------------*/
void
HexahedralElement::
element_level_stiffness_matrix(int numIp, double* coordinates,double* stiffMatrix)
{
  // natural coordinates of Gaussian point
  double parCoords[3];
  // weight of Gaussian points
  double weight, weightK, weightJ, weightI;
  // nodal shape function derivative at Gaussian points, taking the form of [3][8]
  double deriv[24], gradN[24];
  // nodal shape function at Gaussian points
  double shapeFunc[8];
  // determinat of Jacobian matrix at Gaussian point
  double detJac = 0.0;
  // inverse Jacobian matrix
  double iJac[9];

  // loop over all the Gaussian points
  for (int kDim = 0; kDim < numIp; kDim++) {
    weightK = CAFE::GaussianWeight[numIp - 1][kDim];
    parCoords[2] = CAFE::GaussianPoint[numIp - 1][kDim];

    for (int jDim = 0; jDim < numIp; jDim++) {
      weightJ = CAFE::GaussianWeight[numIp - 1][jDim];
      parCoords[1] = CAFE::GaussianPoint[numIp - 1][jDim];
     
      for (int iDim = 0; iDim < numIp; iDim++) {
        weightI = CAFE::GaussianWeight[numIp - 1][iDim];
        parCoords[0] = CAFE::GaussianPoint[numIp - 1][iDim];

        weight = weightK*weightJ*weightI;
        detJac = 0.0;
        for (int i = 0; i < 24; i++)
        {
          deriv[i] = 0.0;
          gradN[i] = 0.0;
          if (i<9)
            iJac[i] = 0.0;
        }
        // get derivative of shape function with respect to parametic coordinates
        derivative_of_shape_function_about_parametic_coords(parCoords, deriv);
        // get inverse Jacob matrix and its determinant
        Jacobian(deriv, coordinates, iJac, detJac);

        // get derivative of shape functin with respect to real coordinates
        derivative_of_shape_fuction_about_real_coords(deriv, iJac, gradN);

        // assemble the stiffness matrix of element
        int count = 0;
        for (int I = 0; I < 8; I++) {
          for (int J = I; J < 8; J++) {
            for (int i = 0; i < 3; i++) {
              stiffMatrix[count] += gradN[i * 8 + I] * gradN[i * 8 + J] * detJac * weight;
            }
            count++;
          }
        } // end I
      } // end iDim
    } // end jDim
  } // end kDim  
}

/*--------------------------------------------------------------------
   integration_over_the_element_domain
   Input: 
     numIp       -  number of integration points in each direction
     nodalCoords -  nodal coordinates
     ipValue     -  functional value at integration points
   Output:
     nodalValue  -  nodal value of the integration
  --------------------------------------------------------------------*/
void
HexahedralElement::integration_over_the_element_domain
(int numIp, double* nodalCoords, double* ipValue, double* nodalValue)
{
  // natural coordinates of Gaussian point
  double parCoords[3];
  // weight of Gaussian points
  double weightK, weightJ, weightI, weight;
  // nodal shape function derivative at Gaussian points, taking the form of [3][8]
  double deriv[24];
  // nodal shape function at Gaussian points
  double shapeFunc[8];
  // determinat of Jacobian matrix at Gaussian point
  double detJac = 0.0;
  // inverse Jacobian matrix
  double iJac[9];
  // ID of the Gaussian points
  int ipID = 0;

  // Confirm 0 intial value in nodalValue
  for (int i = 0; i < 8; i++) {
    nodalValue[i] = 0;
  }

  // loop over all the Gaussian points
  for (int kDim = 0; kDim < numIp; kDim++) {
    weightK = CAFE::GaussianWeight[numIp - 1][kDim];
    parCoords[2] = CAFE::GaussianPoint[numIp - 1][kDim];

    for (int jDim = 0; jDim < numIp; jDim++) {
      weightJ = CAFE::GaussianWeight[numIp - 1][jDim];
      parCoords[1] = CAFE::GaussianPoint[numIp - 1][jDim];

      for (int iDim = 0; iDim < numIp; iDim++) {
        weightI = CAFE::GaussianWeight[numIp - 1][iDim];
        parCoords[0] = CAFE::GaussianPoint[numIp - 1][iDim];

        weight = weightK*weightJ*weightI;
        for (int i = 0; i < 24; i++) {
          deriv[i] = 0.0;
          if (i < 9) {
            iJac[i] = 0.0;           
          }
          if (i < 8) {
            shapeFunc[i] = 0.0;
          }
        }

        // get derivative of shape function with respect to parametic coordinates
        derivative_of_shape_function_about_parametic_coords(parCoords, deriv);

        // get inverse Jacob matrix and its determinant
        Jacobian(deriv, nodalCoords, iJac, detJac);

        // shape function value of each node
        shape_fcn(parCoords, shapeFunc);

        for (int iNode = 0; iNode < 8; iNode++) {
          nodalValue[iNode] += shapeFunc[iNode] * ipValue[ipID] * detJac * weight;
        }

        ipID++;
      }// end iDim
    } // end jDim 
  }
}
/*--------------------------------------------------------------------
   integration_over_the_element_domain
   Input:
     numIp       -  number of integration points in each direction
     nodalCoords -  nodal coordinates

   Output:
     element volume
  --------------------------------------------------------------------*/
double
HexahedralElement::element_volume(int numIp, double *nodalCoords)
{
  double volume = 0.0;
  // natural coordinates of Gaussian point
  double parCoords[3];
  // weight of Gaussian points
  double weightK, weightJ, weightI, weight;
  // nodal shape function derivative at Gaussian points, taking the form of [3][8]
  double deriv[24];
  // determinat of Jacobian matrix at Gaussian point
  double detJac = 0.0;
  // inverse Jacobian matrix
  double iJac[9];

  // loop over all the Gaussian points
  for (int kDim = 0; kDim < numIp; kDim++) {
    weightK = CAFE::GaussianWeight[numIp - 1][kDim];
    parCoords[2] = CAFE::GaussianPoint[numIp - 1][kDim];

    for (int jDim = 0; jDim < numIp; jDim++) {
      weightJ = CAFE::GaussianWeight[numIp - 1][jDim];
      parCoords[1] = CAFE::GaussianPoint[numIp - 1][jDim];

      for (int iDim = 0; iDim < numIp; iDim++) {
        weightI = CAFE::GaussianWeight[numIp - 1][iDim];
        parCoords[0] = CAFE::GaussianPoint[numIp - 1][iDim];

        weight = weightK*weightJ*weightI;
        for (int i = 0; i < 24; i++) {
          deriv[i] = 0.0;
          if (i < 9) {
            iJac[i] = 0.0;           
          }
        }

        // get derivative of shape function with respect to parametic coordinates
        derivative_of_shape_function_about_parametic_coords(parCoords, deriv);

        // get inverse Jacob matrix and its determinant
        Jacobian(deriv, nodalCoords, iJac, detJac);
        
        volume +=  detJac * weight;
      }
    }// end iDim
  } // end jDim

  assert(volume>0.0);
  return volume;
}



/*--------------------------
   element_level_time_step
  --------------------------*/
double
HexahedralElement::element_level_time_step()
{
  double dT;
  // FIXME: element with different \rho*Cp
  return dT;
}

/*=======================================================================
                       Quadrilateral Element
  Declaration of class.
  =======================================================================*/

/*---------------
  - Constructor
  ---------------*/
QuadrilateralElement::QuadrilateralElement()
  :Element()
{
  //nID_ = new int[4];
  //stiffMatrix_ = NULL;
}

/*--------------
  - Destructor
  --------------*/
QuadrilateralElement::~QuadrilateralElement()
{
  if (nID_)
    delete nID_;
  if (stiffMatrix_)
    delete stiffMatrix_;
}

/*-----------------------
  - Jacobian
     coords[8]=[4][2]
     deriv[8] = [2][4]
-------------------------*/
void
QuadrilateralElement::
Jacobian(double* deriv, double* coords, double* iJac, double &detJac)
{
  int count = 0;
  double Jac[2][2];
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      Jac[i][j] = 0.0;

  for (int k = 0; k < 2; k++)
    for (int j = 0; j < 2; j++)
    {
      for (int i = 0; i < 4; i++)
        Jac[k][j] += deriv[k * 4 + i] * coords[i * 2 + j];
    }

  // calculate determinant
  detJac = Jac[0][0] * Jac[1][1] - Jac[0][1] * Jac[1][0];
  assert(detJac >= 0.0);

  // inverse matrix
  iJac[0] =   Jac[1][1] / detJac;
  iJac[1] = - Jac[0][1] / detJac;
  
  iJac[2] = - Jac[1][0] / detJac;
  iJac[3] =   Jac[0][0] / detJac;
  
}
/*----------------------------------------------------------
  - derivative_of_shape_function_about_real_coords
     gradN[8]=[2][4]
  ----------------------------------------------------------*/
void
QuadrilateralElement::derivative_of_shape_fuction_about_real_coords
(double* deriv, double* iJac, double* gradN)
{
  for (int k = 0; k < 4; k++) { // k: node ID
    for (int j = 0; j < 2; j++) { // j: dimension
      for (int i = 0; i < 2; i++) { 
        gradN[j * 4 + k] += iJac[j * 2 + i] * deriv[i * 4 + k];
      }
    }
  }
}

/*------------------------------------------------------------
  - gradient_of_shape_function_about_parametic_coordinates
  ------------------------------------------------------------*/
void
QuadrilateralElement::derivative_of_shape_function_about_parametic_coords
(double* parCoord, double* deriv)
{
  double oneMinusChsi = 1.0 - parCoord[0];
  double onePlusChsi = 1.0 + parCoord[0];
  double oneMinusEta = 1.0 - parCoord[1];
  double onePlusEta = 1.0 + parCoord[1];  

  // with respect to chsi
  deriv[0] = -0.250 * oneMinusEta;
  deriv[2] =  0.250 * onePlusEta ;  
  for (int i = 0; i < 2; i++)
    deriv[i * 2 + 1] = -deriv[i * 2];

  // with respect to eta
  deriv[4] = -0.250 * oneMinusChsi;
  deriv[7] = -deriv[4];
  deriv[5] = -0.250 * onePlusChsi;
  deriv[6] = -deriv[5];
  
}

/*---------------
   shape_fcn
  ---------------*/
void
QuadrilateralElement::shape_fcn(double* parCoord, double* N)
{
  double oneMinusChsi = 1.0 - parCoord[0];
  double onePlusChsi  = 1.0 + parCoord[0];
  double oneMinusEta  = 1.0 - parCoord[1];
  double onePlusEta   = 1.0 + parCoord[1];

  N[0] = 0.25*(oneMinusChsi)*(oneMinusEta);
  N[1] = 0.25*(onePlusChsi) *(oneMinusEta);
  N[2] = 0.25*(onePlusChsi) *(onePlusEta);
  N[3] = 0.25*(oneMinusChsi)*(onePlusEta);
}

/*------------------------
   characteristic_length
  ------------------------*/
double
QuadrilateralElement::charateristic_length(double* coordinates)
{
  double characterLength;
  // to be continued
  return characterLength;
}

/*----------------------------------
    element_level_stiffness_matrix
  ----------------------------------*/
void
QuadrilateralElement::
element_level_stiffness_matrix(int numIp, double* coordinates, double* stiffMatrix)
{
  // natural coordinates of Gaussian point
  double parCoords[2];
  // weight of Gaussian points
  double weightJ, weightI, weight;
  // nodal shape function derivative at Gaussian points, taking the form of [2][4]
  double deriv[8], gradN[8];
  // nodal shape function at Gaussian points
  double shapeFunc[4];
  // determinat of Jacobian matrix at Gaussian point
  double detJac;
  // inverse Jacobian matrix
  double iJac[4];
  
  // loop over all the Gaussian points
  for (int jDim = 0; jDim < numIp; jDim++) {
    weightJ = CAFE::GaussianWeight[numIp - 1][jDim];
    parCoords[1] = CAFE::GaussianPoint[numIp - 1][jDim];
  
    for (int iDim = 0; iDim < numIp; iDim++) {
      weightI = CAFE::GaussianWeight[numIp - 1][iDim];      
      parCoords[0] = CAFE::GaussianPoint[numIp - 1][iDim];
    
      detJac = 0.0;
      weight = weightJ*weightI;
      for (int i = 0; i < 8; i++) {
        deriv[i] = 0.0;
        gradN[i] = 0.0;
        if (i<4)
          iJac[i] = 0.0;
      }    

      // get derivative of shape function with respect to parametic coordinates
      derivative_of_shape_function_about_parametic_coords(parCoords, deriv);
      
      // get inverse Jacob matrix and its determinant
      Jacobian(deriv, coordinates, iJac, detJac);

      // get derivative of shape functin with respect to real coordinates
      derivative_of_shape_fuction_about_real_coords(deriv, iJac, gradN);

      // assemble the stiffness matrix of element
      int count = 0;
      for (int I = 0; I<4; I++) {
        for (int J = I; J < 4; J++) {
          for (int i = 0; i < 2; i++) {
            stiffMatrix[count] += gradN[i * 4 + I] * gradN[i * 4 + J] * detJac*weight;
          }
          count++;
        }
      } // end I
    }// end iDim
  } // end jDim 
}

/*----------------------------------------------------------------------
   integration_over_the_element_domain
   Input:
     numIp       -  number of integration points in each direction
     nodalCoords -  nodal coordinates, taking the form of [4][2]
     ipValue     -  functional value at integration points and shoule 
                    be organized in a specified order
   Output:
     nodalValue  -  nodal value of the integration
  ----------------------------------------------------------------------*/
void
QuadrilateralElement::integration_over_the_element_domain
(int numIp, double* nodalCoords, double* ipValue, double* nodalValue)
{  
  // natural coordinates of Gaussian point
  double parCoords[2];  
  // weight of Gaussian points
  double weightJ, weightI, weight;
  // nodal shape function derivative at Gaussian points, taking the form of [2][4]
  double deriv[8];      
  // nodal shape function at Gaussian points
  double shapeFunc[4];  
  // determinat of Jacobian matrix at Gaussian point
  double detJac = 0.0;  
  // inverse Jacobian matrix
  double iJac[4];
  // ID of the Gaussian points
  int ipID = 0;

  // Confirm 0 intial value in nodalValue
  for (int i = 0; i < 4; i++) {
    nodalValue[i] = 0;
  }

  // loop over all the Gaussian points
  for (int jDim = 0; jDim < numIp; jDim++) {
    weightJ = CAFE::GaussianWeight[numIp - 1][jDim];
    parCoords[1] = CAFE::GaussianPoint[numIp - 1][jDim];

    for (int iDim = 0; iDim < numIp; iDim++) {
      weightI = CAFE::GaussianWeight[numIp - 1][iDim];      
      parCoords[0] = CAFE::GaussianPoint[numIp - 1][iDim];      

      detJac = 0.0;
      weight = weightJ*weightI;
      for (int i = 0; i < 8; i++) {
        deriv[i] = 0.0;        
        if (i < 4) {
          iJac[i] = 0.0;
          shapeFunc[i] = 0.0;
        }
      }

      // get derivative of shape function with respect to parametic coordinates
      derivative_of_shape_function_about_parametic_coords(parCoords, deriv);

      // get inverse Jacob matrix and its determinant
      Jacobian(deriv, nodalCoords, iJac, detJac);      
      
      // shape function value of each node
      shape_fcn(parCoords, shapeFunc);

      for (int iNode = 0; iNode < 4; iNode++) {
        nodalValue[iNode] += shapeFunc[iNode] * ipValue[ipID] * detJac * weight;
      }

      ipID++;
    }// end iDim
  } // end jDim 
}

/*--------------------------
   element_level_time_step
  --------------------------*/
double
QuadrilateralElement::element_level_time_step()
{
  double dT;
  // FIXME: element with different \rho*Cp
  return dT;
}

/*--------------------------------------------------------------------
  integration_over_the_element_domain
  Input:
    numIp       -  number of integration points in each direction
    nodalCoords -  nodal coordinates

  Output:
    element area
  --------------------------------------------------------------------*/
double
QuadrilateralElement::element_volume(int numIp, double *nodalCoords)
{
  double area = 0.0;
  // natural coordinates of Gaussian point
  double parCoords[2];
  // weight of Gaussian points
  double weightJ, weightI, weight;
  // nodal shape function derivative at Gaussian points, taking the form of [2][4]
  double deriv[8];
  // determinat of Jacobian matrix at Gaussian point
  double detJac = 0.0;
  // inverse Jacobian matrix
  double iJac[4];

  // loop over all the Gaussian points
  for (int jDim = 0; jDim < numIp; jDim++) {
    weightJ = CAFE::GaussianWeight[numIp - 1][jDim];
    parCoords[1] = CAFE::GaussianPoint[numIp - 1][jDim];

    for (int iDim = 0; iDim < numIp; iDim++) {
      weightI = CAFE::GaussianWeight[numIp - 1][iDim];
      parCoords[0] = CAFE::GaussianPoint[numIp - 1][iDim];

      detJac = 0.0;
      weight = weightJ*weightI;
      for (int i = 0; i < 8; i++) {
        deriv[i] = 0.0;
        if (i < 4) {
          iJac[i] = 0.0;         
        }
      }

      // get derivative of shape function with respect to parametic coordinates
      derivative_of_shape_function_about_parametic_coords(parCoords, deriv);

      // get inverse Jacob matrix and its determinant
      Jacobian(deriv, nodalCoords, iJac, detJac);

      area += detJac*weight;      
    }// end iDim
  } // end jDim 

  assert(area);
  return area;
}