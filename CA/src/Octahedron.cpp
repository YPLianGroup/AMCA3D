#define M_PI 3.14159265358979323846 
// General includes
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <algorithm>

// Local includes
#include "Octahedron.h"
#include "Grain.h"
#include "Orientation.h"

#define TheodorusConst 1.732050807568877

/*=======================================================================
   Class Definition
     Octahedron - one type of 3D nucleation site
  =======================================================================*/

/*-----------------
   Octahedron
  -----------------*/
Octahedron::Octahedron()
{
  truncationCoefficient_ = sqrtf(3.0);
  // Set values
  /*normalVector[8][3] = { { 1,  1,  1 },{ -1,  1,  1 },{ -1, -1,  1 },{ 1, -1,  1 },
  { 1,  1, -1 },{ -1,  1, -1 },{ -1, -1, -1 },{ 1, -1, -1 } };*/
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 3; j++) {
      normalVector[i][j] = 1.0;
    }
  }

  normalVector[1][0] = -1.0;
  normalVector[2][0] = -1.0;
  normalVector[2][1] = -1.0;
  normalVector[3][1] = -1.0;

  for (int i = 4; i < 8; i++)
    for (int j = 0; j < 3; j++)
    {
      if (j == 2)
        normalVector[i][j] = -1 * normalVector[i - 4][j];
      else
        normalVector[i][j] = normalVector[i - 4][j];
    }

  /*cornerCoord[6][3] = { { length, 0,  0}, { 0,  length, 0},{-length, 0,  0},
    {0, -length, 0}, {0, 0, length},{0, 0, -length} };*/
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 3; j++)
      cornerCoord[i][j] = 0.0;
  cornerCoord[0][0] =  1.0;
  cornerCoord[1][1] =  1.0;
  cornerCoord[2][0] = -1.0;
  cornerCoord[3][1] = -1.0;
  cornerCoord[4][2] =  1.0;
  cornerCoord[5][2] = -1.0;

  /*int surfaceCorner[8][3] = { { 0, 1, 4 },{ 1, 2, 4 },{ 2, 3, 4 },
  { 3, 0, 4 },{ 0, 1, 5 },{ 1, 2, 5 },{ 2, 3, 5 },{ 3, 0, 5 } };*/
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      if (j == 2)
      {
        surfaceCorner[i][j] = 4;
        surfaceCorner[i + 4][j] = 5;
      }
      else
      {
        surfaceCorner[i][j] = i + j;
        surfaceCorner[i + 4][j] = surfaceCorner[i][j];
      }
         
    }
    if (i == 3)
    {
      surfaceCorner[i][1] = 0;
      surfaceCorner[i + 4][1] = 0;
    }
  }
}

/*-----------------
   destructor
  -----------------*/
Octahedron::~Octahedron()
{
  // Do nothing
}

/*-----------------
   cell_is_captured
  -----------------*/
bool
Octahedron::cell_is_captured(const Grain* grain, double* cellCenter)
{
  bool isCaptured = false;
  length = grain->length_;
  if (length == 0) return isCaptured;
  double cGrain[3];
  cGrain[0] = grain->xc_;
  cGrain[1] = grain->yc_;
  cGrain[2] = grain->zc_;
  
  // set coordinate transformation matrixs
  set_coordinate_transformation_matrix(grain->orientation_->EulerAngle_);

  // Compute "grain coordinates" of cell center
  double gcCellCenter[3], cx = 0, cy = 0, cz = 0;
  for (int i = 0; i < 3; i++)
  {
    gcCellCenter[i] = (cellCenter[i] - cGrain[i]);
  
    cx += gcCellCenter[i] * transMatrixGlobal2Local[i];
    cy += gcCellCenter[i] * transMatrixGlobal2Local[i + 3];
    cz += gcCellCenter[i] * transMatrixGlobal2Local[i + 6];
  }

  if (fabs(cx) + fabs(cy) + fabs(cz) <= 1.0 * length)
  {
    isCaptured = true;
  }
  return isCaptured;
}

/*-----------------
   capture_cell
  -----------------*/
void
Octahedron::capture_cell(double * gcCellCenter, double length, double h_, double *newGrainInfo)
{	
  if (length == -1) {
	  for (int j = 0; j < 3; j++)
		  newGrainInfo[j] = 0;
	  return;
  }
  // Compute "grain coordinates" of cell center
  double cx = 0, cy = 0, cz = 0;
  for (int i = 0; i < 3; i++)
  {
    cx += gcCellCenter[i] * transMatrixGlobal2Local[i];
    cy += gcCellCenter[i] * transMatrixGlobal2Local[i + 3];
    cz += gcCellCenter[i] * transMatrixGlobal2Local[i + 6];
  }
  // Cell center should be inside octahedron with half diagonal length
  double interception = fabs(cx) + fabs(cy) + fabs(cz);
  assert(interception <=  1.0 * length);

  // Define several variables for the following calculation
  double cornerCoord_[6][3];
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 3; j++)
      cornerCoord_[i][j] = cornerCoord[i][j]*length;
  
  /* Step 1: Determination of  the face <111> F indicated by the
              quadrant, which has engulfed the cell centre. */
  
  int quad = 0;
  if (cz >= 0)
  {
    if (cx >= 0 && cy >= 0) quad = 0;
    else if (cx >= 0 && cy <  0) quad = 3;
    else if (cx <  0 && cy >= 0) quad = 1;
    else if (cx <  0 && cy <  0) quad = 2;
  }
  else
  {
    if (cx >= 0 && cy >= 0) quad = 4;
    else if (cx >= 0 && cy <  0) quad = 7;
    else if (cx <  0 && cy >= 0) quad = 5;
    else if (cx <  0 && cy <  0) quad = 6;
  }
  
  /* Step 2: calculation of the coordinates of point A, the 
              projection of the cell centre on the face F. */
  double pA[3],pC[3];
  pC[0] = cx; pC[1] = cy; pC[2] = cz;
  double perpendicularDis = length - interception; /* take account 1/TheodorusConst
                                                       later on */
  for (int i = 0; i < 3; i++)
  {
    pA[i] = pC[i] + perpendicularDis * normalVector[quad][i] / 3;
  }

  /* Step 3: Determination of corner S1 on this face F that is closest to 
              the poitn A*/
  double distance[3] = { 0.0, 0.0, 0.0 };
  for (int i = 0; i < 3; i++)
  {
    const int ID = surfaceCorner[quad][i];
    for (int j = 0; j < 3; j++)
    {
      distance[i] += (pA[j] - cornerCoord_[ID][j]) * ((pA[j] - cornerCoord_[ID][j]));
    }
    
  }
  int cornerID[3];
  double minimumDis = distance[0];
  cornerID[0] = surfaceCorner[quad][0];
  for (int i = 1; i < 3; i++)
  {
    if (minimumDis > distance[i])
    {
      minimumDis = distance[i];  // Since you are smaller, let me take the value.
      cornerID[i] = cornerID[0]; 
      cornerID[0] = surfaceCorner[quad][i];  // and your ID
      
    }
    else
      cornerID[i] = surfaceCorner[quad][i];
  }
  
  /* Step 4: Calculation of the coordinates of the two points I and J which are the 
              projections of point A on the two segments S1S2 and S1S3, respectively*/
  // S1 is the closest corner to the point A
  int s1 = cornerID[0], s2 = cornerID[1], s3 = cornerID[2];
  
  double s1a[3],s1s2[3],s1s3[3];
  double norm_s1s2 = 0.0, norm_s1s3 = 0.0;
  double s1a_dot_s1s2 = 0, s1a_dot_s1s3 = 0.0;
  for (int i = 0; i < 3; i++)
  {
    s1a[i] = pA[i] - cornerCoord_[s1][i];
    
    s1s2[i] = cornerCoord_[s2][i] - cornerCoord_[s1][i];
    norm_s1s2 += s1s2[i] * s1s2[i];  // norm of the s1s2 itself
    s1a_dot_s1s2 += s1a[i] * s1s2[i]; // dot product between s1a and s1s2

    s1s3[i] = cornerCoord_[s3][i] - cornerCoord_[s1][i];
    norm_s1s3 += s1s3[i] * s1s3[i];  // norm of the s1s3 itself 
    s1a_dot_s1s3 += s1a[i] * s1s3[i]; // dot probudct between s1a and s1s3
  }

  norm_s1s2 = sqrtf(norm_s1s2);
  norm_s1s3 = sqrtf(norm_s1s3);
  // the length of segments IS1,JS1
  double Is1[2], Js1[2];
  Is1[0] = s1a_dot_s1s2 / norm_s1s2;
  Is1[1] = norm_s1s2 - Is1[0];
  Js1[0] = s1a_dot_s1s3 / norm_s1s3;
  Js1[1] = norm_s1s3 - Js1[0];

  /* Step 5: Calculation of the initial size of the octahedron*/
  double truncatedValue = truncationCoefficient_ * h_;
  if (Is1[0] < truncatedValue && Is1[1] < truncatedValue)
  {
    Is1[0] = norm_s1s2 / 2;
    Is1[1] = norm_s1s3 / 2;
  }

  // Define lenght L1 and L2
  double L1 = 0.5*(std::min(Is1[0], truncatedValue) + std::min(Is1[1], truncatedValue));
  double L2 = 0.5*(std::min(Js1[0], truncatedValue) + std::min(Js1[1], truncatedValue));

  // Compute new grain length
  double newGrainLength = sqrtf(2)*(std::max(L1, L2));

  /* Step 6: calculation of the position of the growth center, and return new grain 
             information, newGrainInfo, where first 3 elements are position and the 
         last one is the new grain length*/

  double magnitude = (length - newGrainLength);
  // transform the local coordinates to global coordinates
  for (int i = 0; i < 3; i++)
  {
    const double component = magnitude * cornerCoord_[s1][i] / (length*TheodorusConst);
    for (int j = 0; j < 3; j++)
      newGrainInfo[j] += component * transMatrixLocal2Global[i + j*3];
  }
  newGrainInfo[3] = newGrainLength;
}

/*---------------------------------------------
    determine coordinate transformation matrix
    set_coordinate_transformation_matrix
  ---------------------------------------------*/
void
Octahedron::set_coordinate_transformation_matrix(double* EulerAngle_)
{
  double alpha = EulerAngle_[0];
  double beta  = EulerAngle_[1];
  double gamma = EulerAngle_[2];
  double C_alpha = cos(alpha);
  double S_alpha = sin(alpha);
  double C_beta = cos(beta);
  double S_beta = sin(beta);
  double C_gamma = cos(gamma);
  double S_gamma = sin(gamma);

  // assemble transformation matrix
  transMatrixLocal2Global[0] =  C_alpha * C_gamma - S_alpha * C_beta * S_gamma;
  transMatrixLocal2Global[1] = -C_alpha * S_gamma - S_alpha * C_beta * C_gamma;
  transMatrixLocal2Global[2] =  S_alpha * S_beta;
  transMatrixLocal2Global[3] =  S_alpha * C_gamma + C_alpha * C_beta * S_gamma;
  transMatrixLocal2Global[4] = -S_alpha * S_gamma + C_alpha * C_beta * C_gamma;
  transMatrixLocal2Global[5] = -C_alpha * S_beta;
  transMatrixLocal2Global[6] =  S_beta * S_gamma;
  transMatrixLocal2Global[7] =  S_beta * C_gamma;
  transMatrixLocal2Global[8] =  C_beta;

  int count = 0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j<3; j++)
     transMatrixGlobal2Local[count++] = transMatrixLocal2Global[i + j*3];
  
}

/*-------------------------------
  set_truncation_coefficient
  -------------------------------*/
void
Octahedron::set_truncation_coefficient(double halfDiagonalLength)
{
  // in terms of the cell spacing
  truncationCoefficient_ = sqrtf(2)*halfDiagonalLength*0.5;
}