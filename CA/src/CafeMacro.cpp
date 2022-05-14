#include "CafeMacro.h"
#include <algorithm>
#include <cmath>
/*-----------------------------------
   cross product of vector: 
   Vector_c = Vector_a x Vector_b
  -----------------------------------*/
void 
CAFE::vector_cross_product(double* a, double* b, double* c)
{
  c[0] = a[1] * b[2] - b[1] * a[2];
  c[1] = b[0] * a[2] - a[0] * b[2];
  c[2] = a[0] * b[1] - b[0] * a[1];
}


/*-------------------------------------------------
    calculate the inverse of a third order matrix
  -------------------------------------------------*/
void 
CAFE::matrix_inverse(double a[][3], double inverseA[][3], double &determinant)
{
  determinant = a[0][0] * a[1][1] * a[2][2] + a[0][1] * a[1][2] * a[2][0] +
    a[1][0] * a[2][1] * a[0][2] - a[0][2] * a[1][1] * a[2][0] -
    a[0][0] * a[1][2] * a[2][1] - a[0][1] * a[1][0] * a[2][2];

  inverseA[0][0] = a[1][1] * a[2][2] - a[1][2] * a[2][1];
  inverseA[0][1] = a[1][0] * a[2][2] - a[1][2] * a[2][0];
  inverseA[0][2] = a[1][0] * a[2][1] - a[1][1] * a[2][0];

  inverseA[1][0] = a[0][1] * a[2][2] - a[0][2] * a[2][1];
  inverseA[1][1] = a[0][0] * a[2][2] - a[0][2] * a[2][0];
  inverseA[1][2] = a[0][0] * a[2][1] - a[0][1] * a[2][0];

  inverseA[2][0] = a[0][1] * a[1][2] - a[0][2] * a[1][1];
  inverseA[2][1] = a[0][0] * a[1][2] - a[0][2] * a[1][0];
  inverseA[2][2] = a[0][0] * a[1][1] - a[0][1] * a[1][0];

  for (unsigned char i = 0; i < 3; i++){
    for (unsigned char j = 0; j < 3; j++) {
      inverseA[i][j] *= pow(-1.0, i + j) / determinant;
    }
  }

  double swap;
  swap = inverseA[0][1];
  inverseA[0][1] = inverseA[1][0];
  inverseA[1][0] = swap;

  swap = inverseA[0][2];
  inverseA[0][2] = inverseA[2][0];
  inverseA[2][0] = swap;

  swap = inverseA[1][2];
  inverseA[1][2] = inverseA[2][1];
  inverseA[2][1] = swap;
}

/*-----------------------------------------
    get the transposition of matrix 
  -----------------------------------------*/
void 
CAFE::matrix_transposition(double matrix[][3], double matrixT[][3])
{
  for (unsigned char i = 0; i < 3; i++)
    for (unsigned char j = 0; j < 3; j++)
      matrixT[i][j] = matrix[j][i];
}

/*-----------------------------------
   calculate vector norm
  -----------------------------------*/
double
CAFE::vector_2nd_norm(double *vector, unsigned int n)
{
  double sum = 0.0;
  for (unsigned int i = 0; i < n; i++)
    sum += vector[i] * vector[i];

  return std::sqrt(sum);
}

/*---------------------------
    inner product of vector
  ---------------------------*/
double 
CAFE::vector_inner_product(double *vector1, double *vector2, unsigned int n)
{
  double sum = 0.0;
  for (unsigned int i = 0; i < n; i++)
    sum += vector1[i] * vector2[i];

  return sum;
}

/*--------------------------
    vector normalization
  --------------------------*/
double 
CAFE::vector_normalization(double* v)
{
  double magnitude = vector_2nd_norm(v, 3);
  double magInverse = 1.0 / magnitude;
  v[0] = v[0] * magInverse;
  v[1] = v[1] * magInverse;
  v[2] = v[2] * magInverse;
  return magnitude;
}

/*--------------------------------
     coordinates_transformation
  --------------------------------*/
void
CAFE::coordinate_transformation
(double transformMatrix[3][3], double inputVector[3], double outputVector[3])
{
  for (int j = 0; j < 3; j++) {
    outputVector[j] = 0.0;
    for (int i = 0; i < 3; i++) {
      outputVector[j] += transformMatrix[j][i] * inputVector[i];
    } // end i      
  } // end j
}