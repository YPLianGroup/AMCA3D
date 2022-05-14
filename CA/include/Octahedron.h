#ifndef OCTAHEDRON_H
#define OCTAHEDRON_H
// Base class
#include "Envelope.h"

// Forward declarations
class Grain;
class Orientation;

class Octahedron : public Envelope
{
public:
  // Constructor and destructor
  Octahedron();
  virtual ~Octahedron();

  virtual bool cell_is_captured(const Grain* grain, double* cellCenter);

  virtual void capture_cell(double * cellCenter, double length, double h_, double * newGrainInfo);
  virtual void set_coordinate_transformation_matrix( double* EulerAngle_);
  void set_truncation_coefficient(double halfDiagonalLength);

  // data members
  double transMatrixLocal2Global[9]; // coordinate transformation matrix, from local to global
  double transMatrixGlobal2Local[9]; // coordinate transformation matrix, from global to local

  double normalVector[8][3];
  double cornerCoord[6][3];
  int surfaceCorner[8][3];
};

#endif

