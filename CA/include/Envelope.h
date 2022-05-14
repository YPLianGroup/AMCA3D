#ifndef ENVELOPE_H
#define ENVELOPE_H

// Forward declarations
class Grain;
class Orientation;

class Envelope
{
public:
  // Constructor and destructor
  Envelope();
  virtual ~Envelope();

  virtual bool cell_is_captured(const Grain* grain, double* cellCenter)=0;

  virtual void capture_cell(double * cellCenter, double length, double h_, double * newGrainInfo)=0;

  virtual void set_coordinate_transformation_matrix(double* angle_)=0;

  double length;

  double truncationCoefficient_;
  virtual void set_truncation_coefficient( double halfDiagonalLength)=0;
};

#endif

