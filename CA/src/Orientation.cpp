#include "Orientation.h"
/*=======================================================================
   Class Definition
     Orientation - orientation of a single grain
  =======================================================================*/
/*-------------
  Constructor
  -------------*/
Orientation::Orientation()
  : id_(-1),
    angle_(0.)
{
  // Nothing for now
}

Orientation::Orientation(int id, double angle)
  : id_(id),
    angle_(angle)
{
  // Nothing for now
}

Orientation::Orientation(int id, int nDim, double* angle)
{
  id_ = id;
  if (nDim == 2)
    angle_ = angle[0];
  else
    for (int i = 0; i < 3; i++)
      EulerAngle_[i] = angle[i];
}

void Orientation::
get_Euler_angles(double EulerAngle[3])
{
  for (int i = 0; i < 3; i++) {
    EulerAngle[i] = EulerAngle_[i];
  }
}