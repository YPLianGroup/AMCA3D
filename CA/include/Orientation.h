#ifndef ORIENTATION_H
#define ORIENTATION_H

class Orientation
{
public:
  // Constructor/destructor
  Orientation();
  Orientation(int id, double angle);
  Orientation(int id, int nDim, double* angle);
  ~Orientation() {};

  int id_;       // unique ID associated with this orientation
  double angle_; // orientation angle
  double EulerAngle_[3]; // orientation angles
  void get_Euler_angles(double EulerAngle[3]);
};

#endif
