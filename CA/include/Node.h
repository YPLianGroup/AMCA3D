#ifndef NODE_H
#define NODE_H
/*=======================================================================
                         Node
    Define the attributes of node for finite element method.
 =======================================================================*/
#include <vector>

class Node
{
public:
  //! Constructor
  Node();
  Node(int nTimes);

  //! Destructor
  virtual ~Node();

  /*-------------------------
      Basic Attributes
    -------------------------*/
  // the nodal mass, which depends on the physical problem.
  double mass_; 
  // the nodal velocity
  double velocity_[3];
  // the current temperature  
  std::vector<double>  theta_; 
  // the coordinates of the point  
  double coordinates_[3];
  // boundary conditions about Theta  
  bool fixTheta_;
  // active state
  bool birth_;
  // elements connected to the node
  std::vector<int> nEle_;
  bool grainAround_;
  //int matID_; // to be deleted soon

 // Node operator=(Node &source);
    // over load operator = since a pointer membera data is included
};

#endif
