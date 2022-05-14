/*=======================================================================
                    Node
   Definition of Node Class.
  =======================================================================*/
#include "Node.h"

Node::Node()
  :mass_(0.0),
   fixTheta_(false),
   birth_(false)
{
  // nothing for now
  for (int i = 0; i < 3; i++) {
    velocity_[i] = 0.0;
  }
}

Node::Node(int nTimes)
{
  birth_ = false;
  fixTheta_ = false;
  theta_.resize(nTimes);
  grainAround_ = false;

  for (int i = 0; i < 3; i++) {
    velocity_[i] = 0.0;
  }
}

Node::~Node()
{
  // nothing for now
}

/*---------------
    =
 
Node
Node::operator=(Node &source)
{
  int nTimes = 0;
  if (source.theta_)
  {
    nTimes = 2;
  }
  Node copyOne(nTimes);
  copyOne.mass_ = source.mass_;
  copyOne.fixTheta_ = source.fixTheta_;
  for (int iCom = 0; iCom < 3; iCom++)
  {
    copyOne.coordinates_[iCom] = source.coordinates_[iCom];
  }
  if (nTimes > 0)
  {
    copyOne.theta_[0] = source.theta_[0];
    copyOne.theta_[1] = source.theta_[1];
  }

  return copyOne;
}
---------------*/