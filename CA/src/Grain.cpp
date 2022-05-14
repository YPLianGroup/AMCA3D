// Local includes
#include "CellularAutomataManager.h"
#include "Grain.h"
#include "Orientation.h"

/*=======================================================================
   Class Definition
     Grain - Information to descripe a grain
  =======================================================================*/
/*--------------
   Constructor
  --------------*/
Grain::Grain(CellularAutomataManager * caManager, bool maintainActiveGrainList)
  : active_(true)
{
  if(maintainActiveGrainList)
    caManager->activeGrainList_.push_back(this);
}
/*--------------
   Destructor
  --------------*/
Grain::~Grain()
{
  // Nothing for now
}
/*--------------
   reset_grain_envelope:
   reset the envelope length to 0, and center to cell center
  --------------*/
void Grain::reset_grain_envelope(const double coor[3]){
    length_ = 0;
    xc_ = coor[0];
    yc_ = coor[1];
    zc_ = coor[2];
}