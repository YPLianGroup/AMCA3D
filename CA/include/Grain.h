#ifndef GRAIN_H
#define GRAIN_H

// Forward declarations
class CellularAutomataManager;
class Orientation;

class Grain
{
public:
  
  // Constructor/destructor
  Grain(CellularAutomataManager * caManager, bool maintainActiveGrainList = true);
  ~Grain();

  int id_;                    // unique grain id
  int cell_;                  // local id of owning cell
  double xc_;                 // x coordinate of grain center
  double yc_;                 // y coordinate of grain center
  double zc_;                 // z coordinate of grain center
  Orientation * orientation_; // grain orientation
  double length_;             // length of grain along diagonal  

  bool is_active() {return active_;}
  void inactivate() {active_ = false;}
  void activate() { active_ = true; }
  void grow(double velocity, double dt) {length_ += dt*velocity;}
  void reset_grain_envelope(const double coor[3]);

protected:
  bool active_;               // flag whether grain is actively growing

};

class BoundaryNucleal
{
	public:
		BoundaryNucleal(double Tc,int cellid) :Tmax_(-1),Tincubation_(0), DeltaTc_(Tc), cellid_(cellid) {}
		~BoundaryNucleal() {};
		double Tmax_;			// maxmium temperature
		double Tincubation_;	// incubation time
		double DeltaTc_;		// critical undercooling
		int cellid_;
};

#endif
