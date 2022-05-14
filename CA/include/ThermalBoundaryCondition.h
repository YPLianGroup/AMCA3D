#ifndef ThermalBoundaryCondition_h
#define ThermalBoundaryCondition_h
#include <vector>
#include <list>
#include <string>
#include "Simulation.h"

class ThermalBoundaryCondition
{
public:
  ThermalBoundaryCondition();
  virtual ~ThermalBoundaryCondition();

  virtual void load(const YAML::Node& node)=0;
  virtual void initialize() = 0;

  std::string type_;
  std::string target_;
  double ambientTemperature_;
  virtual double execute(double temperature) = 0;  
  virtual void get_value(double *parameters) = 0;
};


class FixedBoundaryCondition : public ThermalBoundaryCondition
{
public:
  FixedBoundaryCondition();
  ~FixedBoundaryCondition();
  void load(const YAML::Node& node);
  void initialize();
  double execute(double temperature);
  void get_value(double *parameters);  
  double theta_;
};

class FluxBoundaryCondition : public ThermalBoundaryCondition
{
public:
  FluxBoundaryCondition();
  ~FluxBoundaryCondition();

  void load(const YAML::Node& node);
  void initialize();
  double execute(double temperature);  
  void get_value(double *parameters);

  double flux_;
};


struct ToolPath {
public:
  double time_;
  double coords[3];
  bool turnOn_;
};
class MovingFluxBoundaryCondition : public ThermalBoundaryCondition
{
public:
  MovingFluxBoundaryCondition();
  ~MovingFluxBoundaryCondition();

  void load(const YAML::Node& node);
  void initialize();
  double execute(double distanceSquare);
  double execute_moving_flux(double distanceSquare, double sinTheta);
  void update_moving_flux_info(double time);
  void get_value(double *parameters);

  /*--------------------------
      Gaussian function
    --------------------------*/
  double sourcePower_;    // laser/electron beam power
  double beamRadius_;     // laser/electron beam radius
  double distributionFactor_; // laser/electron beam energy distribution factor
  double layerThickness_;     // layer thickness which is determined experimetally
  double yitaP_;              // the fraction of the laser/electron beam power aborbed
                              //  by the powder in-flight
  double yitaL_;              // the fraction of the laser/electron beam power absorbed
                              //  by the growing layer
  double para1_;        // taking the value of 1/(beamRadius_*sqrt(2*pi))
  double para2_;        // taking the value of 1/(2*beamRadius_^2)
  double supportRadius2_;

  /*--------------------------
      tool path
    --------------------------*/  
  int point0_;   // first point
  int point1_;   // second point
  double center_[3];    // heat sorce center
  double direction_[3]; // heat source direction
  bool birth_;    // heat flux is now on
  int numPoints_; // the number of time spots to define tool path
  std::vector<ToolPath *> toolPathVector_;
  void load_tool_path(std::string fileName);
  void set_time_scale_and_offset(double timeScale, double timeOffset);

  bool considerReflection_; // impose the flux on the side face, which is parallel to the flux direction  
};

class ConvectionBoundaryCondition : public ThermalBoundaryCondition
{
public:
  ConvectionBoundaryCondition();
  ~ConvectionBoundaryCondition();

  void load(const YAML::Node& node);
  void initialize();
  double execute(double temperature);
  void get_value(double *parameters);

  double heatTransferCoefficient_;
};

class RadiationBoundaryCondition : public ThermalBoundaryCondition
{
public:
  RadiationBoundaryCondition();
  ~RadiationBoundaryCondition();

  void load(const YAML::Node& node);
  void initialize();
  double execute(double temperature);  
  void get_value(double *parameters);

  double surfaceEmissivity_;
  double Stefan_BoltzmannConstant_;
  double coefficient_;
};

#endif
