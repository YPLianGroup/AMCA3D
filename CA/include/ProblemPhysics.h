#ifndef PROBLEM_PHYSICS_H
#define PROBLEM_PHYSICS_H

// Forward declarations
class CellularAutomataManager;
class FiniteElementManager;
class Grain;

class ProblemPhysics
{
 public:
  double couplingTimeStepFactor_;  // TimeStepFactor coupling to FEM or CA, value of -1 ( < = 0 ): using CA time step
                                 // value > 0: MC time step = FE time step * couplingTimeStepFactor_

  // Constructor and destructor
  ProblemPhysics(CellularAutomataManager * caManager);
  ProblemPhysics();
  virtual ~ProblemPhysics();

  virtual void load(const YAML::Node& node) = 0;

  virtual void initialize() = 0;
  virtual double compute_growth_velocity(const Grain * grain, double time) = 0;

  virtual double compute_growth_velocity_with_temperature_input(double T) = 0;

  virtual double compute_cell_temperature(int iCell, double time) = 0;

  double get_melting_temperature() { return meltingTemperature_; };
  double get_solidus_temperature() { return solidusTemperature_; };

  double get_minimum_temperature_for_nucleation() 
                                   { return minTemperatureToNucleate_; };

 protected:

  CellularAutomataManager * caManager_;
  double meltingTemperature_;
  double solidusTemperature_;
  double minTemperatureToNucleate_;
  std::string name_;
};


class ProblemPhysics_RappazGandin1993 : public ProblemPhysics
{
 public:

  // Constructor and destructor
  ProblemPhysics_RappazGandin1993(CellularAutomataManager * caManager);
  ~ProblemPhysics_RappazGandin1993();

  virtual void load(const YAML::Node& node);

  void initialize();

  virtual double compute_growth_velocity(const Grain * grain, double time);

  virtual double compute_growth_velocity_with_temperature_input(double T);

  virtual double compute_cell_temperature(int iCell, double time);

 protected:

  double initialTemperature_;
  double Tdot_;
  double a1_,b1_;
  double a2_,b2_;
  double a3_,b3_;
  
};



class Material_HeatConduction
{
public:
  Material_HeatConduction();
  virtual ~Material_HeatConduction();

  std::string name_;
  //FiniteElementManager* femManager_;
  /*---------------------
     physical parameters
    ---------------------*/
private:
  double density_;  
  double Cps_;   // specific thermal capacity for solid
  double Cpl_;   // specific thermal capacity for liquid
  double liquidK_;     // heat conductivity;
  double solidK_;      // heat conductivity;
  double slopeK_;      // slope of K between liquidus temp and solidus temp
  double powderK_;     // heat conductivity;
  double alpha_; // thermal diffusivity = k/(rho*Cp)
  double liquidusTemp_;  // liquidus temperature
  double solidusTemp_;   // solid temperature  
  double specLatentHeat_; // specific latent heat
  double evaporationTemp_;
public:
  /*--------------------------------
      Methods to load the model
    --------------------------------*/
  void load(const YAML::Node & node);
  void initialize();
  /*--------------------------------
     Methods to get some variables
    --------------------------------*/
  double get_thermal_diffusivity(double temperature);
  double get_thermal_capacity(double temperature);  
  double get_heat_conductivity(double temperature);
  double get_liquidus_temperature();
  double get_solidus_temperature();
  double get_evaporation_temperature();
  double get_density(double temperature);

};
#endif
