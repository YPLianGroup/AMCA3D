#include <yaml-cpp/yaml.h>
// Local includes
#include "CellularAutomataManager.h"
#include "FiniteElementManager.h"
#include "ProblemPhysics.h"
#include "CafeParsing.h"
#include <cmath>
#include <algorithm>

/*-----------------
   ProblemPhysics
  -----------------*/
ProblemPhysics::ProblemPhysics(CellularAutomataManager * caManager)
  : caManager_(caManager),
    meltingTemperature_(0.0),
    solidusTemperature_(-1.0e8)
{
  caManager_->set_problem_physics(this);
}

ProblemPhysics::ProblemPhysics()
{
}

ProblemPhysics::~ProblemPhysics()
{
  // Nothing for now
}

/*----------------------------------
  ProblemPhysics_RappazGandin1993
  ----------------------------------*/
ProblemPhysics_RappazGandin1993::
ProblemPhysics_RappazGandin1993(CellularAutomataManager * caManager)
  : ProblemPhysics(caManager),
    initialTemperature_(-0.2),
    Tdot_(-2.3),
    a1_(0.0),
    a2_(0.0),
    a3_(0.0),
    b1_(1.0),
    b2_(2.0),
    b3_(3.0)
{
    initialTemperature_ = 300;
}

ProblemPhysics_RappazGandin1993::
~ProblemPhysics_RappazGandin1993()
{
  
}

/*---------------------------
       load
  ---------------------------*/
void
ProblemPhysics_RappazGandin1993::load(const YAML::Node& node)
{
  // Make some things pretty in the log file
  CafeEnv::self().caOutputP0() << "\n" << "Pysical Model Review" << "\n";
  CafeEnv::self().caOutputP0() << "=============================" << "\n";


  const YAML::Node *RappazGandin = node.FindValue("problem_physics");
  if (RappazGandin)
  {
    get_required(*RappazGandin, "type", name_);
    if (name_ != "RappazGandin")
      throw std::runtime_error("parsing error: realm-problem_physics-RappazGandin");

    get_if_present(*RappazGandin, "initial_temperature",initialTemperature_,initialTemperature_);
    get_if_present(*RappazGandin, "liquidus_temperature", meltingTemperature_, meltingTemperature_);
    get_if_present(*RappazGandin, "melting_temperature", meltingTemperature_, meltingTemperature_);
    get_if_present(*RappazGandin, "solidus_temperature", solidusTemperature_, solidusTemperature_);
    minTemperatureToNucleate_ = solidusTemperature_;
    get_if_present(*RappazGandin, "minimum_temperature_for_nucleation", minTemperatureToNucleate_, minTemperatureToNucleate_);
    get_if_present(*RappazGandin, "t_dot", Tdot_, Tdot_);
    get_if_present(*RappazGandin, "a1", a1_, a1_);
    get_if_present(*RappazGandin, "a2", a2_, a2_);
    get_if_present(*RappazGandin, "a3", a3_, a3_);
    get_if_present(*RappazGandin, "b1", b1_, b1_);
    get_if_present(*RappazGandin, "b2", b2_, b2_);
    get_if_present(*RappazGandin, "b3", b3_, b3_);

    // yaml-cpp Receipt
    CafeEnv::self().caOutputP0() << "Problem_physics details gathered from input file and/or defaults:" << "\n";
    CafeEnv::self().caOutputP0() << "Model: " << name_ << "\n";
    CafeEnv::self().caOutputP0() << "Liquidus Temperature: " << meltingTemperature_ << "\n";
    CafeEnv::self().caOutputP0() << "Solidus Temperature: " << solidusTemperature_ << "\n";
    CafeEnv::self().caOutputP0() << "a1: " << a1_ << "\n";
    CafeEnv::self().caOutputP0() << "a2: " << a2_ << "\n";
    CafeEnv::self().caOutputP0() << "a3: " << a3_ << "\n";
    CafeEnv::self().caOutputP0() << "b1: " << b1_ << "\n";
    CafeEnv::self().caOutputP0() << "b2: " << b2_ << "\n";
    CafeEnv::self().caOutputP0() << "b3: " << b3_ << "\n";
    CafeEnv::self().caOutputP0() << "Tdot: " << Tdot_ << "\n";
    CafeEnv::self().caOutputP0() << "Initial Temperature: " << initialTemperature_ << "\n";
  }
  else
    throw std::runtime_error("parsing error: realm-problem_physics-RappazGandin");
}


/*---------------------------
   initialize
  ---------------------------*/
void
ProblemPhysics_RappazGandin1993::initialize()
{
  // nothing for now
}

/*---------------------------
   compute_growth_velocity
  ---------------------------*/
double
ProblemPhysics_RappazGandin1993::compute_growth_velocity(const Grain * grain,
                                                         double time)
{
  // Note that for this problem physics, temperature doesn't depend on
  // location, so don't bother getting the real cellID
  int iCell = 0;
  double T = compute_cell_temperature(iCell, time);
  double dT = meltingTemperature_ - T;
  double v = 0.0;
  if (dT > 0)
  {
    v = a1_*dT + a2_*dT*dT + a3_*dT*dT*dT;
    if (v < 0.0) {
      v = 0.0;
    }
  }
  return v;
}

/*-------------------------------------------------------------------
   compute_growth_velocity_with_temperature_input
    the temperature is from somewhere else instead of offered by
    the class itself
  -------------------------------------------------------------------*/
double
ProblemPhysics_RappazGandin1993::compute_growth_velocity_with_temperature_input
                                                            (double T)
{
  // Note that for this problem physics, temperature does depend on
  // location, so the grain itself should contain temperature info.
  double dT = meltingTemperature_ - T;
  double v = 0.0;
  if (dT > 0)
  {
    v = a1_*pow(dT,b1_) + a2_* pow(dT, b2_) + a3_* pow(dT, b3_);
    if (v < 0.0) {
      v = 0.0;
    }
  }
  return v;
}

/*--------------------------
  compute_cell_temperature
  --------------------------*/
double
ProblemPhysics_RappazGandin1993::compute_cell_temperature(int iCell,
                                                          double time)
{
    double T = initialTemperature_ + Tdot_ * time;
    return T;
}

/*========================================
    class Material_HeatConduction
  ========================================*/

/*---------------------
    constructor
  ---------------------*/
Material_HeatConduction::Material_HeatConduction()
{
  // nothing for now
}

/*---------------------
   destructor
  ---------------------*/
Material_HeatConduction::~Material_HeatConduction()
{
  // nothing for now
}

/*--------------
       load
  --------------*/
void 
Material_HeatConduction::load(const YAML::Node& node)
{
  // Make some things pretty in the log file
  CafeEnv::self().caOutputP0() << "\n" << "Heat Conduction Model Review" << "\n";
  CafeEnv::self().caOutputP0() << "=============================" << "\n";

  //const YAML::Node *matHeat = node.FindValue("heat_conduction");
  if (node.FindValue("heat_conduction"))
  {   
    get_required(node, "density", density_);   
    get_required(node, "solid_thermal_capacity", Cps_);
    Cpl_ = Cps_;
    get_if_present(node, "liquid_thermal_capacity", Cpl_, Cpl_);

    get_required(node, "solid_thermal_conductivity", solidK_);
    liquidK_ = solidK_;
    powderK_ = solidK_;
    get_if_present(node, "liquid_thermal_conductivity", liquidK_, liquidK_);
    get_if_present(node, "powder_thermal_conductivity", powderK_, powderK_);

    get_required(node, "liquidus_temperature", liquidusTemp_);            
    solidusTemp_ = liquidusTemp_;    
    get_if_present(node, "solidus_temperature", solidusTemp_, solidusTemp_);    
    get_required(node, "evaporation_temperature", evaporationTemp_);

    get_required(node, "specific_latent_heat", specLatentHeat_);
    // yaml-cpp Receipt
    CafeEnv::self().caOutputP0() << "Material details gathered from input file and/or defaults:" << "\n";
    CafeEnv::self().caOutputP0() << "Model: " << name_ << "\n";
    CafeEnv::self().caOutputP0() << "Solid thermal capacity:    " << Cps_ << "\n";
    CafeEnv::self().caOutputP0() << "Liquid thermal capacity:   " << Cpl_ << "\n";
    CafeEnv::self().caOutputP0() << "Solid thermal conductivity:" << solidK_ << "\n";
    CafeEnv::self().caOutputP0() << "Liquid thermal conductivity:" << liquidK_ << "\n";
    CafeEnv::self().caOutputP0() << "Specific latent Heat: " << specLatentHeat_ << "\n";
    CafeEnv::self().caOutputP0() << "Solidus tempearture:     " << solidusTemp_ << "\n";
    CafeEnv::self().caOutputP0() << "Liquidus tempearture:    " << liquidusTemp_ << "\n";
    CafeEnv::self().caOutputP0() << std::endl;
  }
  else
    throw std::runtime_error("parsing error: realm-problem_physics-RappazGandin");
}

/*---------------
    initialize  
  ---------------*/
void
Material_HeatConduction::initialize()
{
  // calculate the slope of K
  double deltaT = liquidusTemp_ - solidusTemp_;
  if (fabsf(deltaT) < 1.0e-8) {
    slopeK_ = 0;
  }
  else {
    slopeK_ = (liquidK_ - solidK_) / (liquidusTemp_ - solidusTemp_);
  }
}

/*-------------------------
   get_thermal_diffusivity
  -------------------------*/
double 
Material_HeatConduction::
get_thermal_diffusivity(double temperature)
{
  if (temperature < liquidusTemp_){
    return solidK_ / (density_*Cps_);
  }
  else {
    return liquidK_ / (density_*Cpl_);
  }
}

/*------------------------
   get_thermal_capacity
  ------------------------*/
double
Material_HeatConduction::
get_thermal_capacity(double temperature)
{
  if (temperature < solidusTemp_)
    return Cps_;
  else if (temperature < liquidusTemp_)
    return (Cps_ + specLatentHeat_ / (liquidusTemp_ - solidusTemp_));
  else
    return Cpl_;
}

double
Material_HeatConduction::
get_heat_conductivity(double temperature)
{
  if (temperature < solidusTemp_) {
    return solidK_;
  }
  else if (temperature < liquidusTemp_){    
    return solidK_ + (temperature - solidusTemp_)*slopeK_;
  }
  else {
    return liquidK_;
  }

}

double
Material_HeatConduction::
get_liquidus_temperature()
{
  return liquidusTemp_;
}

double
Material_HeatConduction::
get_solidus_temperature()
{
  return solidusTemp_;
}

double
Material_HeatConduction::
get_density(double temperature)
{
  /*
  if (temperature >= solidusTemp_){
     return ...
  }
  else {
     return ...
  }
  */
  return density_;
}

double
Material_HeatConduction::
get_evaporation_temperature()
{
  return evaporationTemp_;
}