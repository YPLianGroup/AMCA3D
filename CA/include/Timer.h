#ifndef TIMER_H
#define TIMER_H

#include <mpi.h>
#include <iostream>
#include "CafeEnv.h"

class Timer
{
public:
  // Constructor and destructor
  Timer();
  virtual ~Timer();
  double get_current_time();
  double timerForEntireModeling_;
  double timingForEntireModeling_;
  double timingArray_[10];

  void start_timing_entire_modeling();
  void stop_timing_entire_modeling();

  // For CA
  bool CATimer_;
  void set_CA_timer();
  double timingForCA_Simulation_, timerForCA_Simulation_;
  double timingForNucleation_, timerForNucleation_;
  double timingForCapture_, timerForCapture_;
  double timingForCrossProcCapture_, timerForCrossProcCapture_;
  double timingForCA_temp_, timerForCA_temp_;
  double timingForCA_active_, timerForCA_active_;

  void start_timing_CA_simulation();
  void start_timing_nucleation();
  void start_timing_capture();
  void start_timing_crossProcCapture();
  void start_timing_CA_temp();
  void start_timing_CA_cell_active();

  void stop_timing_CA_simulation();
  void stop_timing_nucleation();
  void stop_timing_capture();
  void stop_timing_crossProcCapture();
  void stop_timing_CA_temp();
  void stop_timing_CA_cell_active();

  // For FEM
  bool FEMTimer_;
  void set_FEM_timer();

  double timingForFEM_Simulation_, timerForFEM_Simulation_;
  double timingForNTI_, timerForNTI_;
  
  void start_timing_FEM_simulation();
  void start_timing_FEM_nodal_time_integration();

  void stop_timing_FEM_simulation();
  void stop_timing_FEM_nodal_time_integration();

  void output_timing_info (int npros, double *maxTimingArray, double *minTimingArray);

  // For MC
  bool MCTimer_;
  void set_MC_timer();
  double timingForMC_Simulation_,timerForMC_Simulation_;
  //void start_timing_MC_simulation();
  //void stop_timing_MC_simulation();
  void start_timing_MC_time_step();
  void stop_timing_MC_time_step();
  
};

#endif // !TIMER_

