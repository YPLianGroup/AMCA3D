#include <fstream>
#include "Timer.h"
/*=======================================================================
   Class Definition
     Timer - manage the timing of individual operation of interest
  =======================================================================*/

/*-----------------------------------------------------------------------
   constructor
  -----------------------------------------------------------------------*/
Timer::Timer()
{
  for (int i = 0; i < 7; i++) {
    timingArray_[i] = 0.0;
  }
  timingForEntireModeling_ = 0.0;

  CATimer_ = false;
  timingForCA_Simulation_ = 0.0;
  timingForNucleation_ = 0.0;
  timingForCapture_ = 0.0;
  timingForCrossProcCapture_ = 0.0;

  FEMTimer_ = false;
  timingForFEM_Simulation_ = 0.0;
  timingForNTI_ = 0.0;

  MCTimer_ = false;
  timingForMC_Simulation_=0.0;
}

/*-----------------------------------------------------------------------
    destructor
  -----------------------------------------------------------------------*/
Timer::~Timer()
{
  // Do nothing
}

/*--------------------------------------
   start/stop timing for entire modeling
  --------------------------------------*/
void 
Timer::start_timing_entire_modeling()
{
  timerForEntireModeling_ = MPI_Wtime();
}

void
Timer::stop_timing_entire_modeling()
{
  double currentTime = MPI_Wtime();
  timingForEntireModeling_ += currentTime - timerForEntireModeling_;
  timingArray_[0] = timingForEntireModeling_;
}

/*--------------------
    set_CA_timer
  --------------------*/
void
Timer::set_CA_timer()
{
  CATimer_ = true;
}
/*--------------------
    set_MC_timer
  --------------------*/
void
Timer::set_MC_timer()
{
  MCTimer_ = true;
}
/*--------------------
    set_FEM_timer
  --------------------*/
void
Timer::set_FEM_timer()
{
  FEMTimer_ = true;
}
/*-----------------------------------------------------------------------
     get current time
  -----------------------------------------------------------------------*/
double
Timer::get_current_time()
{
  double currentTime = MPI_Wtime();
  return currentTime;
}

/*-----------------------------------------------------------------------
   start/stop timing the simulation stage, namely time integration
  -----------------------------------------------------------------------*/
void 
Timer::start_timing_CA_simulation()
{
  timerForCA_Simulation_ = MPI_Wtime();
}

void
Timer::stop_timing_CA_simulation()
{
  double currentTime = MPI_Wtime();
  timingForCA_Simulation_ += currentTime - timerForCA_Simulation_;
  timingArray_[1] = timingForCA_Simulation_;
}

/*-----------------------------------------------------------------------
   start/stop timing the nucleation stage
  -----------------------------------------------------------------------*/
void
Timer::start_timing_nucleation()
{
  timerForNucleation_ = MPI_Wtime();
}

void
Timer::stop_timing_nucleation()
{
  double currentTime = MPI_Wtime();
  timingForNucleation_ += currentTime - timerForNucleation_;
  timingArray_[2] = timingForNucleation_;
}

/*-----------------------------------------------------------------------
   start/stop timing the capture stage except for the cross-process-capture
  -----------------------------------------------------------------------*/
void
Timer::start_timing_capture()
{
  timerForCapture_ = MPI_Wtime();
}

void
Timer::stop_timing_capture()
{
  double currentTime = MPI_Wtime();
  timingForCapture_ += currentTime - timerForCapture_;
  timingArray_[3] = timingForCapture_;
}

/*-----------------------------------------------------------------------
   start/stop timing the cross-process-capture stage
  -----------------------------------------------------------------------*/
void
Timer::start_timing_crossProcCapture()
{
  timerForCrossProcCapture_ = MPI_Wtime();
}

void
Timer::stop_timing_crossProcCapture()
{
  double currentTime = MPI_Wtime();
  timingForCrossProcCapture_ += currentTime - timerForCrossProcCapture_;
  timingArray_[4] = timingForCrossProcCapture_;
}

/*---------------------------------------------
   start/stop timing the temperature evolution
  ---------------------------------------------*/
void
Timer::start_timing_CA_temp()
{
    timerForCA_temp_ = MPI_Wtime();
}

void
Timer::stop_timing_CA_temp()
{
  double currentTime = MPI_Wtime();
  timingForCA_temp_ += currentTime - timerForCA_temp_;
  timingArray_[5] = timingForCA_temp_;
}

/*---------------------------------------------
   start/stop timing the cell active
  ---------------------------------------------*/
void
Timer::start_timing_CA_cell_active()
{
    timerForCA_active_ = MPI_Wtime();
}

void
Timer::stop_timing_CA_cell_active()
{
    double currentTime = MPI_Wtime();
    timingForCA_active_ += currentTime - timerForCA_active_;
    timingArray_[6] = timingForCA_active_;
}

/*---------------------------------------------
   start/stop_timing_FEM_nodal_time_integration
  ---------------------------------------------*/
void
Timer::start_timing_FEM_simulation()
{
    timerForFEM_Simulation_ = MPI_Wtime();
}

void
Timer::stop_timing_FEM_simulation()
{
    double currentTime = MPI_Wtime();
    timingForFEM_Simulation_ += currentTime - timerForFEM_Simulation_;
    timingArray_[7] = timingForFEM_Simulation_;
}

/*----------------------------------------------
   start/stop_timing_FEM_nodal_time_integration
  ----------------------------------------------*/
void
Timer::start_timing_FEM_nodal_time_integration()
{
  timerForNTI_ = MPI_Wtime();
}

void
Timer::stop_timing_FEM_nodal_time_integration()
{
  double currentTime = MPI_Wtime();
  timingForNTI_ += currentTime - timerForNTI_;
  timingArray_[8] = timingForNTI_;
}

/*----------------------------------------------
   Time for one MC step
  ----------------------------------------------*/
void
Timer::start_timing_MC_time_step()
{
    timerForMC_Simulation_ = MPI_Wtime();
}
void
Timer::stop_timing_MC_time_step()
{
    double currentTime = MPI_Wtime();
    timingForMC_Simulation_ += currentTime - timerForMC_Simulation_;
    timingArray_[9] = timingForMC_Simulation_;
}
/*-----------------------------------------------------------------------
   output the time cost information for each stage
  -----------------------------------------------------------------------*/
void 
Timer::output_timing_info(int npros, double *maxTimingArray, double *minTimingArray)
{
  std::ofstream timeFile;
  std::string fileName = "timeLog.dat";
  timeFile.open(fileName.c_str());
  
  timeFile << "---------------------------------------" << std::endl;
  timeFile << "Begin Timer Overview (unit s)" << std::endl;
  timeFile << "-----------------------" << std::endl;
  timeFile << "Number of Processors:  " << npros << std::endl;
  timeFile << std::endl;
  timeFile << "Total Time for Entire Modeling:" << std::endl;
  timeFile << "                               " << maxTimingArray[0] << std::endl;
  timeFile << std::endl;
  if (FEMTimer_) {
    timeFile << "-------------" << std::endl;
    timeFile << "FEM Modeling: " << std::endl;       
    timeFile << "Timing for Simulation: Maximum          Minimum 0\n";
    timeFile << "                    " << maxTimingArray[7] << "      " << minTimingArray[7] << std::endl;
    timeFile << "Timing for Nodal Time Integration: Maximum          Minimum 0\n";
    timeFile << "                    " << maxTimingArray[8] << "      " << minTimingArray[8] << std::endl;
    timeFile << std::endl;
  }

  if (CATimer_) {
    timeFile << "-------------" << std::endl;
    timeFile << "CA Modeling:" << std::endl;    
    timeFile << "Timing for Simulation: Maximum          Minimum 0\n";
    timeFile << "                    " << maxTimingArray[1] << "      " << minTimingArray[1] << std::endl;
    timeFile << std::endl;
    timeFile << "Timing for Nucleatoin: Maximum          Minimum 0\n";
    timeFile << "                    " << maxTimingArray[2] << "      " << minTimingArray[2] << std::endl;
    timeFile << std::endl;
    timeFile << "Timing for Capturing:  Maximum          Minimum 0\n";
    timeFile << "                    " << maxTimingArray[3] << "      " << minTimingArray[3] << std::endl;
    timeFile << std::endl;
    timeFile << "Timing for Cross Capturing: Maximum       Minimum 0\n";
    timeFile << "                    " << maxTimingArray[4] << "      " << minTimingArray[4] << std::endl;
    timeFile << std::endl;
    timeFile << "Timing for temperature evolution: Maximum       Minimum 0\n";
    timeFile << "                    " << maxTimingArray[5] << "      " << minTimingArray[5] << std::endl;
    timeFile << std::endl;
    timeFile << "Timing for cell active: Maximum       Minimum 0\n";
    timeFile << "                    " << maxTimingArray[6] << "      " << minTimingArray[6] << std::endl;
  }
  if (MCTimer_) {
      timeFile << "-------------" << std::endl;
      timeFile << "MC Modeling:" << std::endl;
      timeFile << "Timing for Simulation: Maximum          Minimum 0\n";
      timeFile << "                    " << maxTimingArray[9] << "      " << minTimingArray[9] << std::endl;
  }

  timeFile.close();

  std::cout<< "**************************************************************" << "\n";
  std::cout<< "Simulation Shall Complete                      " << "\n";
  std::cout << "**************************************************************" << "\n";
  std::cout<< "-----------------------" << std::endl;
  std::cout<< "Number of Processors:  " << npros << std::endl;
  std::cout<< std::endl;
  std::cout<< "Total Time for Entire Modeling:" << std::endl;
  std::cout<< "                               " << maxTimingArray[0] << std::endl;
  std::cout << "--------------------------------------------------------------" << std::endl;
  std::cout << std::endl;
}