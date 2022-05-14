#include <mpi.h>
#include <boost/filesystem.hpp>
#include <yaml-cpp/yaml.h>
#include "CellularAutomataManager_3D.h"
#include "VtkManager.h"
#include "CafeParsing.h"
#include "Simulation.h"
#include "CafeEnv.h"
#include "CafeMacro.h"
/*------------------------
   ~ constructor
  ------------------------*/
VtkManager::VtkManager()
  : numOutputStep_(0),
    outputFreq_(0),
    pvtkFileNum_(0),
    timeStampPrecision_(8),
    timeStepWidth_(5),
    vtkWidth_(5),
    pvtkWidth_(5),
    pvdWidth_(0),
    outputBasedonFreq_(true),
    outputTimeThreshold_(1.0e8)
{
  // nothing for now  
}

/*------------------------
   ~ destructor 
  ------------------------*/
VtkManager::~VtkManager()
{
  //nothing for now
}

/*------------------------------------------
  ~ response_to_control_of_output_variables
  ------------------------------------------*/
int 
VtkManager::response_to_output_varialbes_controls(std::string variable)
{
  bool gotIt = false;
  int iD = 0;
  do
  {
    if (variable == CAFE::variableName[iD])
    {
      gotIt = true;
    }
    else
    {
      iD++;
    }
  } while (!gotIt && iD < CAFE::nbNames);

  if (!gotIt)
  {
    CafeEnv::self().caOutputP0() << variable << " is not supported for now by the code base" << std::endl;
    iD = -1;
  }
  return iD;
}

/*--------------------------------------------------------------------
   bool output_result(int numStep, double time); // do we want output result now
  --------------------------------------------------------------------*/
bool VtkManager::output_result(int numStep, double time)
{
  bool output = false;
  if (outputBasedonFreq_) {
    if (numStep % outputFreq_ == 0) {
      output = true;
    }
  }
  else {
    if (time >= outputTimeThreshold_) {
      output = true;
      outputTimeThreshold_ += outputTimeInterval_;
      outputTimeThreshold_ = std::max(int(time / outputTimeInterval_+1)* outputTimeInterval_, outputTimeThreshold_);
    }
  }

  return output;
}
/*--------------------------------------------------------------------
   Initialize pvd File
   pvd_initialization
  --------------------------------------------------------------------*/
void
VtkManager::pvd_initialization()
{
  // Open filestream
  std::string pvdFileName = pvdPath_ + "/" + pvdName_ + ".pvd";
  pvdOut_.open(pvdFileName.c_str());

  // Write Headers
  pvdOut_ << "<?xml version=\"1.0\"?>\n";
  pvdOut_ << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  pvdOut_ << CAFE::SecondEleWidth_ << "<Collection>\n";
}

/*----------------------------------------------------------------------
   Write Time Stamp for Current vtr or vtu File
   pvd_time_stamp
  ---------------------------------------------------------------------*/
void
VtkManager::pvd_time_stamp(int & timeStep, double & timestamp, int & counter)
{
  // Get padded strings for output
  std::string timeStepString;
  std::string pvtkFileNumString;
  int_to_padded_string(timeStepWidth_, timeStep, timeStepString);
  int_to_padded_string(pvtkWidth_, counter, pvtkFileNumString);

  // Set fixed form precision, stream width and stream precision
  pvdOut_.setf(std::ios::fixed);
  pvdOut_.precision(timeStampPrecision_);
  pvdOut_.fill('0');

  // Write time stamp to .pvd file
  pvdOut_ << CAFE::ThirdEleWidth_ << "<DataSet timestep=\"" << std::left <<
    timestamp << "\" group=\"\" part=\"0\" file=\"" <<
    "./" + vtkDir_ << "/timeStep_" << timeStepString <<
    "/" << pvdName_ << pvtkFileNumString << ".p"<<vtkExtension_<<"\"/>\n";
}


/*----------------------------------------------------------------------
  Close pvd File
  pvd_finalization
  ----------------------------------------------------------------------*/
void
VtkManager::pvd_finalization()
{
  // Close Headers
  pvdOut_ << CAFE::SecondEleWidth_ << "</Collection>\n";
  pvdOut_ << "</VTKFile>\n";

  // Close pvd file stream
  pvdOut_.close();
}

/*----------------------------------------------------------------------
    create_path
  ----------------------------------------------------------------------*/
void
VtkManager::create_path(const std::string &filename) const
{
  bool  error_found = false;
  std::ostringstream errmsg;

  if (procID_ == 0)
  {
    // Split filename into a path
    std::size_t found = filename.find_last_of("/");
    std::string dir_path = filename.substr(0, found);

    boost::filesystem::path dir(dir_path);

    if (!boost::filesystem::create_directory(dir))
    {
      //std::cout << "Fail to create the directory: " << dir_path;
      std::cout << "Fail to create the directory";
      std::cout << ", which may exist already! " << std::endl;
    }
  }
}

/*----------------------------------------------------------------------
  Method to create a directory to house all vtk files for a given timestep
  vtkgenerate_time_step_directory
  ----------------------------------------------------------------------*/
void
VtkManager::generate_time_step_directory()
{
  // Create directory to house vtr/vtu files for current timestep
  int fileType = 4;
  std::string timeStepPath;
  assemble_file_name(timeStepPath, fileType);
  create_path(timeStepPath);
}

/*----------------------------------------------------------------------
  Get 0-padded number string for filename creation
  getPaddedFileNumber
  ----------------------------------------------------------------------*/
void
VtkManager::int_to_padded_string(int & padWidth, int & number, std::string & paddedString)
{
  // If pvtr/pvtu or vtr/vtu file, figure out file number
  if (padWidth > 0)
  {  // Get number digits in file number
    std::string fileNumberString = std::to_string(number);
    size_t numDigits = fileNumberString.length();

    // Create string of zeros to satisfy desired width
    std::string zeroString = std::string(padWidth - numDigits, '0');

    // Concatenate and exit
    paddedString = zeroString + fileNumberString;
  }
  // If pvd file, return empty string 
  else
  {
    paddedString = "";
  }
}


/*----------------------------------------------------------------------
  Method to concatenate fileNames with appropriate numbering
  Overload - Unsigned int
  assemble_file_names
  ----------------------------------------------------------------------*/
void
VtkManager::assemble_file_name(std::string & fileName, int & fileType)
{
  // Declare string stream for int-> string conversion
  std::stringstream concatAid;
  std::string timeStepDir;
  std::string paddedFileNameNumber;
  std::string fileExtension;
  std::string numStepString;
  int fileNameNumber = 0;
  int padWidth;

  // Switch on fileType to set variables for concatenation
  switch (fileType)
  {
  case 1:  // vtr/vtu  file
    timeStepDir = "timeStep_";
    fileNameNumber = procID_;
    fileExtension = "." + vtkExtension_;  //".vtr" or ".vtu";
    padWidth = vtkWidth_;
    break;
  case 2:  // pvtr/pvtu file
    timeStepDir = "timeStep_";
    fileNameNumber = pvtkFileNum_;
    fileExtension = ".p" + vtkExtension_; // ".pvtr" or ".pvtu";
    padWidth = pvtkWidth_;
    break;
  case 3:  // pvd  file 
    timeStepDir = "";
    fileNameNumber = 0;
    fileExtension = ".pvd";
    padWidth = pvdWidth_;
    break;
  case 4:  // timeStep directory 
    timeStepDir = "timeStep_";
    fileNameNumber = 0;
    fileExtension = "";
    padWidth = 0;;
    break;
  }

  // Concatenate to create appropriate filename
  int_to_padded_string(timeStepWidth_, numOutputStep_, numStepString);
  int_to_padded_string(padWidth, fileNameNumber, paddedFileNameNumber);
  fileName = pvdPath_ + "/" + vtkDir_ + "/" + timeStepDir +
    numStepString + "/" + pvdName_ +
    paddedFileNameNumber + fileExtension;
}
/*----------------------------------------------------------------------
  Accessor Function: get filename base
  get_base_name
  ----------------------------------------------------------------------*/
std::string
VtkManager::get_base_name() const
{
  return fileBase_;
}

/*---------------------------------------------------------------------- 
  Accessor Function: get frame number
  get_file_number
  ----------------------------------------------------------------------*/
int
VtkManager::get_file_number() const
{
  return pvtkFileNum_;
}

/*----------------------------------------------------------------------
  Accessor Function: get output frequency
  get_output_frequency
  ----------------------------------------------------------------------*/
int
VtkManager::get_output_frequency() const
{
  return outputFreq_;
}

/*----------------------------------------------------------------------
  Mutator Function: set frame number
  set_file_number
  ----------------------------------------------------------------------*/
void
VtkManager::set_file_number(int val)
{
  pvtkFileNum_ = val;
}

/*----------------------------------------------------------------------
  Mutator Function: set output frequency
  set_output_frequency
  ----------------------------------------------------------------------*/
void
VtkManager::set_output_frequency(int freq)
{
  outputFreq_ = freq;
}
