#include <mpi.h>
#include <boost/filesystem.hpp>
#include "VtrBinManager.h"
#include "CellularAutomataManager_3D.h"
#include "CafeParsing.h"
#include "Orientation.h"
#include "CafeMacro.h"
/*=======================================================================
Class Definition
VtrManager - manage parallel vtr output
=======================================================================*/

/*--------------------------------------------------------------------
Constructor
VtrManager
--------------------------------------------------------------------*/
VtrBinManager::VtrBinManager(CellularAutomataManager_3D *caManager)
  :VtkManager(),
  caManager_(caManager)
{
  vtkExtension_ = "vtr";
  fileBase_ = "./Results/Grains";
  vtkDir_ = "grainsTimeStates";
  pvdName_ = "grainFilesIno";
}

/*--------------------------------------------------------------------
Destructor
~VtrManager
--------------------------------------------------------------------*/
VtrBinManager::~VtrBinManager()
{
  // Clear up
  if (x1IDArray)
    delete[] x1IDArray;
  if (x2IDArray)
    delete[] x2IDArray;
  if (y1IDArray)
    delete[] y1IDArray;
  if (y2IDArray)
    delete[] y2IDArray;
  if (z1IDArray)
    delete[] z1IDArray;
  if (z2IDArray)
    delete[] z2IDArray;
}

/*--------------------------------------------------------------------
load relevant data from input file
load
--------------------------------------------------------------------*/
void
VtrBinManager::load(const YAML::Node& node)
{
  // Make some things pretty in the log file
  CafeEnv::self().caOutputP0() << "\n" << "VtrManager Review" << "\n";
  CafeEnv::self().caOutputP0() << "=============================" << "\n";

  const YAML::Node *vtr_output = node.FindValue("output");
  if (vtr_output)
  {
    get_if_present(*vtr_output, "output_frequency", outputFreq_, outputFreq_);
    if (outputFreq_ < 1) {
      outputBasedonFreq_ = false;
    }
    get_if_present(*vtr_output, "output_data_base_name", fileBase_, fileBase_);
    get_if_present(*vtr_output, "output_time_interval", outputTimeInterval_, outputTimeInterval_);
    if (outputTimeInterval_ > 0) {
      outputTimeThreshold_ = 0.0;
    }

    // output variables
    const YAML::Node * varList = vtr_output->FindValue("output_variables");
    if (varList)
    {
      for (size_t iVar = 0; iVar < varList->size(); iVar++)
      {
        std::string name;
        (*varList)[iVar] >> name;
        int ID = response_to_output_varialbes_controls(name);
        if (ID >= 0)
        {
          outputVar_.push_back(ID);
        }
      }
    } // end if varList
    if (outputVar_.size() == 0)
    {
      int ID = response_to_output_varialbes_controls("grain_ID");
      outputVar_.push_back(ID);
    }
  }
  else
    throw std::runtime_error("parser error: realm-output");

}

/*----------------------------------------------------------------------
initialize
----------------------------------------------------------------------*/
void
VtrBinManager::initialize()
{
  // caManager should be full initialized.
  MPI_Comm_rank(caManager_->cartComm_, &procID_);
  MPI_Comm_size(caManager_->cartComm_, &numProcs_);

  h_ = caManager_->h_;

  // set the extent
  coordOffset[0] = caManager_->h_ * caManager_->xStart_ + caManager_->x0_;
  coordOffset[1] = caManager_->h_ * caManager_->yStart_ + caManager_->y0_;
  coordOffset[2] = caManager_->h_ * caManager_->zStart_ + caManager_->z0_;

  // set the mesh info
  nxLocGhost_ = caManager_->nxLocGhost_;
  nyLocGhost_ = caManager_->nyLocGhost_;
  nzLocGhost_ = caManager_->nzLocGhost_;
  nxyLocGhost_ = nxLocGhost_ * nyLocGhost_;

  x1ID = caManager_->xStart_;
  y1ID = caManager_->yStart_;
  z1ID = caManager_->zStart_;
  x2ID = x1ID + caManager_->nxLocal_;
  y2ID = y1ID + caManager_->nyLocal_;
  z2ID = z1ID + caManager_->nzLocal_;

  if (procID_ == 0)
  {
    x1IDArray = new int[numProcs_];
    x2IDArray = new int[numProcs_];
    y1IDArray = new int[numProcs_];
    y2IDArray = new int[numProcs_];
    z1IDArray = new int[numProcs_];
    z2IDArray = new int[numProcs_];
  }

  // Gather each process's extent
  MPI_Gather(&x1ID, 1, MPI_INTEGER,
    x1IDArray, 1, MPI_INTEGER, 0, caManager_->cartComm_);
  MPI_Gather(&x2ID, 1, MPI_INTEGER,
    x2IDArray, 1, MPI_INTEGER, 0, caManager_->cartComm_);
  MPI_Gather(&y1ID, 1, MPI_INTEGER,
    y1IDArray, 1, MPI_INTEGER, 0, caManager_->cartComm_);
  MPI_Gather(&y2ID, 1, MPI_INTEGER,
    y2IDArray, 1, MPI_INTEGER, 0, caManager_->cartComm_);
  MPI_Gather(&z1ID, 1, MPI_INTEGER,
    z1IDArray, 1, MPI_INTEGER, 0, caManager_->cartComm_);
  MPI_Gather(&z2ID, 1, MPI_INTEGER,
    z2IDArray, 1, MPI_INTEGER, 0, caManager_->cartComm_);

  // Generate necessary directories
  std::size_t found = fileBase_.find_last_of("/");
  std::size_t dotPos = fileBase_.find_last_of(".");
  std::string dirPath_ = fileBase_.substr(0, found);
  if (dotPos)
    pvdName_ = fileBase_.substr(found + 1, dotPos);
  else
    pvdName_ = fileBase_.substr(found + 1);

  boost::filesystem::path absolutePath = boost::filesystem::system_complete(dirPath_);
  pvdPath_ = absolutePath.string();

  if (outputFreq_ != 0)
  {
    boost::filesystem::path full_path_(boost::filesystem::current_path());
    CafeEnv::self().caOutputP0() << "Current working path is : " << full_path_ << std::endl;

    std::string pvdFolder = pvdPath_ + "/" + "Dummy";
    create_path(pvdFolder);
    std::string resultsPath = pvdPath_ + "/" + vtkDir_ + "/" + "Dummy";
    create_path(resultsPath);
  }
}
/*----------------------------------------------------------------------
Open filestream to .vtr file and write header lines
vtr_initialization
----------------------------------------------------------------------*/
void
VtrBinManager::vtr_initialization()
{
  // Open file stream to vtr output
  std::string fileName;
  int fileType = 1;
  assemble_file_name(fileName, fileType);
  vtrOut_.open(fileName.c_str());

  // Write Header for .vtu File
  vtrOut_ << "<?xml version=\"1.0\"?>\n";
  vtrOut_ << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  vtrOut_ << CAFE::SecondEleWidth_ << "<RectilinearGrid WholeExtent=\""
    << x1ID << " " << x2ID << " "
    << y1ID << " " << y2ID << " "
    << z1ID << " " << z2ID << "\">\n";

  vtrOut_ << CAFE::ThirdEleWidth_ << "<Piece Extent=\""
    << x1ID << " " << x2ID << " "
    << y1ID << " " << y2ID << " "
    << z1ID << " " << z2ID << "\">\n";
}

/*----------------------------------------------------------------------
Write Grain Scalar: Cell Data
vtr_write_cell_data_scalar
----------------------------------------------------------------------*/
void
VtrBinManager::vtr_write_cell_data()
{
  // Write headers for cell data values
  vtrOut_ << CAFE::FourthEleWidth_ << "<CellData Scalars=\"grain_id\">\n";
  int numOutputVar = outputVar_.size();
  for (int iVar = 0; iVar < numOutputVar; iVar++)
  {
    switch (outputVar_[iVar])
    {
    case -1: // Error
      break;
    case 0:  // temperature
    {
      double * temp = caManager_->cellTemperature_;
      if (temp) {
        vtrOut_ << CAFE::FifthEleWidth_ << "<DataArray type=\"Float64\" Name=\"cell_temperature\" format=\"ascii\">\n";
        for (int k = 1; k<nzLocGhost_ - 1; k++) {
          for (int j = 1; j < nyLocGhost_ - 1; j++) {
            vtrOut_ << CAFE::SixthEleWidth_;
            for (int i = 1; i < nxLocGhost_ - 1; i++) {
              int val = temp[k*nxyLocGhost_ + j*nxLocGhost_ + i];
              vtrOut_ << val << " ";
            }
            vtrOut_ << std::endl;
          }
        }

        // Close headers for scalar values
        vtrOut_ << "\n";
        vtrOut_ << CAFE::FifthEleWidth_ << "</DataArray>" << std::endl;
      }

      break;
    }

    case 1: // temperature_rate
    {
      break;
    }
    case 2: // grain_ID
    {
      int *u = caManager_->cellOrientationID_;
      vtrOut_ << CAFE::FifthEleWidth_ << "<DataArray type=\"Int64\" Name=\"grain_id\" format=\"ascii\">\n";
      for (int k = 1; k<nzLocGhost_ - 1; k++)
      {
        for (int j = 1; j < nyLocGhost_ - 1; j++)
        {
          vtrOut_ << CAFE::SixthEleWidth_;
          for (int i = 1; i < nxLocGhost_ - 1; i++)
          {
            int val = u[k*nxyLocGhost_ + j*nxLocGhost_ + i];
            if (val >= 0 && caManager_->mappedGrainID)
              val = caManager_->mappedGrainID[val];
            vtrOut_ << val << " ";
          }
          vtrOut_ << "\n";
        }
      }
      // Close headers for scalar values
      vtrOut_ << "\n";
      vtrOut_ << CAFE::FifthEleWidth_ << "</DataArray>" << std::endl;

      break;
    }
    case 3: // grain_orientation
    {
      int * cellOrientationID = caManager_->cellOrientationID_;
      vtrOut_ << CAFE::FifthEleWidth_ << "<DataArray type=\"Float64\" Name=\"grain_orientation\" NumberOfComponents=\"3\" format=\"ascii\">\n";
      for (int k = 1; k<nzLocGhost_ - 1; k++)
      {
        for (int j = 1; j < nyLocGhost_ - 1; j++)
        {
          vtrOut_ << CAFE::SixthEleWidth_;
          for (int i = 1; i < nxLocGhost_ - 1; i++)
          {
            int val = cellOrientationID[k*nxyLocGhost_ + j*nxLocGhost_ + i];
            if (val >= 0)
            {
              Orientation * orien = caManager_->orientationVector_[val];
              for (int iAngle = 0; iAngle < 3; iAngle++)
              {
                vtrOut_ << orien->EulerAngle_[iAngle] << " ";
              }
              /*double alpha = orien->EulerAngle_[0];
              double beta = orien->EulerAngle_[1];
              double gamma = orien->EulerAngle_[2];
              double C_alpha = std::cos(alpha);
              double S_alpha = std::sin(alpha);
              double C_beta = std::cos(beta);
              double S_beta = std::sin(beta);
              double C_gamma = std::cos(gamma);
              double S_gamma = std::sin(gamma);
              vtrOut_ << C_alpha * C_gamma - S_alpha * C_beta * S_gamma <<" ";
              vtrOut_ << S_alpha * C_gamma + C_alpha * C_beta * S_gamma <<" ";
              vtrOut_ << S_beta * S_gamma <<" ";*/

            }
            else
            {
              vtrOut_ << "0" << " " << "0" << " " << "0" << " ";
            }

          }
          vtrOut_ << "\n";
        }
      }
      // Close headers for scalar values
      vtrOut_ << "\n";
      vtrOut_ << CAFE::FifthEleWidth_ << "</DataArray>" << std::endl;

      break;
    }

    case 4: // grain_velocity
    {
      break;
    }

    case 6: // melting_flag
    {
      vtrOut_ << CAFE::FifthEleWidth_ << "<DataArray type=\"Int8\" Name=\"melted_flag\" format=\"ascii\">\n";
      for (int k = 1; k<nzLocGhost_ - 1; k++) {
        for (int j = 1; j < nyLocGhost_ - 1; j++) {
          vtrOut_ << CAFE::SixthEleWidth_;
          for (int i = 1; i < nxLocGhost_ - 1; i++) {
            int val = caManager_->hasBeenMelted_[k*nxyLocGhost_ + j*nxLocGhost_ + i];
            vtrOut_ << val << " ";
          }
          vtrOut_ << "\n";
        }
      }
      // Close headers for scalar values
      vtrOut_ << "\n";
      vtrOut_ << CAFE::FifthEleWidth_ << "</DataArray>" << std::endl;

      break;
    }

    }
  }

  // Write ending for cell data values
  vtrOut_ << CAFE::FourthEleWidth_ << "</CellData>\n";
}

/*----------------------------------------------------------------------
Write Point Coordinates
vtr_write_coordinates
----------------------------------------------------------------------*/
void
VtrBinManager::vtr_write_coordinates()
{
  // Write Coordinate Headers  
  vtrOut_ << CAFE::FourthEleWidth_ << "<Coordinates>\n";

  // X_COORDINATES
  vtrOut_ << CAFE::FifthEleWidth_ << "<DataArray type=\"Float64\" name=\"X_Coordinates\" format=\"ascii\">\n";
  vtrOut_ << CAFE::SixthEleWidth_;
  for (int i = 0; i < nxLocGhost_ - 1; i++)
  {
    vtrOut_ << coordOffset[0] + i*h_ << " ";
  }
  vtrOut_ << CAFE::FifthEleWidth_ << "</DataArray>\n";

  // Y_COORDINATES
  vtrOut_ << CAFE::FifthEleWidth_ << "<DataArray type=\"Float64\" name=\"Y_Coordinates\" format=\"ascii\">\n";
  vtrOut_ << CAFE::SixthEleWidth_;
  for (int j = 0; j < nyLocGhost_ - 1; j++)
  {
    vtrOut_ << coordOffset[1] + j*h_ << " ";
  }
  vtrOut_ << CAFE::FifthEleWidth_ << "</DataArray>\n";

  // Z_COORDINATES
  vtrOut_ << CAFE::FifthEleWidth_ << "<DataArray type=\"Float64\" name=\"Z_Coordinates\" format=\"ascii\">\n";
  vtrOut_ << CAFE::SixthEleWidth_;
  for (int k = 0; k < nzLocGhost_ - 1; k++)
  {
    vtrOut_ << coordOffset[2] + k*h_ << " ";
  }
  vtrOut_ << CAFE::FifthEleWidth_ << "</DataArray>\n";
  vtrOut_ << CAFE::FourthEleWidth_ << "</Coordinates>\n\n";
}

/*----------------------------------------------------------------------
Close Current vtr File
vtr_finalization
----------------------------------------------------------------------*/
void
VtrBinManager::vtr_finalization()
{
  // Close Data Categories
  vtrOut_ << CAFE::ThirdEleWidth_ << "</Piece>\n";
  vtrOut_ << CAFE::SecondEleWidth_ << "</RectilinearGrid>\n";
  vtrOut_ << "</VTKFile>\n";

  // Close file stream
  vtrOut_.close();
}

/*----------------------------------------------------------------------
Create pvtr file
pvtr_generator
----------------------------------------------------------------------*/
void
VtrBinManager::pvtr_generator()
{
  // Open file stream to vtu output
  std::string pvtrPath;
  int fileType = 2;
  assemble_file_name(pvtrPath, fileType);
  pvtrOut_.open(pvtrPath.c_str());

  // Write Header for .vtu File
  pvtrOut_ << "<?xml version=\"1.0\"?>\n";
  pvtrOut_ << "<VTKFile type=\"PRectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";


  pvtrOut_ << CAFE::SecondEleWidth_ << "<PRectilinearGrid  WholeExtent=\""
    << 0 << " " << caManager_->nx_ << " "
    << 0 << " " << caManager_->ny_ << " "
    << 0 << " " << caManager_->nz_ << "\" \n";
  pvtrOut_ << "                  " << "GhostLevel = \"0\">\n";

  pvtrOut_ << CAFE::ThirdEleWidth_ << "<PCellData Scalars=\"grain_id\">\n";

  int numOutputVar = outputVar_.size();
  for (int iVar = 0; iVar < numOutputVar; iVar++)
  {
    switch (outputVar_[iVar])
    {
    case -1: // Error
      break;
    case 0:  // temperature
    {
      pvtrOut_ << CAFE::FourthEleWidth_ << "<PDataArray type=\"Float64\" Name=\"cell_temperature\"  format=\"ascii\"/>\n";
      break;
    }
    case 2:  // grain _ID
    {
      pvtrOut_ << CAFE::FourthEleWidth_ << "<PDataArray type=\"Int64\" Name=\"grain_id\" format=\"ascii\"/>\n";
      break;
    }
    case 3: // grain_orientation
    {
      pvtrOut_ << CAFE::FourthEleWidth_ << "<PDataArray type=\"Float64\" Name=\"grain_orientation\" NumberOfComponents=\"3\" format=\"ascii\"/>\n";
      break;
    }
    case 4: //
    {
      break;
    }
    case 6: // melted_flag
    {
      pvtrOut_ << CAFE::FourthEleWidth_ << "<PDataArray type=\"Int8\" Name=\"melted_flag\" format=\"ascii\"/>\n";
      break;
    }
    }
  }
  pvtrOut_ << CAFE::ThirdEleWidth_ << "</PCellData>\n";

  pvtrOut_ << CAFE::ThirdEleWidth_ << "<PCoordinates>\n";
  pvtrOut_ << CAFE::FourthEleWidth_ << "<PDataArray type=\"Float64\" format=\"ascii\"/>\n";
  pvtrOut_ << CAFE::FourthEleWidth_ << "<PDataArray type=\"Float64\" format=\"ascii\"/>\n";
  pvtrOut_ << CAFE::FourthEleWidth_ << "<PDataArray type=\"Float64\" format=\"ascii\"/>\n";
  pvtrOut_ << CAFE::ThirdEleWidth_ << "</PCoordinates>\n";

  // Write out list of vtr files
  std::string paddedProcString;
  for (int proc = 0; proc < caManager_->numProcs_; ++proc)
  {
    int_to_padded_string(vtkWidth_, proc, paddedProcString);

    pvtrOut_ << CAFE::ThirdEleWidth_ << "<Piece Extent=\""
      << x1IDArray[proc] << " " << x2IDArray[proc] << " "
      << y1IDArray[proc] << " " << y2IDArray[proc] << " "
      << z1IDArray[proc] << " " << z2IDArray[proc] << "\"\n";
    pvtrOut_ << CAFE::ThirdEleWidth_ << "      " << "Source=\""
      << pvdName_ << paddedProcString << ".vtr\" />\n";
  }

  // Close pvtu file
  pvtrOut_ << CAFE::SecondEleWidth_ << "</PRectilinearGrid>\n";
  pvtrOut_ << "</VTKFile>\n";
  pvtrOut_.close();

  // Update .pvtu file count
  int newFrame = pvtkFileNum_ + 1;
  set_file_number(newFrame);
}


/*----------------------------------------------------------------------
Coordinate VtuManager methods to output particle data
particleSnapShot
----------------------------------------------------------------------*/
void
VtrBinManager::execute(int numStep, double time)
{
  // Do we want to output grain data at any point in simulation?
  if (outputFreq_ != 0)
  {
    // If so, does the current time step satisfy the set output frequency?
    if (numStep % outputFreq_ == 0)
    {
      // Create path to current time step
      generate_time_step_directory();

      // Throw barrier so no files are written to a directory which doesn't yet exist
      MPI_Barrier(MPI_COMM_WORLD);

      // Push back time stamp to times vector
      std::pair<int, double> timePair(numOutputStep_, time);
      timeVec_.push_back(timePair);

      // Initialize .vtu file and write headers
      vtr_initialization();

      // Write out Cell data to .vtr file
      vtr_write_cell_data();

      // Write out Coordinates to .vtr file
      vtr_write_coordinates();

      // Close vtu headers and close filestream
      vtr_finalization();

      if (caManager_->procID_ == 0)
      {
        // Create pvtu file
        pvtr_generator();

        // Initialize .pvd file and write headers
        pvd_initialization();

        int counter = 0;
        // Write Time stamp to .pvd file
        for (std::vector<std::pair<int, double> >::iterator iter = timeVec_.begin(); iter != timeVec_.end(); ++iter)
        {
          pvd_time_stamp(iter->first, iter->second, counter);
          counter++;
        }
        // Close pvd headers and close filestream
        pvd_finalization();
      }
      numOutputStep_++;
    }
  }
}