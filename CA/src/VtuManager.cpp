#include <mpi.h>
#include <boost/filesystem.hpp>
#include "VtuManager.h"
#include "FiniteElementManager.h"
#include "CafeParsing.h"
#include "Node.h"
#include "Element.h"
#include "CafeMacro.h"

/*=======================================================================
  Class Definition
  VtrManager - manage parallel vtr output
=======================================================================*/

/*--------------------------------------------------------------------
Constructor
VtrManager
--------------------------------------------------------------------*/
VtuManager::VtuManager(FiniteElementManager *feManager)
  :VtkManager(),
   feManager_(feManager)
{
  vtkExtension_ = "vtu";
  fileBase_ = "./Results/FEM";
  vtkDir_ = "finiteElementTimeStates";
  pvdName_ = "finiteElementFilesIno";  
}

/*--------------------------------------------------------------------
  Destructor
  ~VtuManager
  --------------------------------------------------------------------*/
VtuManager::~VtuManager()
{
  // nothing for now
}

/*--------------------------------------------------------------------
  load relevant data from input file
  load
  --------------------------------------------------------------------*/
void
VtuManager::load(const YAML::Node& node)
{
  // Make some things pretty in the log file
  CafeEnv::self().caOutputP0() << "\n" << "VtuManager Review" << "\n";
  CafeEnv::self().caOutputP0() << "=============================" << "\n";

  const YAML::Node *vtu_output = node.FindValue("output");
  if (vtu_output)
  {
    get_if_present(*vtu_output, "output_frequency", outputFreq_, outputFreq_);
    if (outputFreq_ < 1) {
      outputBasedonFreq_ = false;
    }
    get_if_present(*vtu_output, "output_data_base_name", fileBase_, fileBase_);
    outputTimeInterval_ = -1;
    get_if_present(*vtu_output, "output_time_interval", outputTimeInterval_, outputTimeInterval_);
    if (outputTimeInterval_ > 0) {
      outputTimeThreshold_ = 0.0;
    }

    // Default variables for output
    {
      int ID = response_to_output_varialbes_controls("temperature");
      outputVar_.push_back(ID);
    }
    // output variables
    const YAML::Node * varList = vtu_output->FindValue("output_variables");
    if (varList)
    {
      for (size_t iVar = 0; iVar < varList->size(); iVar++)
      {
        std::string name;
        (*varList)[iVar] >> name;
        int ID = response_to_output_varialbes_controls(name);
        if (ID >= 0 && ID != 0 && ID != 5)
        {
          outputVar_.push_back(ID);
        }
      }
    } // end if varList
    
  }
  else
    throw std::runtime_error("parser error: realm-output");

}

/*----------------------------------------------------------------------
  initialize
  ----------------------------------------------------------------------*/
void
VtuManager::initialize()
{
  MPI_Comm_rank(MPI_COMM_WORLD, &procID_);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs_);
  //Notation: since we are using a special decomposition strategy, let
  numProcs_ = 1;

  // FIXME: this may cause problem
  numPoints_ = feManager_->numNodes_;
  numCells_  = feManager_->numElements_;
  numCellsOutput_ = numCells_;

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
  // If we want elements data output, generate necessary directories
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
  Open filestream to .vtu file and write header lines
  vtr_initialization
----------------------------------------------------------------------*/
void
VtuManager::vtu_initialization()
{
  // Open file stream to vtr output
  std::string fileName;
  int fileType = 1;
  assemble_file_name(fileName, fileType);
  vtuOut_.open(fileName.c_str());

  // Write Header for .vtu File
  vtuOut_ << "<?xml version=\"1.0\"?>\n";
  vtuOut_ << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  vtuOut_ << "<UnstructuredGrid>\n";
  vtuOut_ << "<Piece NumberOfPoints=\"" << numPoints_
              << "\" NumberOfCells=\""  << numCellsOutput_ << "\">\n";

}


/*----------------------------------------------------------------------
   Write Point Coordinates
   vtu_write_coordinates
----------------------------------------------------------------------*/
void
VtuManager::vtu_write_coordinates()
{
  vtuOut_ << "<Points>\n";
  // Write Coordinate Headers  
  vtuOut_ << "<DataArray type=\"Float64\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\">\n";

  for (int i=0; i<numPoints_; i++)
  {
    Node * node = &feManager_->nodeVector_[i];
    double * tempCoords = node->coordinates_;
    vtuOut_ << tempCoords[0] << " ";
    vtuOut_ << tempCoords[1] << " ";
    vtuOut_ << tempCoords[2] << "\n";
  }

  // Close Coordinates Section
  vtuOut_ << "\n";
  vtuOut_ << "</DataArray>\n";
  vtuOut_ << "</Points>\n\n";
}

/*----------------------------------------------------------------------
   Write Point Scalar
    vtu_write_point_scalar_data
  ----------------------------------------------------------------------*/
void
VtuManager::vtu_write_point_data()
{
  // Write headers for point (nodal) data
  vtuOut_ << "<PointData Scalars=\"temperature\">\n";

  int numOutputVar = outputVar_.size();
  for (int iVar = 0; iVar < numOutputVar; iVar++)
  {
    switch (outputVar_[iVar])
    {
      case -1: // Error
        break;
      case 0:  // temperature
      {
        vtuOut_ << "<DataArray type=\"Float64\" Name=\"temperature\" format=\"ascii\">\n";
        for (int i = 0; i < numPoints_; i++)
        {
          double tempTheta = feManager_->thetaArray_[i];
          vtuOut_ << tempTheta << " ";
          if (i && i % 100 == 0)
            vtuOut_ << "\n";
        }
        vtuOut_ << "\n";
        vtuOut_ << "</DataArray>\n";

        break;
      }
      case 1: //temperature_rate
      {
        break;
      }
      case 10: // velocity
      {
        vtuOut_ << "<DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (int i = 0; i < numPoints_; i++)
        {
          Node * node = &feManager_->nodeVector_[i];
          vtuOut_ << node->velocity_[0] << " " << node->velocity_[1] << " " << node->velocity_[2] << " ";
          if (i && i % 100 == 0){
            vtuOut_ << "\n";
          }
        }
        vtuOut_ << "\n";
        vtuOut_ << "</DataArray>\n";
      }
    }
  } 
  // Close headers for point (nodal) data
  vtuOut_ << "</PointData>\n";
}

/*-------------------------------------------------
   Write Cell data
   vtu_write_cell_data_scalar
  -------------------------------------------------*/
void 
VtuManager::vtu_write_cell_data_scalar()
{
  // Write headers for scalar values
  vtuOut_ << "<CellData Scalars=\"element_birth\">\n";
  /* if there is no cell data do not write out blank DataArray*/
  int numOutputVar = outputVar_.size();
  for (int iVar = 0; iVar < numOutputVar; iVar++)
  {
    switch (outputVar_[iVar]) {    
    case 5:  // element_birth
    {
      vtuOut_ << "<DataArray type=\"Int8\" Name=\"element_birth\" format=\"ascii\">\n";
      std::vector<Element *> *elementVector = &feManager_->elementVector_;
      for (int i = 0; i < numCells_; i++)
      {
        int birth = (*elementVector)[i]->birth_;
        if (!birth) {
          continue;
        }
        vtuOut_ << birth << " ";
        if (i && i % 100 == 0)
          vtuOut_ << "\n";
      }
      vtuOut_ << "\n";
      vtuOut_ << "</DataArray>\n";

      break;
    }    
    }
  }
  vtuOut_ << "</CellData>\n";
}
/*----------------------------------------------------------------------
   Write connectivity, Offsets and Types
   vtu_set_connect 
  ----------------------------------------------------------------------*/
void
VtuManager::vtu_set_connect()
{
  // Write cell data headers
  vtuOut_ << "<Cells>\n";
  vtuOut_ << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

  // Connectivity
  for (unsigned int iCell = 0; iCell < numCells_; iCell++)
  {
    Element * ele = feManager_->elementVector_[iCell];
    if (!ele->birth_){
      continue;
    }

    for (int i=0; i<8; i++)
      vtuOut_ << ele->nID_[i] <<" ";
    vtuOut_ << "\n";
  }
  vtuOut_ << "</DataArray>\n\n";

  // Offsets
  vtuOut_ << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  for (unsigned int iCell = 0; iCell < numCellsOutput_; iCell++)
  {    
    vtuOut_ << (iCell+1)*8 << " ";
    if (iCell && iCell % 100 == 0)
      vtuOut_ << "\n";
  }
  vtuOut_ << "</DataArray>\n\n";

  // Types
  vtuOut_ << "<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n";
  for (unsigned int iCell = 0; iCell < numCellsOutput_; iCell++)
  {   
    vtuOut_ << 12 << " ";
    if (iCell && iCell % 100 == 0)
      vtuOut_ << "\n";
  }
  vtuOut_ << "</DataArray>\n\n";

  // Close Cells 
  vtuOut_ << "</Cells>\n";
}

/*----------------------------------------------------------------------
  Close Current vtu File
  vtu_finalization
----------------------------------------------------------------------*/
void
VtuManager::vtu_finalization()
{
  // Close Data Categories
  vtuOut_ << CAFE::ThirdEleWidth_ << "</Piece>\n";
  vtuOut_ << CAFE::SecondEleWidth_ << "</UnstructuredGrid>\n";
  vtuOut_ << "</VTKFile>\n";

  // Close file stream
  vtuOut_.close();
}

/*----------------------------------------------------------------------
    Create pvtu file
    pvtu_generator
----------------------------------------------------------------------*/
void
VtuManager::pvtu_generator()
{
  // Open file stream to vtu output
  std::string pvtuPath;
  int fileType = 2;
  assemble_file_name(pvtuPath, fileType);
  pvtuOut_.open(pvtuPath.c_str());

  // Write Header for .vtu File
  pvtuOut_ << "<?xml version=\"1.0\"?>\n";
  pvtuOut_ << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  pvtuOut_ << "<PUnstructuredGrid GhostLevel=\"0\">\n";
  pvtuOut_ << "<PPointData  Scalars=\"temperature\">\n";

  int numOutputVar = outputVar_.size();
  for (int iVar = 0; iVar < numOutputVar; iVar++)
  {
    switch (outputVar_[iVar])
    {
    case -1: // Error
      break;
    case 0:  // temperature
    {
      pvtuOut_ << "<DataArray type=\"Float64\" Name=\"temperature\" format=\"ascii\"/>\n";
      break;
    }
    case 1: //temperature_rate
    {
      break;
    }
    case 10: // velocity
    {
      pvtuOut_ << "<DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\"/>\n";
    }
    }
  }
  pvtuOut_ << "</PPointData>\n";
  pvtuOut_ << "<PPoints>\n";
  pvtuOut_ << "<PDataArray type=\"Float64\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\"/>\n";
  pvtuOut_ << "</PPoints>\n";
  pvtuOut_ << "<PCells>\n";
  pvtuOut_ << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"/>\n";
  pvtuOut_ << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\"/>\n";
  pvtuOut_ << "<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\"/>\n";
  pvtuOut_ << "</PCells>\n";

  // Write out list of vtu files
  std::string paddedProcString;
  for (int proc = 0; proc < numProcs_; ++proc)
  {
    int_to_padded_string(vtkWidth_, proc, paddedProcString);
    pvtuOut_ << "<Piece Source=\"" << pvdName_ << paddedProcString << ".vtu\" />\n";
  }

  // Close pvtu file
  pvtuOut_ << "</PUnstructuredGrid>\n";
  pvtuOut_ << "</VTKFile>\n";
  pvtuOut_.close();

  // Update .pvtu file count
  int newFrame = pvtkFileNum_ + 1;
  set_file_number(newFrame);
}

/*----------------------------------------------------------------------
  Coordinate VtuManager methods to output particle data
  finite element snapshot
----------------------------------------------------------------------*/
void
VtuManager::execute(int numStep, double time)
{
  // Do we want to output grain data at any point in simulation?
    // If so, does the current time step satisfy the set output frequency?
    if (output_result(numStep, time))
    {
      // get the number of the active elements
      numCellsOutput_ = 0;
      for (int iCell = 0; iCell < numCells_; iCell++) {
        Element * ele = feManager_->elementVector_[iCell];
        if (ele->birth_) {
          numCellsOutput_++;
        }
      }
      // Create path to current time step
      generate_time_step_directory();

      // Throw barrier so no files are written to a directory which doesn't yet exist
      if (numProcs_ > 1)
        MPI_Barrier(MPI_COMM_WORLD);

      // Push back time stamp to times vector
      std::pair<int, double> timePair(numOutputStep_, time);
      timeVec_.push_back(timePair);

      // Initialize .vtu file and write headers
      vtu_initialization();

      // Write out Point data to .vtu file
      vtu_write_point_data();

      // Write cell data .vtu file
      vtu_write_cell_data_scalar();

      // Write out Coordinates to .vtu file
      vtu_write_coordinates();

      // Write out cell connectivity 
      vtu_set_connect();

      // Close vtu headers and close filestream
      vtu_finalization();

      if (procID_ == 0)
      {
        // Create pvtu file
        pvtu_generator();

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
