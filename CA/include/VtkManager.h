#ifndef VtkManager_H
#define VtkManger_H
#include <vector>
#include <string>
#include <fstream>

class VtkManager
{
public:
  // Constructor
  VtkManager();
  // Destructor
  virtual ~VtkManager();

  /*------------------------
     high-level data member
    ------------------------*/
  std::string vtkExtension_;   // could be vtr,vtu, ...

  /*-------------------------
     Parallel info
    -------------------------*/
  int procID_;
  int numProcs_;

  /*----------------------------
      Output variables controls 
    ----------------------------*/
  int numOutputStep_;       // the number of output steps
  int outputFreq_;          // the step interval for output
  double outputTimeThreshold_;  // Time to output
  double outputTimeInterval_;   // time interval for output
  bool outputBasedonFreq_;  // flag for output interval

  std::vector<int> outputVar_;      // list of output variables
  virtual int response_to_output_varialbes_controls(std::string variable);
  std::vector<std::pair<int, double>> timeVec_; // time stamp
  bool output_result(int numStep, double time); // do we want output result now
  /*-------------------
      Path Variables
    -------------------*/
  std::string fileBase_;
  std::string vtkDir_;
  std::string pvdPath_;
  std::string pvdName_;

  /*-------------------
     Output widths
    -------------------*/
  int timeStampPrecision_;
  int vtkWidth_;
  int pvtkWidth_;
  int pvdWidth_;
  int timeStepWidth_;

  /*----------------------------
      Pvd Methods
    ----------------------------*/
  int pvtkFileNum_;
  std::ofstream pvdOut_;

  void pvd_initialization();
  // pvd_initialization - open filestream to .pvd file and write file headers
  void pvd_time_stamp(int & timeStep, double & timestamp, int & counter);
  // pvd_time_stamp - write a timestap in .pvd file
  void pvd_finalization();
  
  /*--------------------------------------------
     create directory and prepare the file name
    --------------------------------------------*/
  virtual void create_path(const std::string &filename) const;
    // create_path - Create directory to house target file
  virtual void generate_time_step_directory();
    // generate_time_step_directory - create directory for a given timestep to house vtr/pvtr files.
  virtual void int_to_padded_string(int &padWidth, int &number, std::string &paddedString);
    // in_to_padded_string - convert an integer to a string with leading zeros
  virtual void assemble_file_name(std::string &fileName, int &fileType);
  // assemble_file_name - concatenate paths, a filename and integers itno absolute filename

  /*-------------------
     Getters/Setters
    -------------------*/
  std::string get_base_name() const;
    // get_base_name - get base name of output files
  int get_file_number() const;
    // get_file_number - get frame number
  int get_output_frequency() const;
    // get_output_frequency - get output frequency
  void set_file_number(int val);
    // set_file_number - set frame number
  void set_output_frequency(int freq);
    // set_output_frequency - set output frequency
};

#endif
