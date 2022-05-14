/*-------------------------------------------------------------------------*\   
You can redistribute this code and/or modify this code under the 
terms of the GNU General Public License (GPL) as published by the  
Free Software Foundation, either version 3 of the License, or (at 
your option) any later version. see <http://www.gnu.org/licenses/>.

Please see our website for relavant literature:
  Xiong f. et al, Journal of Materials Processing Technology, 2022 - https://linkinghub.elsevier.com/retrieve/pii/S0924013622000504
  Xiong f. et al, Materials & Design, 2021 - https://linkinghub.elsevier.com/retrieve/pii/S0264127520309461
  Lian y. et al, Materials & Design, 2019 - https://www.sciencedirect.com/science/article/pii/S0264127519301091
  Lian y. et al, Computational Mechanics, 2018 - https://link.springer.com/article/10.1007%2Fs00466-017-1535-8

For further information please contact us by email:
  Prof. Yanping Lian:   yanping.lian@bit.edu.cn

Description
    Main program
\*---------------------------------------------------------------------------*/
// mpi
#include <mpi.h>

// boost for input params
#include <boost/program_options.hpp>

// yaml for parsing
#include <yaml-cpp/yaml.h>

// STL
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>

// ameps
//#include "CafeParsing.h"
#include "Simulation.h"
#include "CafeEnv.h"

int main(int argc, char ** argv)
{
  // start up MPI
  if (MPI_SUCCESS != MPI_Init(&argc, &argv))
  {
    throw std::runtime_error("MPI_Init failed");
  }

  // CaEnv singleton
  CafeEnv &cafeEnv = CafeEnv::self();

  // command line iptions
  std::string inputFileName, logFileName;

  boost::program_options::options_description desc("Cafe supported options");
  desc.add_options()
    ("help,h", "Help message")
    ("version,v", "Code Version 1.0")
    ("input-deck,i", boost::program_options::value<std::string>(&inputFileName)->default_value("cafe.i"),
      "Analysis input file")
    ("log-file,o", boost::program_options::value<std::string>(&logFileName),
      "Analysis log file");

  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);

  boost::program_options::notify(vm);

  // deal with some default parameters
  if (vm.count("help"))
  {
    if (!cafeEnv.parallel_rank())
      std::cerr << desc << std::endl;
    return 0;
  }

  if (vm.count("version"))
  {
    if (!cafeEnv.parallel_rank())
      std::cerr << "Version: Cafe1.0" << std::endl;
    return 0;
  }

  std::ifstream fin(inputFileName.c_str());
  if (!fin.good())
  {
    if (!cafeEnv.parallel_rank())
      std::cerr << "Input file is not specified or does not exist: user specified or (default) name= "
      << inputFileName << std::endl;
    return 0;
  }
    
  // deal with logfile name; if none supplied, go with inputFileName.log
  if (vm.count("log-file"))
  {
    int dotPos = logFileName.rfind(".");
    if (-1 == dotPos)
    {
      // lacking extension
      logFileName = inputFileName + ".log";
    }
    else
    {
      // with extension; swap with .log
      logFileName = logFileName.substr(0, dotPos) + ".log";
    }
  }

  // deal with log file stream
  cafeEnv.set_log_file_stream(logFileName);

  // proceed with reading input file "document" from YAML
  YAML::Parser parser(fin);
  YAML::Node doc;

  try
  {
    parser.GetNextDocument(doc);
  }
  catch (YAML::ParserException &e)
  {
    std::cout << e.what() << std::endl;
  }

  Simulation sim(doc);
  sim.load(doc);
  sim.initialize();
  sim.run();

  // close log file stream 
  cafeEnv.close_log_file_stream();
    //will be moved to the destructor of CaEnv class

  MPI_Finalize();
    // will be moved to the destructor of CaEnv class

  // all done
  return 0;

}
