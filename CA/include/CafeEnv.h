#ifndef CAFEENV_H
#define CAFEENV_H

#include <mpi.h>
#include <fstream>
#include <streambuf>

class CAEmptyStreamBuffer : public std::filebuf
{
public:
  int overflow(int c) { return c; }
};

class CafeEnv
{
public:
  CafeEnv();
  virtual ~CafeEnv();

  static CafeEnv &self();

  MPI_Comm parallelCommunicator_;
  int pSize_;
  int pRank_;
  std::ostream *caLogStream_;
  std::ostream *caParallelStream_;

  CAEmptyStreamBuffer caEmptyStreamBuffer_;
  std::filebuf caStreamBuffer_;

  std::ostream & caOutputP0();
  std::ostream & caOutput();

  MPI_Comm parallel_comm();
  int parallel_size();
  int parallel_rank();
  void set_log_file_stream(std::string caLogName);
  void close_log_file_stream();
};

#endif

