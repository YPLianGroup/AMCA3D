#include <iostream>
#include "CafeEnv.h"
/*=======================================================================
   Class Definition
     CaEnv - manage parallel output in CA
           (will include parallel in future)
  =======================================================================*/

/*-----------------------------------------------------------------------
   constructor
  -----------------------------------------------------------------------*/
CafeEnv::CafeEnv()
  : parallelCommunicator_(MPI_COMM_WORLD),
    pSize_(-1),
    pRank_(-1),
    caLogStream_(&std::cout),
    caParallelStream_(&std::cout)
{
    // initialize
  MPI_Comm_size(parallelCommunicator_, &pSize_);
  MPI_Comm_rank(parallelCommunicator_, &pRank_);
}

/*-----------------------------------------------------------------------
   self
  -----------------------------------------------------------------------*/
CafeEnv &
CafeEnv::self()
{
  static CafeEnv s;
  return s;
}

/*-----------------------------------------------------------------------
   caOutputP0
  -----------------------------------------------------------------------*/
std::ostream &
CafeEnv::caOutputP0()
{
  return *caLogStream_;
}

/*-----------------------------------------------------------------------
   caOutput
  -----------------------------------------------------------------------*/
std::ostream &
CafeEnv::caOutput()
{
  return *caParallelStream_;
}

/*-----------------------------------------------------------------------
   parallel_size
  -----------------------------------------------------------------------*/
int
CafeEnv::parallel_size()
{
  return pSize_;
}

/*-----------------------------------------------------------------------
   parallel_rank
  -----------------------------------------------------------------------*/
int
CafeEnv::parallel_rank()
{
  return pRank_;
}

/*-----------------------------------------------------------------------
   parallel_comm
  -----------------------------------------------------------------------*/
MPI_Comm
CafeEnv::parallel_comm()
{
  return parallelCommunicator_;
}

/*-----------------------------------------------------------------------
   set_log_file_stream
  -----------------------------------------------------------------------*/
void
CafeEnv::set_log_file_stream(std::string caLogName)
{
  if (pRank_ == 0)
  {
    caStreamBuffer_.open(caLogName.c_str(), std::ios::out);
    caLogStream_->rdbuf(&caStreamBuffer_);
  }
  else
    caLogStream_->rdbuf(&caEmptyStreamBuffer_);
}

/*-----------------------------------------------------------------------
   close_log_file_stream
  -----------------------------------------------------------------------*/
void
CafeEnv::close_log_file_stream()
{
  if (pRank_ == 0)
    caStreamBuffer_.close();
}

/*-----------------------------------------------------------------------
   destructor
  -----------------------------------------------------------------------*/
CafeEnv::~CafeEnv()
{
  //close_log_file_stream();
}