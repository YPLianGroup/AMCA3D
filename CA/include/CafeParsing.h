#ifndef CAFEPARSING_H
#define CAFEPARSING_H

// TPL & STL
#include <yaml-cpp/yaml.h>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "CafeEnv.h"
#include "ParsingHelper.h"

// our data types
struct Coordinates
{
  double x_, y_, z_;
  Coordinates()
    :x_(0.0), y_(0.0), z_(0.0)
  {}
};

struct LateralSizes
{
  double xL_, yL_, zL_;
  LateralSizes()
    :xL_(1.0), yL_(1.0),zL_(1.0)
  {}
};

struct CellSizes
{
  double dCellx_, dCelly_, dCellz_;
  CellSizes()
    :dCellx_(0.01),dCelly_(0.01),dCellz_(0.01)
  {}
};

struct NumberCells
{
  int numbX_, numbY_, numbZ_;
  NumberCells()
    :numbX_(100),numbY_(100),numbZ_(100)
  {}
};

struct SiteDensity
{
  double siteDensity_;
  SiteDensity()
    :siteDensity_(5.5e1)
  {}
};

struct Mean
{
  double meanUndercooling_;
  Mean()
    :meanUndercooling_(6.0)
  {}
};

struct StandardDeviation
{
  double standardDeviation_;
  StandardDeviation()
    :standardDeviation_(0.1)
  {}
};

struct RandomType
{
  double rejectionRandom_;
  RandomType()
    :rejectionRandom_(false)
  {}
};

// Now the extraction operators for these types
void operator >> (const YAML::Node& node, Coordinates& X);
void operator >> (const YAML::Node& node, LateralSizes& L);
void operator >> (const YAML::Node& node, CellSizes& dCell);
void operator >> (const YAML::Node& node, NumberCells& numbCell);
void operator >> (const YAML::Node& node, SiteDensity& grainDensity);
void operator >> (const YAML::Node& node, Mean& deltaTmax);
void operator >> (const YAML::Node& node, StandardDeviation& deltaTsigma);
void operator >> (const YAML::Node& node, RandomType& randomType);

void operator >> (const YAML::Node& node, double * variable);
void operator >> (const YAML::Node& node, int * variable);

/* ----------------------------------------------------------------------------
   Set @param result if the @param key is present in the @param node, else set 
    it to the given default value
   ----------------------------------------------------------------------------*/
template<typename T>
void get_if_present(const YAML::Node& node, const std::string& key, T& result, 
                                        const T& default_if_not_present = T())
{
  const YAML::Node *value = node.FindValue(key);
  if (value)
    *value >> result;
  else
    result = default_if_not_present;
}

/*------------------------------------------------------------------------------
   this version doesn't change @param result unless the @param key is present
    in the @param node
  ------------------------------------------------------------------------------*/
template<typename T>
void get_if_present_no_default(const YAML::Node& node, const std::string& key, T& result)
{
  const YAML::Node *value = node.FindValue(key);
  if (value)
    *value >> result;
}

/*------------------------------------------------------------------------------
   this version requires the @param key to be present
  ------------------------------------------------------------------------------*/
template<typename T>
void get_required(const YAML::Node& node, const std::string& key, T& result)
{
  const YAML::Node *value = node.FindValue(key);
  if (value)
    *value >> result;
  else
  {
    if (!CafeEnv::self().parallel_rank())
    {
      std::ostringstream err_msg;
      err_msg << "\n\n Error: parsing missing required key: " << key
              << " at " << ParsingHelper::line_info(node)
              << " for node= " << std::endl;
      ParsingHelper::emit(err_msg, node);
      std::cout << err_msg.str() << std::endl;
    }
    throw std::runtime_error("Error: parsing missing required key: " + key);
  }
}

/*----------------------------------------------------------------------------
    these can be used to check and ensure a type of yaml node is as expected
  ----------------------------------------------------------------------------*/
const YAML::Node *
expect_type(const YAML::Node& node, const std::string& key, YAML::NodeType::value type, bool optional=false);

const YAML::Node *
expect_null(const YAML::Node& node, const std::string& key, bool optional=false);

const YAML::Node *
expect_scalar(const YAML::Node& node, const std::string& key, bool optional=false);

const YAML::Node *
expect_sequence(const YAML::Node& node, const std::string& key, bool optional=false);

const YAML::Node *
expect_map(const YAML::Node& node, const std::string& key, bool optional=false);

#endif
