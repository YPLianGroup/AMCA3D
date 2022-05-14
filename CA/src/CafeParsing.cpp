#include "CafeParsing.h"
#include "CafeEnv.h"
#include "Simulation.h"

// now the extraction operators for these types
void
operator >> (const YAML::Node& node, Coordinates& X)
{
  node[0] >> X.x_;
  node[1] >> X.y_;
  if (node.size() > 2)
    node[2] >> X.z_;
}

void
operator >> (const YAML::Node& node, LateralSizes& L)
{
  node[0] >> L.xL_;
  node[1] >> L.yL_;
  if (node.size() > 2)
    node[2] >> L.zL_;
}

void
operator >> (const YAML::Node& node, CellSizes& dCell)
{
  node[0] >> dCell.dCellx_;
  node[1] >> dCell.dCelly_;
  if (node.size() > 2)
    node[2] >> dCell.dCellz_;
}

void
operator >> (const YAML::Node& node, NumberCells& numbCell)
{
  node[0] >> numbCell.numbX_;
  node[1] >> numbCell.numbY_;
  if (node.size() > 2)
    node[2] >> numbCell.numbZ_;
}

void
operator >> (const YAML::Node& node, SiteDensity& grainDensity)
{
  node >> grainDensity.siteDensity_;
}

void
operator >> (const YAML::Node& node, Mean& deltaTmax)
{
  node >> deltaTmax.meanUndercooling_;
}

void
operator >> (const YAML::Node& node, StandardDeviation& deltaTsigma)
{
  node >> deltaTsigma.standardDeviation_;
}

void
operator >> (const YAML::Node& node, RandomType& randomType)
{
  node >> randomType.rejectionRandom_;
}

void
operator >> (const YAML::Node& node, double * variable)
{
  node[0] >> variable[0];
  node[1] >> variable[1];
  if (node.size() > 2)
    node[2] >> variable[2];
}

void
operator >> (const YAML::Node& node, int * variable)
{
  node[0] >> variable[0];
  node[1] >> variable[1];
  if (node.size() > 2)
    node[2] >> variable[2];
}

/*================================================================================*/

const YAML::Node *
expect_type(const YAML::Node& node, const std::string& key, YAML::NodeType::value type, bool optional)
{
  static std::string types[] = {"Null", "Scalar", "Sequence", "Map"};
  const YAML::Node *value = node.FindValue(key);
  std::ostringstream err_msg;
  if (!optional && !value)
  {
    if (!CafeEnv::self().parallel_rank())
    {
      err_msg << "Error: parsing expected required value " << key << " but it was not found at"
        << ParsingHelper::line_info(node)
        << " for Node= " << std::endl;
      ParsingHelper::emit(err_msg, node);
      std::cout << err_msg.str() << std::endl;
    }
    throw std::runtime_error("Error:: parsing");
  }

  if (value && (value->Type() != type))
  {
    if (!CafeEnv::self().parallel_rank())
    {
      err_msg << "Error: parsing expected type " << types[type] << " got type ="
        << types[value->Type()]
        << "  for key= " << key
        << " at " << ParsingHelper::line_info(node)
        << " node= " << std::endl;
      ParsingHelper::emit(err_msg, node);
      err_msg << "Check indentation of input file.";
      std::cout << err_msg.str() << std::endl;
    }
    throw std::runtime_error("Error: parsing - Check indentation of input file.");
  }
  return value;
}

const YAML::Node *
expect_null(const YAML::Node& node, const std::string& key, bool optional)
{
  return expect_type(node, key, YAML::NodeType::Null, optional);
}

const YAML::Node *
expect_scalar(const YAML::Node& node, const std::string& key, bool optional)
{
  return expect_type(node, key, YAML::NodeType::Scalar, optional);
}

const YAML::Node *
expect_sequence(const YAML::Node& node, const std::string& key, bool optional)
{
  return expect_type(node, key, YAML::NodeType::Sequence, optional);
}

const YAML::Node *
expect_map(const YAML::Node& node, const std::string& key, bool optional)
{
  return expect_type(node, key, YAML::NodeType::Map, optional);
}