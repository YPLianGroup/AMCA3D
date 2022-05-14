#ifndef PARSINGHELPER_H
#define PARSINGHELPER_H

#include <yaml-cpp/yaml.h>
#include <boost/lexical_cast.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "Simulation.h"

class ParsingHelper 
{
public:
  /*-----------------------------------------------------------------------------------
      emit
    -----------------------------------------------------------------------------------*/
  static void emit(YAML::Emitter& emout, const YAML::Node& node)
  {
    //recursive depth first
    YAML::NodeType::value type = node.Type();
    std::string out;
    switch (type)
    {
    case YAML::NodeType::Scalar:
      node >> out;
      emout << out;
      break;
      
    case YAML::NodeType::Sequence:
      emout << YAML::BeginSeq;
      for (unsigned int i = 0; i < node.size(); i++)
      {
        const YAML::Node & subnode = node[i];
        emit(emout, subnode);
      }
      emout << YAML::EndSeq;
      break;

    case YAML::NodeType::Map:
      emout << YAML::BeginMap;
      for (YAML::Iterator i = node.begin(); i != node.end(); i++)
      {
        const YAML::Node & key = i.first();
        const YAML::Node & value = i.second();
        key >> out;
        emout << YAML::Key << out;
        emout << YAML::Value;
        emit(emout, value);
      }
      emout << YAML::EndMap;
      break;

    case YAML::NodeType::Null:
      emout << " (empty) ";
      break;

    default:
      std::cerr << "Warning: emit: unknown/unsupported node type " << std::endl;
      break;
    }
  }

  /*--------------------------------------------------------------------------------
     emit - uses Emitter to print node to stream
    --------------------------------------------------------------------------------*/
  static void emit(std::ostream& sout, const YAML::Node& node)
  {
    YAML::Emitter out;
    emit(out, node);
    sout << out.c_str() << std::endl;
  }

  /*--------------------------------------------------------------------------------
     line_info
    --------------------------------------------------------------------------------*/
  static std::string line_info(const YAML::Node & node)
  {
    std::ostringstream sout;
    sout << "(pos, line, colum) = ("
         << node.GetMark().pos << ", "
         << node.GetMark().line << ","
         << node.GetMark().column << ")";
    return sout.str();
  }

  /*-------------------------------------------------------------------------------
      info
    -------------------------------------------------------------------------------*/
  static std::string info(const YAML::Node& node)
  {
    std::ostringstream sout;
    sout << "Node at " << line_info(node) << " => \n";
    emit(sout, node);
    return sout.str();
  }

  /*-------------------------------------------------------------------------------
    traverse - just prints nodes depth-first
    -------------------------------------------------------------------------------*/
  static void traverse(std::ostream& sout, const YAML::Node& node, unsigned int depth = 0)
  {
    //recursive depth first
    YAML::NodeType::value type = node.Type();
    std::string indent((size_t)depth*2, ' ');
    std::string out;
    switch (type)
    {
    case YAML::NodeType::Scalar:
      node >> out;
      sout << indent << "Scalar: " << out << std::endl;
      break;

    case YAML::NodeType::Sequence:
      sout << indent << "Sequence: " << std::endl;
      for (unsigned int i = 0; i < node.size(); i++)
      {
        const YAML::Node& subnode = node[i];
        sout << indent << "[" << i << "]" << std::endl;
        traverse(sout, subnode, depth + 1);
      }
      break;

    case YAML::NodeType::Map:
      sout << indent << "Map: " << std::endl;
      for (YAML::Iterator i = node.begin(); i != node.end(); i++)
      {
        const YAML::Node& key = i.first();
        const YAML::Node& value = i.second();

        key >> out;
        sout << indent << "Key: " << out << std::endl;
        sout << indent << "Value: " << std::endl;
        traverse(sout, value, depth + 1);
      }
      break;

    case YAML::NodeType::Null:
      sout << indent << "(empty)" << std::endl;
      break;

    default:
      std::cerr << "Warning: traverse: unknown/unsupported node type" << std::endl;
      break;
    }
  }

  /*-------------------------------------------------------------------------------
     find_nodes_given_key
     - returns a vector of nodes fthat match the given key (depth first traversal)
    -------------------------------------------------------------------------------*/
  static void find_nodes_given_key(const std::string& key, const YAML::Node& node, 
                                          std::vector<const YAML::Node *>& result)
  {
    // recursive depth first
    YAML::NodeType::value type = node.Type();
    if (type != YAML::NodeType::Scalar && type != YAML::NodeType::Null)
    {
      const YAML::Node *value = node.FindValue(key);
      if (value)
        result.push_back(&node);
    }

    switch (type)
    {
    case YAML::NodeType::Scalar:
      break;

    case YAML::NodeType::Sequence:
      for (unsigned int i = 0; i < node.size(); ++i)
      {
        const YAML::Node& subnode = node[i];
        find_nodes_given_key(key, subnode, result);
      }
      break;

    case YAML::NodeType::Map:
      for (YAML::Iterator i = node.begin(); i != node.end(); i++)
      {
        const YAML::Node& value = i.second();
        find_nodes_given_key(key, value, result);
      }
      break;

    case YAML::NodeType::Null:
      break;

    default:
      break;
    }
  }


};

#endif

