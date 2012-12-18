#ifndef INIVISITOR_H
#define INIVISITOR_H

#include "params.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <iostream>

/// Simple handler for ini formatted files
/// Parameter object must have at most one level (section.x = val)
class INIHandler : public ParamFileHandler {
public:
  INIHandler(const std::string& _root) : ParamFileHandler(_root) {}
  void read(std::istream& file, ptree& pt) const
    { boost::property_tree::ini_parser::read_ini(file, pt); }
  void write(std::ostream& file, const ptree& pt) const
    { boost::property_tree::ini_parser::write_ini(file, pt); }
  static std::string getPath(const std::string& root, const std::string& name, const std::string& path = "")
    { return ParamFileHandler::getPath(root, name, path); }
};

#endif // INIVISITOR_H
