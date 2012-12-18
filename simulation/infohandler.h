#ifndef INFOVISITOR_H
#define INFOVISITOR_H

#include "params.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <iostream>

/// Simple handler for info formatted files
class INFOHandler : public ParamFileHandler {
public:
  explicit INFOHandler(const std::string& _root) : ParamFileHandler(_root) {}
  void read(std::istream& file, ptree& pt) const
    { boost::property_tree::info_parser::read_info(file, pt); }
  void write(std::ostream& file, const ptree& pt) const
    { boost::property_tree::info_parser::write_info(file, pt); }
  static std::string getPath(const std::string& root, const std::string& name, const std::string& path = "")
    { return ParamFileHandler::getPath(root, name, path); }
};

#endif // INFOVISITOR_H
