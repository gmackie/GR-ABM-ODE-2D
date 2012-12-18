#ifndef JSONHANDLER_H
#define JSONHANDLER_H

#include "params.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <iostream>

/// Simple handler for json formatted files
class JSONHandler : public ParamFileHandler {
public:
  explicit JSONHandler(const std::string& _root) : ParamFileHandler(_root) {}
  void read(std::istream& file, ptree& pt) const
    { boost::property_tree::json_parser::read_json(file, pt); }
  void write(std::ostream& file, const ptree& pt) const
    { boost::property_tree::json_parser::write_json(file, pt); }
  static std::string getPath(const std::string& root, const std::string& name, const std::string& path = "")
    { return ParamFileHandler::getPath(root, name, path); }
};

#endif // JSONHANDLER_H
