#ifndef XMLHANDLER_H
#define XMLHANDLER_H

#include "params.h"
#include <boost/property_tree/xml_parser.hpp>


/// More complicated handler for XML files
class XMLHandler : public ParamFileHandler {
  boost::property_tree::xml_parser::xml_writer_settings<std::string::value_type> settings;
  void write_xml(std::istream& file, ptree& pt);
public:
  explicit XMLHandler(const std::string& root) : ParamFileHandler(root), settings(' ', 4) {}
  XMLHandler(const std::string& root, boost::property_tree::xml_parser::xml_writer_settings<std::string::value_type>& _settings) : ParamFileHandler(root), settings(_settings) {}

  void read(std::istream& file, ptree& pt) const {
    boost::property_tree::xml_parser::read_xml(file, pt);
  }

  static std::string getPath(const std::string& root, const std::string& name, const std::string& path) {
    return root + '.' + (path.empty() ? std::string() : (path + ".")) + "<xmlattr>." + name;
  }

  std::string getFullPath(const std::string& root, const std::string &name, const std::string &path) const {
    return XMLHandler::getPath(root, name, path);
  }

  void write(std::ostream& file, const ptree& pt) const {
    file<<"<?xml version=\"1.0\" encoding=\"utf-8\"?>"<<std::endl  ; //Add this to ensure proper encoding and validation
    this->write_xml(file, pt);    //Don't use the boost version, we reimplement our own for proper formatting
  }
  void write_attrs(std::ostream& file, const ptree& pt, size_t indent) const;
  bool write_xml(std::ostream& file, const ptree& pt, size_t indent=0) const;
};

inline void XMLHandler::write_attrs(std::ostream& file, const ptree & pt, size_t indent) const {
  using namespace boost::property_tree;
  ptree::const_iterator end = pt.end();
  const char eol = (pt.size() > 3 ? '\n' : ' ');
  file << eol;
  for(ptree::const_iterator it = pt.begin(); it!=end; it++)
  {
    if(eol == '\n')
      for(size_t i=0;i<indent*settings.indent_count;i++) file << settings.indent_char;
    file << it->first << " = \"" << it->second.data() << '"' << eol;
  }
}

inline bool XMLHandler::write_xml(std::ostream& file, const ptree & pt, size_t indent) const {
  using namespace boost::property_tree;
  boost::optional<const ptree&> attrs = pt.get_child_optional("<xmlattr>");
  if(!!attrs) { write_attrs(file, *attrs, indent); }
  bool empty = true;
  ptree::const_iterator end = pt.end();
  for(ptree::const_iterator it = pt.begin(); it!=end; it++) {
    if(it->first == "<xmlattr>") continue;  //Already handled in write_attrs
    if(empty && indent != 0) {
        if(attrs)
          for(size_t i=0;i<indent*settings.indent_count;i++) file << settings.indent_char;
        file << ">" << std::endl; empty = false;
    }
    for(size_t i=0;i<indent*settings.indent_count;i++) file << settings.indent_char;
    if(it->first == "<xmlcomment>")   //Handle comments
      file << "<!--" << it->second.data() << "-->" << std::endl;
    else if(it->first == "<xmltext>") { //Handle inner text
      file << it->second.data();
    }
    else {
      file << '<' << it->first;
      if(!write_xml(file, it->second, indent+1))  //Recurse!
      {
        for(size_t i=0;i<indent*settings.indent_count;i++) file << settings.indent_char;
        file << "</" << it->first << '>' << std::endl;
      }
    }
  }
  if(empty && indent != 0) {
      for(size_t i=0;i<indent*settings.indent_count;i++) file << settings.indent_char;
      file << "/>" <<std::endl;
  }
  return empty;
}

#endif // XMLHANDLER_H
