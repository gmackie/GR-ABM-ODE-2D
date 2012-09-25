/*
 * Serialization.cpp
 *
 *  Created on: Jun 17, 2011
 *      Author: pwolberg
 */

#include "serialization.h"

const std::string Serialization::_HeaderSuffix = "_Start";
const std::string Serialization::_FooterSuffix = "_End";

Serialization::Serialization()
{

}

Serialization::~Serialization()
{
}

void Serialization::writeHeader(std::ostream& out, std::string className)
{
  out << className << _HeaderSuffix << std::endl;
}

void Serialization::writeFooter(std::ostream& out, std::string className)
{
  out << className << _FooterSuffix << std::endl;
}

bool Serialization::readHeader(std::istream& in, std::string className)
{
  return readHeaderFooter(in, className, className + _HeaderSuffix, "header");
}

bool Serialization::readFooter(std::istream& in, std::string className)
{
  return readHeaderFooter(in, className, className + _FooterSuffix, "footer");
}

bool Serialization::readHeaderFooter(std::istream& in, std::string className, std::string expectedText, std::string type)
{
  std::string readText;

  try
    {
      in >> readText;
    }
  catch(std::exception& e)
    {
      std::cerr << "Caught an exception on reading " << type << " " << expectedText << std::endl;
      std::cerr<<e.what()<<std::endl;
    }

  if (readText != expectedText)
    {
      std::cerr << "Error deserializing " << className << ". The deserialized " << type << ", '" << readText << "'"
                " does not match the expected " << type << " of '" << expectedText << "'" << std::endl;
      std::cerr << "Check the " << className << " class to be sure that required members are serialized and deserialized"
                << std::endl
                << "and that each member which is serialized is deserialized and in the same order as it was serialized."
                << std::endl;
      throw std::ios_base::failure("Failed to deserialize");
      return false;
    }

  return true;
}
