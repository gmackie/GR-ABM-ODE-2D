/*
 * Params.cpp
 *
 *  Created on: Jun 3, 2011
 *      Author: pwolberg
 */

#include "params.h"
#include <iostream>
#include <fstream>
#include <stdio.h>

Params* Params::_pInstance = NULL;

Params::Params(bool ode)
	: ParamsBase(ode)
{
	// Make sure _pInstance is defined when instantiating a Params object,
	// even if not using getInstance or reinit. This could happen if
	// instantiating a sub-class of Params.
	Params::_pInstance = this;
}

Params::~Params()
{
}

bool Params::reinit(const char* filename)
{
	Params* pNewInstance = new Params(false);
	if (!pNewInstance->fromXml(filename))
	{
		delete pNewInstance;
		return false;
	}
	else
	{
		delete _pInstance;
		_pInstance = pNewInstance;
		return true;
	}
}

bool Params::fromXml(const char* filename)
{
	if(!ParamsBase::fromXml(filename)) return false;  //Failed to read from file

	computeParams();

  return checkParams(); //Failed? to validate
}

bool Params::readParam(const TiXmlElement* pElement, const TiXmlAttribute* pAttrib,  ParamDoubleType param)
{
	double* pVar = _doubleParam + param;

	if(pAttrib->QueryDoubleValue(pVar) == TIXML_WRONG_TYPE)
	{
		std::cerr << "Value '" << pAttrib->Value() << "' of attribute '" << pElement->Value() << "/@"
			<< pAttrib->Name() << "' must be a double. " << std::endl;
		return false;
	}

	if (_description[param].probPos && !(0 <= *pVar && *pVar <= 1))
	{
		std::cerr << "Value '" << pAttrib->Value() << "' of attribute '" << pElement->Value() << "/@"
			<< pAttrib->Name() << "' must be in the range [0,1]" << std::endl;
		return false;
	}

	return true;
}

bool Params::readParam(const TiXmlElement* pElement, const TiXmlAttribute* pAttrib,  ParamIntType paramIntIndex)
{
	int* pVar = _intParam + paramIntIndex;
	bool pos = _description[paramIndex(paramIntIndex)].probPos;

	if(pAttrib->QueryIntValue(pVar) == TIXML_WRONG_TYPE)
	{
		std::cerr << "Value '" << pAttrib->Value() << "' of attribute '" << pElement->Value() << "/@"
			<<  pAttrib->Name()  << "' must be a" << (pos ? "n " : " positive ")
			<< "integer" << std::endl;
		return false;
	}

	if (pos && *pVar < 0)
	{
		std::cerr << "Value '" << pAttrib->Value() << "' of attribute '" << pElement->Value() << "/@"
			<<  pAttrib->Name()  << "' must be a positive integer" << std::endl;
		return false;
	}

	return true;
}
