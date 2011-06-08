/*
 * lhs.cpp
 *
 *  Created on: 07-jan-2010
 *      Author: M. El-Kebir
 */

#include "lhs.h"
#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stack>

namespace po = boost::program_options;

/*
 * The lhs program calls lhs::init which calls Params::fromXml to read an LHS parameter file.
 * The calls to readParam functions in fromXml and its sub-functions call the readParam functions
 * here, which read parameter ranges (i.e. [min,max]) rather than individual parameters. The ranges
 * are stored in the _lhsDoubleParam and _lhsIntParam arrays. Any parameters not present in the
 * LHS parameter file have a range of [0,0].
 *
 * Then function Lhs:performLhs is called which uses the range arrays to generate the specified
 * regular parameter files. Function Params::toXml is called for each parameter file to be created.
 * Function toXml only writes parameters to the generated parameter file that had a range explicitly
 * specified in the LHS parameter file and that are known to the Params class as defined parameters
 * (i.e. have an entry in the Params::_description array).
 */

Lhs::Lhs(int nSamples, bool ode)
	: ParamsBase(ode)
	, _nSamples(nSamples)
	, _lhsDoubleParam()
	, _lhsIntParam()
{
}

bool Lhs::readParam(const TiXmlElement* pElement, const TiXmlAttribute* pAttrib,  ParamDoubleType param)
{
	LhsDoubleParam& range = _lhsDoubleParam[param];
	const char* str  = pAttrib->Value();

	char c;
	if (sscanf(str, "[%lf, %lf]%c", &range._min, &range._max, &c) != 2)
	{
		std::cerr << "Value '" << str << "' of attribute '" << pElement->Value() << "/@"
			<< pAttrib->Name() << "' must be a double range ([%lf, %lf])" << std::endl;
		return false;
	}

	if (_description[param].probPos && !(0 <= range._min && range._min <= 1 && 0 <= range._max && range._max <= 1))
	{
		std::cerr << "Values '" << str << "' of attribute '" << pElement->Value() << "/@"
			<< pAttrib->Name() << "' must be in the range [0,1]" << std::endl;
		return false;
	}

	return true;
}

bool Lhs::readParam(const TiXmlElement* pElement, const TiXmlAttribute* pAttrib,  ParamIntType param)
{
	LhsIntParam& range = _lhsIntParam[intIndex(param)];
	const char* str  = pAttrib->Value();

	char c;
	bool pos = _description[param].probPos;
	if (sscanf(str, "[%d, %d]%c", &range._min, &range._max, &c) != 2)
	{
		std::cerr << "Value '" << str << "' of attribute '" << pElement->Value() << "/@"
			<< pAttrib->Name() << "' must be a" << (pos ? " positive" : "n") <<
			" integer in the range ([%d, %d])" << std::endl;
		return false;
	}

	if (pos && !(0 <= range._min && 0 <= range._max))
	{
		std::cerr << "Value '" << str << "' of attribute '" << pElement->Value() << "/@"
			<< pAttrib->Name() << "' must be a positive integer in the range ([%d, %d])" << std::endl;
		return false;
	}

	return true;
}

bool Lhs::init(const char* filename)
{
	if (!fromXml(filename))
	{
		return false;
	}

	// Initialize the range arrays for those parameters that were not present in the parameter file
	// and for which default values are to be used. Without this parameter checks on a generated
	// parameter file might fail.
	for (int i = 0; i < _PARAM_COUNT; i++)
	{
		if (_description[i].useDefault && !_paramsRead[i])
		{
			if (isDouble(i))
			{
				_lhsDoubleParam[i]._min = _description[i].doubleDefault;
				_lhsDoubleParam[i]._max = _description[i].doubleDefault;
			}
			else
			{
				_lhsIntParam[intIndex(i)]._min = _description[i].intDefault;
				_lhsIntParam[intIndex(i)]._max = _description[i].intDefault;
			}
		}
	}

	return true;
}

// Enforce that the CCL5 secretion rate is the same as the CCL2 secretion rate,
// and that the CXCL9 secretion rate is twice the CCL2 rate.
void Lhs::updateParamDouble(ParamDoubleType param, double val)
{
	switch (param)
	{
	case PARAM_MAC_SEC_RATE_CCL2:
		setParam(param, val);
		setParam(PARAM_MAC_SEC_RATE_CCL5, val);
		setParam(PARAM_MAC_SEC_RATE_CXCL9, 2 * val);
		break;
	case PARAM_MAC_SEC_RATE_CCL5:
	case PARAM_MAC_SEC_RATE_CXCL9:
		// ignore
		break;
	default:
		setParam(param, val);
	}
}

void Lhs::performLhs()
{
	const int totalParamCount = PARAM_DOUBLE_COUNT + PARAM_INT_COUNT;

	// first index = param index
	std::vector<std::vector<int> > bins(totalParamCount, std::vector<int>(_nSamples));

	for (int i = 0; i < PARAM_DOUBLE_COUNT; i++)
	{
		for (int j = 0; j < _nSamples; j++)
		{
			bins[i][j] = j;
		}
	}

	for (int i = 0; i < PARAM_INT_COUNT; i++)
	{
		for (int j = 0; j < _nSamples; j++)
		{
			int count = _nSamples / (_lhsIntParam[i]._max - _lhsIntParam[i]._min + 1);
			if (count == 0)
			{
				count = (_lhsIntParam[i]._max - _lhsIntParam[i]._min + 1) / _nSamples;
				int min = j * count + _lhsIntParam[i]._min;
				int max = (j + 1 == _nSamples) ? _lhsIntParam[i]._max + 1 : (j + 1) * count + _lhsIntParam[i]._min;

				bins[i + PARAM_DOUBLE_COUNT][j] = g_Rand.getInt(max, min);
			}
			else
			{
				int val = j / count + _lhsIntParam[i]._min;

				if (val <= _lhsIntParam[i]._max)
				{
					bins[i + PARAM_DOUBLE_COUNT][j] = val;
				}
				else
				{
					bins[i + PARAM_DOUBLE_COUNT][j] =
						g_Rand.getInt(_lhsIntParam[i]._max + 1, _lhsIntParam[i]._min);
				}
			}
		}
	}

	for (int i = 0; i < _nSamples; i++)
	{
		for (int j = 0; j < totalParamCount; j++)
		{
			// pick a random bin, and remove from the vector of available bins
			int binNr = g_Rand.getInt(bins[j].size());
			int k = bins[j][binNr];
			bins[j].erase(bins[j].begin() + binNr);

			if (j < PARAM_DOUBLE_COUNT)
			{
				double a = _lhsDoubleParam[j]._min +
					k * (_lhsDoubleParam[j]._max - _lhsDoubleParam[j]._min) / _nSamples;

				double b = _lhsDoubleParam[j]._min +
					(k + 1) * (_lhsDoubleParam[j]._max - _lhsDoubleParam[j]._min) / _nSamples;

				ParamDoubleType param = (ParamDoubleType) j;
				updateParamDouble(param, g_Rand.getReal(a, b));
			}
			else
			{
				ParamIntType param = (ParamIntType) intIndex(j);
				setParam(param, k);
			}
		}

		bool res = checkParams();
		if (!res)
		{
			std::cerr << "Parameter check failed for generated parameter file " << (i+1) << "." << std::endl;
		}

		// Number the generated parameter files from 1, not 0.
		// Ex. 1.xml, 2.xml,... rather than 0.xml, 1.xml,...
		char buf[1024];
		sprintf(buf, "%d.xml", i+1);
		toXml(buf);
	}
}

bool Lhs::checkParams() const
{
	bool res = true;
	res &= ParamsBase::checkParams();

	// Check that each generated parameter value is in its proper range.
	for (int i = 0; i < _PARAM_COUNT; i++)
	{
		if (isDouble(i))
		{
			if (_doubleParam[i] < _lhsDoubleParam[i]._min || _doubleParam[i] > _lhsDoubleParam[i]._max)
			{
				std::cerr << "Generated value of " << _doubleParam[i] << " for parameter " << _description[i].name
						  << " is outside the specified range of ["  << _lhsDoubleParam[i]._min << ", " << _lhsDoubleParam[i]._max << "]." << std::endl;
				res = false;
			}
		}
		else
		{
			if (_intParam[intIndex(i)] < _lhsIntParam[intIndex(i)]._min || _intParam[intIndex(i)] > _lhsIntParam[intIndex(i)]._max)
			{
				std::cerr << "Generated value of " << _intParam[intIndex(i)] << " for parameter " << _description[i].name
						  << " is outside the specified range of ["  << _lhsIntParam[intIndex(i)]._min << ", " << _lhsIntParam[intIndex(i)]._max << "]." << std::endl;
				res = false;
			}
		}
	}

	return res;
}

void printUsage(char* pArgv0, po::options_description& desc)
{
	std::cout << "Usage: " << pArgv0 << " [options]\n" << desc << std::endl;
}

void printVersion()
{
	std::cout << "Version: " << GR_VERSION << std::endl;
}

int main(int argc, char** argv)
{
	std::string inputFileName;
	int nSamples;
	unsigned long seed;
	bool ode;

	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "Help message")
		("input-file,i", po::value<std::string>(&inputFileName), "Input file name")
		("seed,s", po::value<unsigned long>(&seed)->default_value(1), "Seed")
		("samples,n", po::value<int>(&nSamples), "Number of samples")
		("ode", "Use integrated lymph node ODE for recruitment")
		("version,v", "Version number");

	try
	{
		po::positional_options_description p;
		p.add("input-file", -1);
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
		po::notify(vm);

		if (vm.count("version"))
		{
			printVersion();
			return 0;
		}

		if (vm.count("help"))
		{
			printUsage(argv[0], desc);
			return 0;
		}

		ode = vm.count("ode");

		if (!vm.count("samples"))
		{
			std::cerr << "Missing number of samples" << std::endl;
			return 1;
		}

		if (!vm.count("input-file"))
		{
			std::cerr << "Missing input file" << std::endl;
			return 1;
		}
	}
	catch (std::exception& e)
	{
		std::cerr << e.what() << std::endl;
		return 1;
	}

	g_Rand.setSeed(seed);

	Lhs lhs(nSamples);
	if (!lhs.init(inputFileName.c_str()))
	{
		return 1;
	}

	lhs.performLhs();

	return 0;
}
