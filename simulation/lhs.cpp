/*
 * lhs.cpp
 *
 *  Created on: 07-jan-2010
 *      Author: M. El-Kebir
 */

#include "lhs.h"
#include "params.h"
#include <boost/program_options.hpp>
#include <stdio.h>

namespace po = boost::program_options;

Lhs::Lhs(int nSamples, bool ode)
	: Params(false, ode)
	, _nSamples(nSamples)
	, _lhsDoubleParam()
	, _lhsIntParam()
{
}

bool Lhs::readParam(const TiXmlElement* pElement, ParamDoubleType param, bool prob)
{
	const char* paramName = getName(param);
	const char* str = pElement->Attribute(paramName);
	if (!str)
	{
		std::cerr << "Expected attribute '" << pElement->Value() << "/@"
			<< paramName << "'" << std::endl;
		return false;
	}

	LhsDoubleParam& range = _lhsDoubleParam[param];

	char c;
	if (sscanf(str, "[%lf, %lf]%c", &range._min, &range._max, &c) != 2)
	{
		std::cerr << "Value of attribute '" << pElement->Value() << "/@"
			<< paramName << "' must be a double range ([%lf, %lf])" << std::endl;
		return false;
	}

	if (prob && !(0 <= range._min && range._min <= 1 && 0 <= range._max && range._max <= 1))
	{
		std::cerr << "Values of attribute '" << pElement->Value() << "/@"
			<< paramName << "' must be in the range [0,1]" << std::endl;
		return false;
	}

	return true;
}

// Ignore unspecified parameters with default values.
bool Lhs::readParam(const TiXmlElement* pElement, ParamDoubleType param, double defaultVal, bool prob)
{
	const char* paramName = getName(param);
	if (pElement->Attribute( paramName ))
	{
		bool res = readParam(pElement, param, prob);
		return res;
	}
	else
	{
		return true;
	}
}

bool Lhs::readParam(const TiXmlElement* pElement, ParamIntType param, bool pos)
{
	const char* paramName = getName(param);
	const char* str = pElement->Attribute(paramName);
	if (!str)
	{
		std::cerr << "Expected attribute '" << pElement->Value() << "/@"
			<< paramName << "'" << std::endl;
		return false;
	}

	LhsIntParam& range = _lhsIntParam[param];

	char c;
	if (sscanf(str, "[%d, %d]%c", &range._min, &range._max, &c) != 2)
	{
		std::cerr << "Value of attribute '" << pElement->Value() << "/@"
			<< paramName << "' must be a" << (pos ? " positive" : "n") <<
			" integer in the range ([%d, %d])" << std::endl;
		return false;
	}

	if (pos && !(0 <= range._min && 0 <= range._max))
	{
		std::cerr << "Value of attribute '" << pElement->Value() << "/@"
			<< paramName << "' must be a" << (pos ? " positive" : "n") <<
			" integer in the range ([%d, %d])" << std::endl;
		return false;
	}

	return true;
}

// Ignore unspecified parameters with default values.
bool Lhs::readParam(const TiXmlElement* pElement, ParamIntType param, int defaultVal, bool pos)
{
	const char* paramName = getName(param);
	if (pElement->Attribute( paramName ))
	{
		bool res = readParam(pElement, param, pos);
		return res;
	}
	else
	{
		return true;
	}
}

bool Lhs::init(const char* filename)
{
	return fromXml(filename);
}

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
				ParamIntType param = (ParamIntType) (j - PARAM_DOUBLE_COUNT);
				setParam(param, k);
			}
		}

		if (!getUseRecruitmentWeights())
			updateRecruitmentWeights();

		// Number the parameter files wrttin from 1, not 0.
		// Ex. 1.xml, 2.xml,... rather than 0.xml, 1.xml,...
		char buf[1024];
		sprintf(buf, "%d.xml", i+1);
		toXml(buf);
	}
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
