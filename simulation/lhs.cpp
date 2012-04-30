/*
 * lhs.cpp
 *
 *  Created on: 07-jan-2010
 *      Author: M. El-Kebir
 */

#include "lhs.h"
#include "rand.h"
#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <stdio.h>

Rand g_Rand(1337);
namespace po = boost::program_options;

/*
 * The lhs program implements the Latin Hypercube Sampling algorithm.
 * This algorithm generates parameter sets by choosing parameter values
 * at random from specified probability distributions for each parameter
 * to be varied. See Mckay, M.D., Beckman, R.J., Conover, W.J., 1979.
 * "Comparison of 3 methods for selecting values of input variables in
 * the analysis of output from a computer code." Technometrics 21 (2), 239–245,
 * http://www.jstor.org/stable/1268522
 *
 * The general algorithm places no restrictions on the probability
 * distributions used and a different distribution can be used for each parameter.
 * This implementation uses a uniform distribution for each for each parameter.
 *
 * The lhs program has 2 command line arguments, the number of samples (number of
 * regular parameter files to be generated) and the name of an LHS parameter file.
 *
 * An LHS parameter file is an XML file that contains ranges for some or all parameters.
 * A range is specified as [min,max], ex. thresholdApoptosisTNF = "[0.0525,0.0725]".
 * If a parameter has a single value specified instead of a range then that value
 * is used as the min and max for the range.
 * If a parameter does not appear in an LHS parameter file then its range is [0,0].
 *
 * For each parameter, its range is divided into N sub-ranges, where N is the number
 * of samples. For each parameter file to be created, for each parameter a sub-range for the
 * parameter is chosen at random and then a value is chosen at random from within that
 * sub-range. That sub-range for that parameter is then no longer used (sampling without
 * replacement). The sub-ranges are chosen at random rather than iterating over them
 * because we don't want to create parameter files such that the first one has values
 * only from the first sub-range for each parameter, the second parameter file only
 * has values from the 2nd sub-range for each parameter, etc. We want a random mix
 * of values from among the sub-ranges.
 *
 * Consider the following simple example, where we are only using 2 parameters,
 * have a sample size of 4 and have the following LHS parameter file (the XML tags are
 * left out for simplicity):
 *
 * p1 = "[1,100]
 * p2 = "[1000,2000]"
 *
 * This produces the following set of sub-ranges for the parameters.
 *
 * p1	1-25		26-50		51-75		76-100
 * p2	1000-1250	1251-1500	1501-1750	1751-2000
 *
 * 4 parameter files are created, as follows:
 *
 * Parameter File 1:
 *
 * For p1, pick a sub-range at random, which is equivalent to choosing a random number
 * between 1 and 4. Say its 2, so the sub-range is 26-50. Then pick a number at random
 * from this sub-range, say 37. We are no longer going to use this sub-range for
 * generating subsequent parameter files, so our sets of sub-ranges are:
 * p1	1-25		51-75		76-100
 * p2	1000-1250	1251-1500	1501-1750	1751-2000
 *
 * Now do the same for parameter p2. Say the sub-range is 3 and the value is 1508.
 * The remaining sub-ranges are:
 * p1	1-25		51-75		76-100
 * p2	1000-1250	1251-1500	1751-2000
 *
 * The first generated parameter file is:
 * p1 = "37"
 * p2 = "1508"
 *
 * We repeat the process for the next parameter file, but now we choose a random number between
 * 1 and 3, rather than 1 and 4, to choose the sub-range for each parameter. Say we choose
 * sub-range 1 for p1, and value 13 from that sub-range, and sub-range 3 for p2, and value
 * 1974 from that sub-range. The remaining sub-ranges and the generated parameter file are:
 *
 * p1	51-75		76-100
 * p2	1000-1250	1251-1500
 *
 * The 2nd generated parameter file is:
 * p1 = "13"
 * p2 = "1974"
 *
 * We would repeat this process twice more to generate the remaining 2 parameter files.
 *
 * If we were to iterate over the sub-ranges, rather than pick them at random for each parameter,
 * we might end up with parameter files like:
 *
 * p1 = "13"
 * p2 = "1089"
 *
 * p1 = "42"
 * p2 = "1422"
 *
 * p1 = "67"
 * p2 = "1598"
 *
 * p1 = "77"
 * p2 = "1905"
 *
 * The first parameter file would have parameters with values all taken from the first sub-range,
 * the 2nd parameter file with values all taken from the 2nd sub-range, etc. This would be a less
 * random sampling of the parameter space.
 *
 * The lhs program calls lhs::init which calls Params::fromXml to read an LHS parameter file.
 * The calls to readParam functions in fromXml and its sub-functions call the readParam functions
 * here, in the Lhs class, which read parameter ranges (i.e. [min,max]) rather than individual
 * parameters. The ranges are stored in the _lhsDoubleParam and _lhsIntParam arrays. Any
 * parameters not present in the LHS parameter file have a range of [0,0].
 *
 * Then function Lhs:performLhs is called which uses the range arrays to generate the specified
 * regular parameter files. Function Params::toXml is called for each parameter file to be written.
 * Params::toXml only writes parameters to the generated parameter file that had a value explicitly
 * specified in the LHS parameter file and that are known to the Params class as defined parameters
 * (i.e. have an entry in the Params::_description array).
 */

Lhs::Lhs(int nSamples, Pos dim)
	: ParamsBase(dim)
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
		double val;
		if (sscanf(str, "%lf", &val) == 1)
		{
			range._min = val;
			range._max = val;
                        
		}
		else
		{
			std::cerr << "Value '" << str << "' of attribute '" << pElement->Value() << "/@"
				<< pAttrib->Name() << "' must be a double range ([%lf, %lf])" << std::endl;
			return false;
		}
	}

	if (range._min == range._max)
	{
		range._isRange = false;
	}
	else
	{
		range._isRange = true;
	}

	if (range._min > range._max)
	{
		std::cerr << "Values '" << str << "' of attribute '" << pElement->Value() << "/@"
			<< pAttrib->Name() << "' have a min of " << range._min << " which is greater than the max of " << range._max << std::endl;
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

bool Lhs::readParam(const TiXmlElement* pElement, const TiXmlAttribute* pAttrib,  ParamIntType paramIntIndex)
{
	LhsIntParam& range = _lhsIntParam[paramIntIndex];
	const char* str  = pAttrib->Value();

	char c;
	bool pos = _description[paramIndex(paramIntIndex)].probPos;

	if (sscanf(str, "[%d, %d]%c", &range._min, &range._max, &c) != 2)
	{
		int val;
		if (sscanf(str, "%d", &val) == 1)
		{
			range._min = val;
			range._max = val;
            
		}
		else
		{
			std::cerr << "Value '" << str << "' of attribute '" << pElement->Value() << "/@"
				<< pAttrib->Name() << "' must be a" << (pos ? " positive" : "n") <<
				" integer in the range ([%d, %d])" << std::endl;
			return false;
		}
	}

	if (range._min == range._max)
	{
		range._isRange = false;
	}
	else
	{
		range._isRange = true;
	}

	if (range._min > range._max)
	{
		std::cerr << "Values '" << str << "' of attribute '" << pElement->Value() << "/@"
			<< pAttrib->Name() << "' have a min of " << range._min << " which is greater than the max of " << range._max << std::endl;
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
	// Need to initialize _isRange here for parameters that don't appear in the LHS parameter file.
	// They will have a single value defined, from the default value in paramsbase.cpp.
	for (int i = 0; i < PARAM_DOUBLE_COUNT; ++i)
	{
		_lhsDoubleParam[i]._isRange = false;
	}

	for (int i = 0; i < PARAM_INT_COUNT; ++i)
	{
		_lhsIntParam[i]._isRange = false;
	}

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

/*
 * Perform the LHS
 *
 * See the comments at the start of this file for a general description of the LHS
 * algorithm and a reference.
 *
 * This implementation has some details beyond the general algorithm, based on the
 * particular data structures used. In a particular, parameters with integer values are
 * handled in a special way.
 *
 * The bins data structure is a dynamic 2D array, implemented as a vector of vectors of ints,
 * with one row per parameter and one column per LHS sample. It is initialized so that each
 * element of the array contains its sample ordinal, i.e. each column of the array has
 * the sample ordinal as its value - all elements of column 0 are 0, of column 1 are 1, etc.
 * Vectors are used because elements of the array will be deleted as they are used.
 *
 * The bins for all the integer valued parameters are defined as follows. If the number of
 * samples is less than the number of integers between the min and max for that parameter,
 * then each sample is assigned an equal sub-range and a random value from that sub-range
 * is assigned to that bin. For example, if a parameter's range is [1, 100] and
 * _nSamples = 10, then the 1st sub-range is [1,10] and the first bin for this parameter
 * will be a random number between 1 and 10. The 2nd bin for this parameter is a random value
 * in the sub-range [11,20], etc.
 *
 * If _nSamples >= the number of integers between the min and max for a parameter then
 * several samples must be assigned the same integer. For example, if a parameter's
 * range is [1,100] and _nSamples = 1000, then the number of samples that will have the
 * same parameter value = 1000/100 = 10. The samples 1 to 10 all have a parameter
 * value of 1, the samples 11 to 20 all have a parameter value of 2, etc.
 *
 * Once this initialization of the bins is complete, one parameter file is defined
 * for each sample. For each parameter for a sample a value is determined and stored in the
 * parameter data of the Lhs class base class (ParamsBase), using the ParamsBase::setParam
 * function. Then checkParams is used to check the generated parameters to be sure that
 * each parameter value is within its specified range and that any other constraints
 * are met, and then ParamBase::toXml is used to write the parameters to a file.
 * The file name is the sample number with an extension of ".xml", ex. 1.xml, 2.xml, etc.
 *
 * The value for a parameter is determined as follows. A bin is chosen at random from
 * that parameter's bin vector. Note that if constructing the parameters for sample i, the
 * parameter's ith bin is NOT necessarily used. Otherwise sample i would have values
 * only from bin i for all its parameters, which would not be  as well randomized.
 * If the parameter is an integer parameter then the value from the selected bin is used.
 * If the parameter is not an integer parameter, then the bin contains an ordinal value k
 * which is used to define the kth sub-range of the range for the parameter. A random
 * value is chosen from the sub-range. The bin is deleted from the vector of bins for this
 * parameter so it is not reused (random selection without replacement).
 *
 * A couple of further things to note.
 *
 * This function processes each parameter defined in ParamsBase, whether or not that
 * parameter was present in the LHS parameter file. If a parameter is not present in
 * the LHS parameter file its range is [0,0].
 *
 * This function also does not restrict itself to those parameters that have a range where the min < max,
 * that is an actual range with more than one value. If a parameter has a range where min = max,
 * (because that is how it was specified in the parameter file or because it was not specified so min = max = 0)
 * then the generated values for that parameter will be that specific value for each generated
 * parameter file.
 *
 * Despite generating a value for each parameter defined in ParamsBase, parameter values are written to a
 * generated parameter file for only those parameters that were present in the LHS parameter file. This
 * is how ParamsBase::toXml works.
 *
 */
void Lhs::performLhs(bool logscale)
{
	const int totalParamCount = PARAM_DOUBLE_COUNT + PARAM_INT_COUNT;

	// first index = param index
	// bins is a vector of vectors. Essentially it is a dynamic 2D array with one row per
	// parameter and one column per LHS sample. Initialize each element of the bin 2D array
	// to be its sample ordinal.
	std::vector<std::vector<int> > bins(totalParamCount, std::vector<int>(_nSamples));

	// Initialize the bins for each parameter to the sample ordinal:
	// bin 0 for a parameter is 0, bin 1 is 1, etc.
	for (int i = 0; i < PARAM_DOUBLE_COUNT; i++)
	{
		for (int j = 0; j < _nSamples; j++)
		{
			bins[i][j] = j;
		}
	}

	for (int i = 0; i < PARAM_INT_COUNT; i++)
	{
		int paramindex = paramIndex((ParamIntType) i);

		for (int j = 0; j < _nSamples; j++)
		{
			if (!_lhsIntParam[i]._isRange)
			{
				bins[paramindex][j] = _lhsIntParam[i]._min;
				continue;
			}

			int count = _nSamples / (_lhsIntParam[i]._max - _lhsIntParam[i]._min + 1);
			if (count == 0)
			{
				// The number of samples < number of integers in the parameter range.
				// Divide the parameter range into sub-ranges and pick a random value
				// from each sub-range. Note that max is one greater than the max
				// for a sub-range because getInt(max, min) chooses a value in the range [min, max-1].
				count = (_lhsIntParam[i]._max - _lhsIntParam[i]._min + 1) / _nSamples;
				int min = j * count + _lhsIntParam[i]._min;
				int max = (j + 1 == _nSamples) ? _lhsIntParam[i]._max + 1 : (j + 1) * count + _lhsIntParam[i]._min;

				bins[paramindex][j] = g_Rand.getInt(max, min);
			}
			else
			{
				// The number of samples >= number of integers in the parameter range.
				// Repeat each integer value for count number of samples.
				int val = j / count + _lhsIntParam[i]._min;

				if (val <= _lhsIntParam[i]._max)
				{
					bins[paramindex][j] = val;
				}
				else
				{
					bins[paramindex][j] =
						g_Rand.getInt(_lhsIntParam[i]._max + 1, _lhsIntParam[i]._min);
				}
			}
		}
	}

	// Create a parameter file for each sample.
	for (int i = 0; i < _nSamples; i++)
	{
		for (int j = 0; j < totalParamCount; j++)
		{
			// pick a random bin, and remove from the vector of available bins
			int binNr = g_Rand.getInt(bins[j].size());
			int k = bins[j][binNr];
            
//            std::cout << binNr << "   " << k << std::endl;
            
			bins[j].erase(bins[j].begin() + binNr);

			if (j < PARAM_DOUBLE_COUNT)
			{
				ParamDoubleType param = (ParamDoubleType) j;

				if (_lhsDoubleParam[j]._isRange)
				{
					// A non-integer variable. Pick a random value between the min and max for the sub-range
					// for the selected bin.
                    // If the log scale option is on then the random value is chosen on a log scale
        
                    
                    if (logscale) {
                        
                        if (_lhsDoubleParam[j]._min <= 0)
                        {
                            std::cerr << "Minimum Value of a Paramater is either negative or zero. Log-scale cannot be used in this case " << std::endl;
                            exit(1);
                        }
                        
                        double log_min = log10(_lhsDoubleParam[j]._min);
                        double log_max = log10(_lhsDoubleParam[j]._max);
                        double a = log_min + k * (log_max - log_min) / _nSamples;
                        double b =log_min + (k + 1) * (log_max - log_min) / _nSamples;
                        double convert = pow(10,g_Rand.getReal(a, b));

                        updateParamDouble(param, convert);
                    }
                    
                    
                    else
                    {
                        double a = _lhsDoubleParam[j]._min + k * (_lhsDoubleParam[j]._max - _lhsDoubleParam[j]._min) / _nSamples;
                        double b = _lhsDoubleParam[j]._min + (k + 1) * (_lhsDoubleParam[j]._max - _lhsDoubleParam[j]._min) / _nSamples;
                        
                        updateParamDouble(param, g_Rand.getReal(a,b));
                        
                    }

				}
				else
				{
					// This isn't a range, i.e. min == max, so just use the single value for this parameter.
					updateParamDouble(param, _lhsDoubleParam[j]._min);
				}
			}
			else
			{
				// An integer parameter. Use the integer value from the randomly selected bin.
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
			// Don't check parameters calculated in updateParamDouble since they might not match any range specified.
			if ( ((i != PARAM_MAC_SEC_RATE_CCL5) && (i != PARAM_MAC_SEC_RATE_CXCL9))  && (_doubleParam[i] < _lhsDoubleParam[i]._min || _doubleParam[i] > _lhsDoubleParam[i]._max) )
			{
				std::cerr << "Generated value of " << _doubleParam[i] << " for parameter " << _description[i].name
						  << " is outside the specified range of ["  << _lhsDoubleParam[i]._min << ", " << _lhsDoubleParam[i]._max << "]." << std::endl;
				res = false;
			}
		}
		else
		{
			if ((_intParam[intIndex(i)] < _lhsIntParam[intIndex(i)]._min || _intParam[intIndex(i)] > _lhsIntParam[intIndex(i)]._max))
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
	int dim;
	bool logscale = false;

	/* set seed to current time, in case not specified */
	time_t curTime;
	time(&curTime);

	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "Help message")
		("input-file,i", po::value<std::string>(&inputFileName), "Input file name")
		("seed,s", po::value<unsigned long>(&seed)->default_value((unsigned long) curTime), "Seed")
		("samples,n", po::value<int>(&nSamples), "Number of samples")
		("dim,d", po::value(&dim)->default_value(100), "Size of simulation grid")
		("version,v", "Version number")
        ("log-scale", "Use log scale intervals for double-type parameters");

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
        if (vm.count("log-scale")) {
            logscale = 1;
        }
	}
	catch (std::exception& e)
	{
		std::cerr << e.what() << std::endl;
		return 1;
	}

	g_Rand.setSeed(seed);

	Lhs lhs(nSamples, Pos(dim, dim));
	if (!lhs.init(inputFileName.c_str()))
	{
		return 1;
	}

	lhs.performLhs(logscale);

	return 0;
}
