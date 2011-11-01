/*
 * lhs.h
 *
 *  Created on: 07-jan-2010
 *      Author: M. El-Kebir
 */

#ifndef LHS_H
#define LHS_H

#include "paramsbase.h"
#include <vector>

struct LhsDoubleParam
{
	double _min;
	double _max;
};

struct LhsIntParam
{
	int _min;
	int _max;
};

class Lhs : public ParamsBase
{
private:

	int _nSamples;
	LhsDoubleParam _lhsDoubleParam[PARAM_DOUBLE_COUNT];
	LhsIntParam _lhsIntParam[PARAM_INT_COUNT];

	bool readParam(const TiXmlElement* pElement, const TiXmlAttribute* pAttrib,  ParamDoubleType param);
	bool readParam(const TiXmlElement* pElement, const TiXmlAttribute* pAttrib,  ParamIntType param);

	void updateParamDouble(ParamDoubleType param, double val);

	bool checkParams() const;

public:
	Lhs(int nSamples, bool ode = false);
	bool init(const char* filename);
	void performLhs();

};

#endif /* LHS_H */
