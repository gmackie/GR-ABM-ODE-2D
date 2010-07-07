/*
 * lhs.h
 *
 *  Created on: 07-jan-2010
 *      Author: M. El-Kebir
 */

#ifndef LHS_H
#define LHS_H

#include "tinyxml/tinyxml.h"
#include "params.h"

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

class Lhs : public Params
{
private:
	int _nSamples;
	LhsDoubleParam _lhsDoubleParam[PARAM_DOUBLE_COUNT];
	LhsIntParam _lhsIntParam[PARAM_INT_COUNT];
	bool readParam(const TiXmlElement* pElement, ParamDoubleType param, bool prob);
	bool readParam(const TiXmlElement* pElement, ParamIntType param, bool pos);
	void updateParamDouble(ParamDoubleType param, double val);

public:
	Lhs(int nSamples, bool ode = false);
	bool init(const char* filename);
	void performLhs();
};

#endif /* LHS_H */
