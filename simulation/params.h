/*
 * Params.h
 *
 *  Created on: Jun 3, 2011
 *      Author: pwolberg
 */

#ifndef PARAMS_H_
#define PARAMS_H_

#include "paramsbase.h"

// introduce a convenient shorthand
#define _PARAM(a) Params::getInstance()->getParam((a))

class Params : public ParamsBase
{
public:

	static Params* getInstance();
	static bool reinit(const char* filename);

	bool fromXml(const char* filename);

protected:
	Params();
	virtual ~Params();

	bool readParam(const TiXmlElement* pElement, const TiXmlAttribute* pAttrib,  ParamDoubleType param);
	bool readParam(const TiXmlElement* pElement, const TiXmlAttribute* pAttrib,  ParamIntType param);

private:
	static Params* _pInstance;

};

inline Params* Params::getInstance()
{
	if (!_pInstance)
		_pInstance = new Params();

	return _pInstance;
}

#endif /* PARAMS_H_ */
