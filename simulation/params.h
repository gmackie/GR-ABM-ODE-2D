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

  static Params* getInstance(const Pos& dim);
  static Params* getInstance();
  static bool reinit(const char* filename);

  bool fromXml(const char* filename);

protected:
  Params(const Pos& dim);
  virtual ~Params();

  bool readParam(const TiXmlElement* pElement, const TiXmlAttribute* pAttrib,  ParamDoubleType param);
  bool readParam(const TiXmlElement* pElement, const TiXmlAttribute* pAttrib,  ParamIntType param);

private:
  static Params* _pInstance;

};

inline Params* Params::getInstance()
{
  assert(_pInstance != NULL);
  return _pInstance;
}

inline Params* Params::getInstance(const Pos& dim)
{
  if (!_pInstance)
    {
      assert(dim.x>0 && dim.y>0);
      _pInstance = new Params(dim);
    }

  return _pInstance;
}

#endif /* PARAMS_H_ */
