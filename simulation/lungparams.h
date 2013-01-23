#ifndef LUNGPARAM_H
#define LUNGPARAM_H

#include "params.h"

/// Singleton version of Params
class LungParam : public Params {
private:
  LungParam(const Pos& dim) : Params(dim) {}      //Do not override
  LungParam(const LungParam&);      //Do not define
  void operator=(const LungParam&); //Do not define

  static LungParam* _pInstance;

public:
  static LungParam* getInstance(const Pos& dim=Pos(-1,-1));

};

// Keep a pointer to the singleton rather than a static local variable
// because we need to pass in the simulation dimensions for defining
// some parameters based on dimension, in function computeParams.
inline LungParam* LungParam::getInstance(const Pos& dim)
{
  if (!_pInstance) {
    assert(dim.x>0 && dim.y>0);
    _pInstance = new LungParam(dim);
  }

  return _pInstance;
}

#define _PARAM(E) (LungParam::getInstance()->get##E())
#define _PARAMDESC(E) (LungParam::getInstance()->get##E##Desc())


#endif // LUNGPARAM_H
