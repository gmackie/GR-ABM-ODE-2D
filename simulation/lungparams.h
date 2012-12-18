#ifndef LUNGPARAM_H
#define LUNGPARAM_H
#include "params.h"
/// Singleton version of Params
class LungParam : public Params {
private:
  LungParam() {}                    //Do not override
  LungParam(const LungParam&);      //Do not define
  void operator=(const LungParam&); //Do not define
public:
  static LungParam& getInstance()
  {
    static LungParam instance;
    return instance;
  }
};

#define _PARAM(E) (LungParam::getInstance().get##E())

#endif // LUNGPARAM_H
