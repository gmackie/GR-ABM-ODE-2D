#include "params.h"
Params::Params() :
#define P(type, n, p, defau, d, u, min, max)  \
  _##p##_##n(type(defau)), \
  _##p##_##n##Desc(type(min), type(max), #n, boost::optional<type>(defau), std::string(d), std::string(u), make_path(#p)),
#include "params.def"
  dummy()
{

}
/* Called from constructor */
/*void Params::initialize() {
#define P(type, n, p, defau, d, u, min, max)  \
 { if(type(min) != type(max)) \
   _##p##_##n##Desc.range = Range<type>(type(min), type(max));             \
 _##p##_##n##Desc.name = std::string(#n);           \
 _##p##_##n##Desc.def = boost::optional<type>(defau); \
 _##p##_##n##Desc.desc = std::string(d);           \
 _##p##_##n##Desc.units = std::string(u);         \
 _##p##_##n##Desc.path = make_path(#p);           \
 _##p##_##n = !_##p##_##n##Desc.def ? type() : *(_##p##_##n##Desc.def); }
#include "params.def"
}*/

