#ifndef VERSION_H_
#define VERSION_H_

#define _NUM(x) 0##x
#define NUM(x) _NUM(x)

#if NUM(SVN_VERSION) == 0
#undef SVN_VERSION  //Invalid svn_version
#endif

#ifndef SVN_VERSION
#define GR_VERSION "12022010"
#else

#define _STRINGIFY(x) #x
#define STRINGIFY(x) _STRINGIFY(x)

#define GR_VERSION "r" STRINGIFY(SVN_VERSION)
#endif

#endif //VERSION_H_
