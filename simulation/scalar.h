#ifndef SCALAR_H_
#define SCALAR_H_

#if defined(USE_FLOAT)
typedef float Scalar;
#define Scalar_PRECISION (FLT_DIG+2)
#else
typedef double Scalar;
#define Scalar_PRECISION (DBL_DIG+2)
#endif //USE_FLOAT

#endif //SCALAR_H_
