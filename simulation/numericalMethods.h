#pragma once
#ifndef NUMERICALMETHODS_H
#define NUMERICALMETHODS_H
#include <valarray>
#include <vector>
#include <assert.h>
#include <stdarg.h>
#include <stdexcept>
#include <limits>
using namespace std;

namespace ODESolvers {
typedef valarray<double> ODEState;
typedef valarray<double> Derivative;

//typedef void (*DerivativeFunc)(const ODEState& i, double t, double dt, const Derivative& d, Derivative& out);

struct DerivativeFunc {
  virtual void operator()(const ODEState& /*s*/, double /*t*/, Derivative& /*out*/, void* /*params*/) const {}
  virtual size_t dim() const { return 0; }  //it's own null object
  virtual ~DerivativeFunc() {}
};

//*** TESTING METHODS ***///
#if 0
// Testing Function based on spring constant
// x'' = -k*x - b*x'
struct SpringDeriv : DerivativeFunc {
  //Mapped struct to map valarray to names
  /*virtual*/ void operator()(const ODEState& s, double /*t*/, Derivative& out, void* /*params*/) const {
    const Scalar K = 0.001, B = 0.01;
    out[0] = s[1]; //The value of dx
    out[1] = -K * (s[0])  /* -K * x, Hooke's Law */
              -B * (s[1]); /* -B * v, Dampening term */
  }
  /*virtual*/ size_t dim() const { return 2; }
#if 0
  //Analytical solution (for testing) for x(0) = 5, x'(0) = 6
  //Without dampening - wolframalpha.com
  double operator()(double t) const {
    return 189.737*sin(0.0316228*t)+5*cos(0.0316228*t);
  }
#else
  //Analytical solution (for testing) for x(0) = 5, x'(0) = 6
  //With dampening - wolframalpha.com
  double operator()(double t) const {
    return exp(-0.005*t)*(192.954*sin(0.031225*t)+5*cos(0.031225*t));
  }
#endif
};
#endif
// Testing Function for dependent variable functions
// u' = -2*u + t + 4
struct EqnDeriv : DerivativeFunc {
  //Mapped struct to map valarray to names
  /*virtual*/ void operator()(const ODEState& s, double t, Derivative& out, void* /*params*/) const {

    out[0] = -2 * s[0] + t + 4;
  }
  /*virtual*/ size_t dim() const { return 1; }

  //Analytical solution (for testing)
  double operator()(double t) const {
    // For u(0) = 1
    return 0.25*(2*t - 3*exp(-2*t) + 7);
  }
};

///*** END TESTING METHODS ***///

// ************
// * Steppers *
// ************

//Stepper interface
struct Stepper {
  Stepper(size_t /*dim*/) {}
  virtual void step(ODEState& i, DerivativeFunc& fn, double t, double dt, double& lasttimestep, Derivative& error, void* params=0) = 0;
  virtual size_t order() const = 0;
  virtual Stepper* clone() const = 0;
  virtual ~Stepper() {}
};

//Runge-Kutta 3rd order. \ref http://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Kutta.27s_third-order_method
struct RK2Stepper : Stepper {
  Derivative a, b;
  RK2Stepper(size_t dim)
    : Stepper(dim), 
     a(dim),
     b(dim) {}

  /*virtual*/ void step(ODEState& i, DerivativeFunc& fn, double t, double dt, double& /*lasttimestep*/, Derivative& /*error*/, void* params=0) {
    fn(i, t, a, params);
    a*=dt;
    fn(i+0.5*a, t+dt*0.5, b, params);
    b*=dt;
    i += b;
  }
  /*virtual*/ size_t order() const { return 2; }
  /*virtual*/ Stepper* clone() const { return new RK2Stepper(*this); }
};

//Runge-Kutta 3rd order. \ref http://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Kutta.27s_third-order_method
struct RK3Stepper : Stepper {
  Derivative a, b, c;
  RK3Stepper(size_t dim)
    : Stepper(dim), 
     a(dim),
     b(dim),
     c(dim) {}

  /*virtual*/ void step(ODEState& i, DerivativeFunc& fn, double t, double dt, double& /*lasttimestep*/, Derivative& /*error*/, void* params=0) {
    fn(i, t, a, params);
    a*=dt;
    fn(i+0.5*a, t+dt*0.5, b, params);
    b*=dt;
    fn(i+2.0*b-a, t+dt, c, params);
    c*=dt;
    i += (1.0 / 6.0) * (a + 4.0 * b + c);
  }
  /*virtual*/ size_t order() const { return 3; }
  /*virtual*/ Stepper* clone() const { return new RK3Stepper(*this); }
};

//Runge-Kutta 4th order. \ref http://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Classic_fourth-order_method
struct RK4Stepper : Stepper {
  Derivative a, b, c, d, tmp;
  RK4Stepper(size_t dim)
    : Stepper(dim), 
     a(dim),
     b(dim),
     c(dim),
     d(dim),
     tmp(dim) {}

  /*virtual*/ void step(ODEState& i, DerivativeFunc& fn, double t, double dt, double& /*lasttimestep*/, Derivative& /*error*/, void* params=0) {
    fn(i, t+0.0, a, params);
    a*=dt;
    tmp=i+0.5*a;
    fn(tmp, t+dt*0.5, b, params);
    b*=dt;
    tmp=i+0.5*b;
    fn(tmp, t+dt*0.5, c, params);
    c*=dt;
    tmp=i+c;
    fn(tmp, t+dt, d, params);
    d*=dt;
    i += (1.0 / 6.0) * (a + 2.0 * (b + c) + d);

//    for (size_t j=0; j<fn.dim();j++) {
//        std::cout << i[j] << ",";
//    }

//    std::cout << dt << std::endl;

  }
  /*virtual*/ size_t order() const { return 4; }
  /*virtual*/ Stepper* clone() const { return new RK4Stepper(*this); }
};

//Forward-Euler. \ref http://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Forward_Euler
struct ForwardEulerStepper : Stepper {
  Derivative a;
  ForwardEulerStepper(size_t dim)
    : Stepper(dim), 
     a(dim) {}

  /*virtual*/ void step(ODEState& i, DerivativeFunc& fn, double t, double dt, double& /*lasttimestep*/, Derivative& /*error*/, void* params=0) {
    fn(i, t+0.0, a, params);
    a *= dt;  //Multiply itself here to prevent temporary variable creation
    i += a;
  }
  /*virtual*/ size_t order() const { return 1; }
  /*virtual*/ Stepper* clone() const { return new ForwardEulerStepper(*this); }
};

//Runge-Kutta 2nd order (x=1). \ref http://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Generic_second-order_method
struct EulerPCStepper : Stepper {
  Derivative a, b;
  EulerPCStepper(size_t dim)
    : Stepper(dim),
     a(dim),
     b(dim) {}
  /*virtual*/ void step(ODEState& i, DerivativeFunc& fn, double t, double dt, double& /*lasttimestep*/, Derivative& /*error*/, void* params=0) {
    fn(i, t, a, params);
    a*=dt;
    fn(i+a, t, b, params);
    b*=dt;
    i += 0.5 * (a + b);
  }
  /*virtual*/ size_t order() const { return 2; }
  /*virtual*/ Stepper* clone() const { return new EulerPCStepper(*this); }
};

//Heun-Euler. http://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Heun.E2.80.93Euler
struct HEStepper : Stepper {
  Derivative a, b, tmp;
  HEStepper(size_t dim)
    : Stepper(dim),
     a(dim),
     b(dim),
     tmp(dim) {}
  /*virtual*/ void step(ODEState& i, DerivativeFunc& fn, double t, double dt, double& /*lasttimestep*/, Derivative& error, void* params=0) {
    fn(i, t, a, params);
    a*=dt;
    tmp=i+a;
    fn(tmp, t+dt, b, params);
    b*=dt;
    i += ((1.0/2.0) * a + (1.0/2.0)*b );
    /* 1st order - 2nd order */
    error = ((1.0/2.0 - 1.0) * a + (1.0/2.0)*b );
  }
  /*virtual*/ size_t order() const { return 2; }
  /*virtual*/ Stepper* clone() const { return new HEStepper(*this); }
};

//Runge-Kutta Cash-Karp. \ref http://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Cash-Karp
struct RKCKStepper : Stepper {
  Derivative a, b, c, d, e, f, tmp;
  RKCKStepper(size_t dim)
    : Stepper(dim),
     a(dim),
     b(dim),
     c(dim),
     d(dim),
     e(dim),
     f(dim),
     tmp(dim) {}
  /*virtual*/ void step(ODEState& i, DerivativeFunc& fn, double t, double dt, double& /*lasttimestep*/, Derivative& error, void* params=0) {
    fn(i, t+0.0, a, params);
    a*=dt;
    tmp = i+(1.0/5.0)*a;
    fn(tmp, t+dt*(1.0/5.0), b, params);
    b*=dt;
    tmp = i+(3.0/40.0)*a + (9.0/40.0)*b;
    fn(tmp, t+dt*(3.0/10.0), c, params);
    c*=dt;
    tmp = i+(3.0/10.0)*a + (-9.0/10.0)*b + (6.0 / 5.0)*c;
    fn(tmp, t+dt*(3.0/5.0), d, params);
    d*=dt;
    tmp = i+(-11.0/54.0)*a + (5.0/2.0)*b + (-70.0 / 27.0)*c + (35.0 / 27.0)*d;
    fn(tmp, t+dt, e, params);
    e*=dt;
    tmp = i+(1631.0/55296.0)*a + (175.0/512.0)*b + (575.0 / 13824.0)*c + (44275.0 / 110592.0)*d + (253.0/4096.0)*e;
    fn(tmp, t+dt*(7.0/8.0), f, params);
    f*=dt;
    i += ((37.0/378.0) * a /*+ 0*b */ + (250.0/621.0) * c + (125.0/594.0)*d /*+ 0*e */ + (512.0/1771.0)*f);

    /* 4th order - 5th order */
    error = ((37.0/378.0 - 2825.0 / 27648.0) * a /*+ 0*b */ + (250.0/621.0 - 18575.0 / 48384.0) * c + (125.0/594.0 - 13525.0 / 55296.0)*d + (-277.0/14336.0)*e  + (512.0/1771.0 - 0.25)*f);
  }
  /*virtual*/ size_t order() const { return 4; }
  /*virtual*/ Stepper* clone() const { return new RKCKStepper(*this); }
};

//Runge-Kutta Fehlberg. \ref http://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Fehlberg
struct RKFStepper : Stepper {
  Derivative a, b, c, d, e, f;
  RKFStepper(size_t dim)
    : Stepper(dim),
     a(dim),
     b(dim),
     c(dim),
     d(dim),
     e(dim),
     f(dim) {}
  /*virtual*/ void step(ODEState& i, DerivativeFunc& fn, double t, double dt, double& /*lasttimestep*/, Derivative& error, void* params=0) {
    fn(i, t+0.0, a, params);
    a*=dt;
    fn(i+0.25*a, t+dt*0.25, b, params);
    b*=dt;
    fn(i+(3.0/32.0)*a + (9.0/32.0)*b,t+dt*(3.0/8.0), c, params);
    c*=dt;
    fn(i+(1932.0/2197.0)*a + (-7200.0/2197.0)*b + (7296.0 / 2197.0)*c, t+dt*(12.0/13.0), d, params);
    d*=dt;
    fn(i+(439.0/216.0)*a + (-8.0)*b + (3680.0 / 513.0)*c + (-845.0 / 4104.0)*d, t+dt, e, params);
    e*=dt;
    fn(i+(-8.0/27.0)*a + (2.0)*b + (-3544.0 / 2565.0)*c + (1859.0 / 4104.0)*d + (-11.0/40.0)*e, t+dt*0.5, f, params);
    f*=dt;
    i += ((25.0/216.0) * a /*+ 0*b */ + (1408.0/2565.0) * c + (2197.0/4104.0)*d + (-1.0/5.0)*e);
    /* 4th order - 5th order */
    error = ((25.0/216.0 - 16.0/135.0) * a /*+ 0*b */ + (1408.0/2565.0 - 6656.0/12825.0) * c + (2197.0/4104.0 - 28561.0/56430.0)*d + (-1.0/5.0 - -9.0/50.0)*e + (0.0 - 2.0/55.0)*f);
  }
  /*virtual*/ size_t order() const { return 4; }
  /*virtual*/ Stepper* clone() const { return new RKFStepper(*this); }
};

//Bogacki-Shampine. \ref http://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Bogacki.E2.80.93Shampine
struct BSStepper : Stepper {
  Derivative a, b, c, d, tmp;
  BSStepper(size_t dim)
    : Stepper(dim),
     a(dim),
     b(dim),
     c(dim),
     d(dim),
     tmp(dim) {}
  /*virtual*/ void step(ODEState& i, DerivativeFunc& fn, double t, double dt, double& /*lasttimestep*/, Derivative& error, void* params=0) {
    fn(i, t+0.0, a, params);
    a*=dt;
    tmp = i+0.5*a;
    fn(tmp, t+dt*0.5, b, params);
    b*=dt;
    tmp = i+(3.0/4.0)*b;
    fn(tmp, t+dt*(3.0/4.0), c, params);
    c*=dt;
    tmp = i+(2.0/9.0)*a + (1.0/3.0)*b + (4.0/9.0)*c;
    fn(tmp, t+dt, d, params);
    d*=dt;
    i += (2.0/9.0) * a + (1.0/3.0)*b  + (4.0/9.0) * c;
    /* 3rd order - 2nd order */
    error = (2.0/9.0 - 7.0/24.0) * a + (1.0/3.0 - 1.0/4.0)*b  + (4.0/9.0 - 1.0/3.0) * c + (0 - 1.0/8.0) * d;
  }
  /*virtual*/ size_t order() const { return 3; }
  /*virtual*/ Stepper* clone() const { return new BSStepper(*this); }
};

#if 1
////TODO: Verify the adaptive stepper
//struct AdaptiveStepper : Stepper {
//  Stepper* stepper;
//  valarray<double> dy, yscal, tmp;
//  double TINY, MAX_STEPS, PSHRINK, PGROW, EPS, SAFETY, ERRCON;
//  //~AdaptiveStepper() { delete stepper; }
//  AdaptiveStepper(size_t dim, Stepper* _s, double _EPS, double _TINY=std::numeric_limits<double>::min(), double _MAX_STEPS=10000, double _PSHRINK=-0.25, double _PGROW=-0.2, double _SAFETY=0.9)
//    : Stepper(dim), stepper(_s), dy(dim), yscal(dim), tmp(dim),
//      TINY(_TINY), MAX_STEPS(_MAX_STEPS), PSHRINK(_PSHRINK), PGROW(_PGROW), EPS(_EPS), SAFETY(_SAFETY), ERRCON(pow(5.0/_SAFETY, 1.0/_PGROW))
//      {}
//  /*virtual*/ void step(ODEState& i, DerivativeFunc& fn, double t, double dt, Derivative& err, void* params=0) {
//    tmp = i;
//    const double t1 = t, t2 = t+dt;
//    for(size_t st=0;st<MAX_STEPS;st++) {
//      fn(i, t, dy, params);

//      for (int j=0; j <= (int)dy.size(); j++)
//      {
//          std::cout << "i: " << i[j] << std::endl;
//      }

//      for (int j=0; j <= (int)dy.size(); j++)
//      {
//          std::cout << "dy: " << dy[j] << std::endl;
//      }

////      yscal = abs(i) + abs(dy)*dt + TINY;
//      yscal = abs(i) + abs(dy)*dt;


//      for (int j=0; j <= (int)yscal.size(); j++)
//      {
//          std::cout << "yscal: " << yscal[j] << std::endl;
//      }

//      if((t+dt - t2) * (t+dt - t1) > 0.0) //Finish up the last timestep
//        dt = t2 - t;
//      for(double h=dt;;) {
//          std::cout << "CurrentStep: " << h << std::endl;
//        stepper->step(tmp, fn, t, h, err, params);

//        for (int j=0; j <= (int)err.size(); j++)
//        {
//            std::cout << "err: " << err[j] << std::endl;
//        }

//        const double errmax = abs(err/yscal).max() / EPS;

//        std::cout << "Max Error: " << errmax << std::endl;

//        std::cin.get();

//        if(errmax > 1.0) {
//          double htmp = SAFETY*h*pow(errmax, PSHRINK);
//          h = (h > 0.0 ? max(htmp, 0.1*h) : min(htmp, 0.1*h));
//          tmp = i;
//        }
//        else {
//          if(errmax > ERRCON) dt = SAFETY*h*pow(errmax, PGROW);
//          else dt = 5.0*h;
//          t += h;
//          i = tmp;
//          std::cout << "NExtStep: " << dt << std::endl;
//          break;
//        }
//      }
//      if((t - t2)*(t2 - t1) >= 0.0)
//        return;
//    }
//    assert(!"Took too many steps");
//    //Warning! - required too many steps, increase max_steps, safety, or decrease pgrow
//  }
//  /*virtual*/ size_t order() const { return stepper->order(); }
//  /*virtual*/ Stepper* clone() const {
//    AdaptiveStepper* ret = new AdaptiveStepper(*this);
//    ret->stepper = stepper->clone();
//    return ret;
//  }

//};


struct AdaptiveStepper : Stepper {
  Stepper* stepper;
  valarray<double> dy, yscal, tmp;
  double TINY, MAX_STEPS, PSHRINK, PGROW, EPS, SAFETY, ERRCON, ABSTOL, RELTOL, MAXSCALE, MINSCALE;
  //~AdaptiveStepper() { delete stepper; }
  AdaptiveStepper(size_t dim, Stepper* _s, double _EPS, double _TINY=std::numeric_limits<double>::min(), double _MAX_STEPS=10000, double _PSHRINK=-0.25,
                  double _PGROW=-0.2, double _SAFETY=0.95, double _ABSTOL=1e-6, double _RELTOL = 1e-5, double _MAXSCALE = 10.0, double _MINSCALE = 0.2)
    : Stepper(dim), stepper(_s), dy(dim), yscal(dim), tmp(dim),
      TINY(_TINY), MAX_STEPS(_MAX_STEPS), PSHRINK(_PSHRINK), PGROW(_PGROW), EPS(_EPS), SAFETY(_SAFETY), ERRCON(pow(5.0/_SAFETY, 1.0/_PGROW)), ABSTOL(_ABSTOL),
      RELTOL(_RELTOL), MAXSCALE(_MAXSCALE), MINSCALE(_MINSCALE)
      {}
  /*virtual*/ void step(ODEState& i, DerivativeFunc& fn, double t, double dt, double& lasttimestep, Derivative& err, void* params=0) {
    tmp = i; // Write initial to tmp

    const double t1 = t, t2 = t+dt; // define start time and end time
    const double maxlasttimestep = (2.0/3.0)*dt; // Don't be too greedy
    const size_t dim = fn.dim(); // Use fn.dim() as little as possible
    double sk,error=0.0, errornorm, scale;
    dt = lasttimestep;

    for(size_t st=0;st<MAX_STEPS;st++) {

        bool reject=false;

        if((t+dt - t2) * (t+dt - t1) > 0.0) //Finish up the last timestep
            dt = t2 - t;

        for(;;) {
            double h=dt; // Set the temperary timestep to its intial value passed in from function call or its last value predicted
            stepper->step(tmp, fn, t, h, lasttimestep, err, params); // Attempt a step

            // Evaluate accuracy of attempted step
            // Ignore absolute tolerance since we have vast differences in magnitude

            for (size_t j=0; j<dim; j++) {
                sk = ABSTOL + (RELTOL*std::max(abs(i[j]+TINY),abs(tmp[j]+TINY)));
                error += (err[j]/sk)*(err[j]/sk);
            }

            // Compute the norm of the error
            errornorm = sqrt(error/dim);

            // Do not use PI error control
            if (errornorm <= 1.0) {
                scale = SAFETY*pow(errornorm,PGROW);
                scale = scale < MINSCALE ? MINSCALE : scale;
                scale = scale > MAXSCALE ? MAXSCALE : scale;

                if (reject)
                    dt = h*std::min(scale,1.0);
                else
                {
                    dt = h*scale;
                    reject=false;
                    t += h;
                    i = tmp;
                    lasttimestep = std::min(dt, maxlasttimestep);
                break;
                }
            }
            else {
                scale = std::max(SAFETY*pow(errornorm,PGROW), MINSCALE);
                dt = h*scale;
                tmp = i;
                reject=true;
            }
        }
        if((t - t2)*(t2 - t1) >= 0.0)
            return;
    }
        assert(!"Took too many steps");
  }
  /*virtual*/ size_t order() const { return stepper->order(); }
  /*virtual*/ Stepper* clone() const {
    AdaptiveStepper* ret = new AdaptiveStepper(*this);
    ret->stepper = stepper->clone();
    return ret;
  }

};



#endif

enum ODEMethod {
  FEuler = 0,
  EulerPC,
  RK2,
  RK3,
  RK4,
  HeunEuler,
  RKCK,
  RKF,
  BogShamp,
  NMethods
};

inline Stepper* StepperFactory(ODEMethod method, size_t dim, bool adaptive=false, double eps=1e-5) {
  //if(adaptive) throw std::runtime_error("Not yet implemented");
  Stepper* ret = NULL;
  switch(method) {
  case FEuler:    ret = new ForwardEulerStepper(dim); break;
  case EulerPC:   ret = new EulerPCStepper     (dim); break;
  case RK2:       ret = new RK2Stepper         (dim); break;
  case RK3:       ret = new RK3Stepper         (dim); break;
  case RK4:       ret = new RK4Stepper         (dim); break;
  case HeunEuler: ret = new HEStepper          (dim); break;
  case RKCK:      ret = new RKCKStepper        (dim); break;
  case RKF:       ret = new RKFStepper         (dim); break;
  case BogShamp:  ret = new BSStepper          (dim); break;
  default:
    throw std::runtime_error("Unknown method");
  }
  if(adaptive)
    ret = new AdaptiveStepper(dim, ret, eps);
  return ret;
}

inline std::ostream& operator<<(std::ostream& s, ODEMethod method) {
  switch(method) {
  case FEuler:    s<<"FEuler"; break;
  case EulerPC:   s<<"EulerPC"; break;
  case RK2:       s<<"RK2"; break;
  case RK3:       s<<"RK3"; break;
  case RK4:       s<<"RK4"; break;
  case HeunEuler: s<<"HeunEuler"; break;
  case RKCK:      s<<"RKCK"; break;
  case RKF:       s<<"RKF"; break;
  case BogShamp:  s<<"BogShamp"; break;
  default:
    throw std::runtime_error("Unknown method");
  }
  return s;
}

} //ODESolvers
#endif //NUMERICALMETHODS_H
