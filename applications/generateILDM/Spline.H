#ifndef SPLINE_H
#define SPLINE_H

#include <vector>
#include <map>
#include <memory>

#include "octave/config.h"
#include "octave/dColVector.h"

#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"

class Spline
{

    int n_;
    std::auto_ptr<double> x_;
    std::vector<gsl_interp_accel*> accels;
    std::vector<gsl_spline*> splines;
    
    void create
    (
        int spnum,
        const std::map<double, const double*>& data
    );

public:
    Spline
    (
        const std::map<double, ColumnVector>& data
    );

    Spline
    (
        int spnum,
        const std::map<double, const double*>& data
    );
    ~Spline();

    ColumnVector interpolate(double x);
    ColumnVector firstDerivative(double x);
    ColumnVector secondDerivative(double x);

    inline double x(int i) const { return x_.get()[i]; }
  
};

#endif
