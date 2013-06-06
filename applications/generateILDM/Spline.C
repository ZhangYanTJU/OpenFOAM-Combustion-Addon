#include "Spline.H"

using namespace std;


Spline::Spline
(
    const std::map<double, ColumnVector>& data
)
{
    std::map<double, const double*> data2;

    int spnum=data.begin()->second.length();;

    for (std::map<double, ColumnVector>::const_iterator it=data.begin();
         it!=data.end(); it++)
    {
        data2[it->first]=it->second.data();
    }

    create(spnum, data2);
}

Spline::Spline
(
    int spnum,
    const std::map<double, const double*>& data
)
{
    create(spnum, data);
}

void Spline::create
(
    int spnum,
    const std::map<double, const double*>& data
)
{
    int np=data.size();
    n_=np;
    x_.reset(new double[np]);
    int j=0;
    for (std::map<double, const double*>::const_iterator it=data.begin();
         it!=data.end();it++)
    {
        x_.get()[j++]=it->first;
    }

    for (int i=0; i<spnum; i++)
    {       
        accels.push_back(gsl_interp_accel_alloc ());
        splines.push_back(gsl_spline_alloc (gsl_interp_cspline, np));

        double y[np];
        j=0;
        for (std::map<double, const double*>::const_iterator it=data.begin();
             it!=data.end();it++)
        {
            y[j++]=it->second[i];
        }

        gsl_spline_init(splines[i], x_.get(), y, np);
    }
}

Spline::~Spline()
{
    for (int i=0; i<splines.size(); i++)
    {
        gsl_spline_free(splines[i]);
        gsl_interp_accel_free(accels[i]);
    }
}

ColumnVector Spline::interpolate(double x)
{
    x=min( (&(*x_))[n_-1], max((&(*x_))[0], x));

    ColumnVector y(splines.size(), 0.0);

    for (int i=0; i<splines.size(); i++)
    {
        y(i)=gsl_spline_eval(splines[i], x, accels[i]);
    }

    return y;
}


ColumnVector Spline::firstDerivative(double x)
{
    ColumnVector dydt(splines.size(), 0.0);

    for (int i=0; i<splines.size(); i++)
    {
        dydt(i)=gsl_spline_eval_deriv(splines[i], x, accels[i]);
    }

    return dydt;
}

ColumnVector Spline::secondDerivative(double x)
{
    ColumnVector dydt(splines.size(), 0.0);

    for (int i=0; i<splines.size(); i++)
    {
        dydt(i)=gsl_spline_eval_deriv2(splines[i], x, accels[i]);
    }

    return dydt;
}
