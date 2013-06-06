#ifndef DASSL2_H
#define DASSL2_H

#include "octave/DASSL.h"

template<class ODEIntegrand>
class DASSL2
{

    static ODEIntegrand* i_;
    static int n_;
    
    static ColumnVector integrand_adaptor
    (
        const ColumnVector& y,
        const ColumnVector& ydot,
        double t,
        int &ires
    );

public:

    ColumnVector integrate
    (
        ODEIntegrand* integrand,
        const ColumnVector& y0,
        double endTime
    )
        {
            i_=integrand;

            DAEFunc func(&integrand_adaptor);
            DASSL dassl(y0, 0, func);
            return dassl.do_integrate(endTime);
            /*
            if (idid<0) 
            {
                cout<<"IDID="<<idid<<endl;
                throw std::string("DASSL failed!");
            }
            else
            {
                return y;
            }
                */
        }

};

template<class ODEIntegrand>
ODEIntegrand* DASSL2<ODEIntegrand>::i_=NULL;

template<class ODEIntegrand>
int DASSL2<ODEIntegrand>::n_=0;


double lastt=0.0;

template<class ODEIntegrand>
ColumnVector DASSL2<ODEIntegrand>::integrand_adaptor
(
    const ColumnVector& y,
    const ColumnVector& ydot,
    double t,
    int &ires
)
{
    if (i_!=NULL)
    {

        if (t!=lastt)
            cout<<"t="<<t<<endl;
        lastt=t;
        // call integrand
        return ydot - i_->integrand(y, t);
    }
    else
    {
        std::cout<<"Error in DASSL::integrand_adaptor():"<<std::endl
            <<"Pointer to ODEIntegrator object not set!"<<std::endl;
        return ColumnVector();
    }
}

#endif
