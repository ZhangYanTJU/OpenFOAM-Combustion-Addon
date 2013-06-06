#ifndef EULEX_H
#define EULEX_H

typedef 
void (*EULEX_FCN)
    (
        int* n,
        double* t,
        double* y,
        double* dy
    );

extern "C" void eulex_
(
    int* n,
    EULEX_FCN,
    double* t,
    double* y,
    double* tend,
    double* tol,
    double* hmax,
    double* h,
    int* kflag
);


template<class ODEIntegrand>
class EULEX
{

    static ODEIntegrand* i_;
    
    static void integrand_adaptor
    (
        int* n,
        double* t,
        double* y,
        double* dy
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

            ColumnVector y=y0;
            int n=y0.length();
            double startTime=0;

            double tol=1e-6, hmax=1e-3, h=1e-9;
            int kflag=1;
            
            eulex_
                (
                    &n,
                    &integrand_adaptor,
                    &startTime,
                    y.fortran_vec(),
                    &endTime,
                    &tol,
                    &hmax,
                    &h,
                    &kflag
                );
            
            if (kflag<0) 
            {
                std::string errmsg;
                switch (kflag)
                {
                    case -1:
                        errmsg="endTime smaller than startTime";
                        break;
                    case -2:
                        errmsg="more than NSTMAX basic integration steps per interval have been performed";
                        break;
                    case -3:
                        errmsg="more than JRMAX stepsize reductions occurred per basic tntegration step";
                        break;
                    case -4:
                        errmsg="stepsize proposal for next basic integration too small";
                        break;
                    default:
                        errmsg="unknown error";
                }
                throw errmsg;
            }
            else
            {
                return y;
            }

        }

};

template<class ODEIntegrand>
ODEIntegrand* EULEX<ODEIntegrand>::i_=NULL;

template<class ODEIntegrand>
void EULEX<ODEIntegrand>::integrand_adaptor
(
    int* n,
    double* t,
    double* y,
    double* dy
)
{
    if (i_!=NULL)
    {
        ColumnVector vy(*n, 0.0);
        for (int i=0; i<*n; i++) vy(i)=y[i];
        
        ColumnVector dydt = i_->integrand(vy, *t);
        
        for (int i=0; i<*n; i++) dy[i]=dydt(i);
    }
    else
    {
        std::cout<<"Error in EULEX::integrand_adaptor():"<<std::endl
            <<"Pointer to ODEIntegrator object not set!"<<std::endl;
    }

}



#endif
