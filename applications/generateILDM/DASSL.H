#ifndef DASSL_H
#define DASSL_H

typedef void (*DASSL_FCN)
  (
      double *t,
      double *y,
      double *dydt,
      double *delta,
      int *ires,
      double *rpar,
      int *ipar
);



extern "C" void ddassl_
(
 DASSL_FCN,
 int *n,
 double *time,
 double* y,
 double* dydt,
 double* t_End,
 int *info,
 double* rTol,
 double* aTol,
 int *idid,
 double* rwork,
 int *lrwk,
 int *iwrk,
 int *liwrk,
 double* rpar,
 int* ipar,
 void *jac
 );


template<class ODEIntegrand>
class DASSL
{

    static ODEIntegrand* i_;
    static int n_;
    
    static void integrand_adaptor
    (
      double *t,
      double *y,
      double *dydt,
      double *delta,
      int *ires,
      double *rpar,
      int *ipar
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

            int n=y0.length();
            n_=n;

            ColumnVector y=y0;
            ColumnVector dydt=i_->integrand(y, 0);

            double t=0.0, tout=endTime;

            int JBW=3;

            int liw=20+n;
            int iwrk[liw];
            cout<<liw<<endl;
            memset(iwrk, 0, sizeof(iwrk));

            //int lrw=40+(5+4)*n+n*n;
            int lrw=40+(5+4)*n+(2*JBW+JBW+1)*n+2*(n/(JBW+JBW+1)+1);
            cout<<lrw<<endl;
            double rwrk[lrw];
            memset(rwrk, 0, sizeof(rwrk));

            int info[15];
            memset(info, 0, sizeof(info));

            info[0]=0;
            info[1]=0;
            info[2]=0;
            info[3]=0;
            info[4]=0;
            info[5]=1; // banded jacobian
            iwrk[0]=JBW; // lower bandwidth
            iwrk[1]=JBW; // upper bandwidth
            info[6]=0;
            info[7]=0;
            info[8]=0;
            info[9]=0;
            info[10]=0; // consistent start values set

            double rtol    = 1e-9;
            double atol    = 1e-9;

            double rpar=0.0;
            int ipar=0;

            int idid=0;
            ddassl_
                (
                    &integrand_adaptor,
                    &n,
                    &t,
                    y.fortran_vec(),
                    dydt.fortran_vec(),
                    &tout,
                    info,
                    &rtol,
                    &atol,
                    &idid,
                    rwrk,
                    &lrw,
                    iwrk,
                    &liw,
                    &rpar,
                    &ipar,
                    NULL
                );
            
 
            if (idid<0) 
            {
                cout<<"IDID="<<idid<<endl;
                throw std::string("DASSL failed!");
            }
            else
            {
                return y;
            }

        }

};

template<class ODEIntegrand>
ODEIntegrand* DASSL<ODEIntegrand>::i_=NULL;

template<class ODEIntegrand>
int DASSL<ODEIntegrand>::n_=0;

template<class ODEIntegrand>
void DASSL<ODEIntegrand>::integrand_adaptor
(
    double *t,
    double *y,
    double *dydt,
    double *delta,
    int *ires,
    double *rpar,
    int *ipar
)
{
    if (i_!=NULL)
    {

            cout<<"start integ"<<endl;
        // copy into ColumnVector
        ColumnVector vy(n_, 0.0);
        for (int i=0; i<n_; i++) vy(i)=y[i];
        
        // call integrand
        ColumnVector rhs=i_->integrand(vy, *t);
        
        // copy back
        for (int i=0; i<n_; i++) delta[i]=rhs(i)-dydt[i];
    }
    else
    {
        std::cout<<"Error in DASSL::integrand_adaptor():"<<std::endl
            <<"Pointer to ODEIntegrator object not set!"<<std::endl;
    }
    return 0;
}

#endif
