#ifndef LIMEX_H
#define LIMEX_H

typedef void (*LIMEX_FCN)
  (
   int *n,
   int *nz,
   double *t,
   double *y,
   double *f,
   double *b,
   int *ir,
   int *ic,
   int *FcInfo
);



extern "C" void limex_
(
 int* intn,
 LIMEX_FCN,
 void *,
 double* t_begin,
 double* t_End,
 double* y,
 double* ys,
 double* rTol,
 double* aTol,
 double* h,
 int *Iopt,
 double* Ropt,
 int* IPos,
 int* IFail
 );


template<class ODEIntegrand>
class LIMEX
{

    static ODEIntegrand* i_;
    
    static void integrand_adaptor
    (
        int *n,
        int *nz,
        double *t,
        double *y,
        double *f,
        double *b,
        int *ir,
        int *ic,
        int *FcInfo
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
            ColumnVector y=y0;
            ColumnVector ys(n, 0.0);

            int Iopt[30], IPos[n], IFail[3];
            double Ropt[5];
            Iopt[0]=1;
            Iopt[1]=0;
            Iopt[2]=1;
            Iopt[3]=0;
            Iopt[4]=0;
            Iopt[5]=0;//1;
            Iopt[6]=0;
            Iopt[7]=0;
            Iopt[8]=0;
            Iopt[9]=0;
            Iopt[10]=0;
            Iopt[11]=0;
            Iopt[12]=0;
            Iopt[13]=0;
            Iopt[14]=0;
            Iopt[15]=0;
            Iopt[16]=0;
            Iopt[17]=0;
            Iopt[18]=0;
            Iopt[19]=0;
            Iopt[20]=0;
            Iopt[21]=0;
            Iopt[22]=0;
            Iopt[23]=0;
            Iopt[24]=0;
            Iopt[25]=0;
            Iopt[26]=0;
            Iopt[27]=0;
            Iopt[28]=0;
            Iopt[29]=0;

            Ropt[0]=0;
            Ropt[1]=0;
            Ropt[2]=0;
            Ropt[3]=0;
            Ropt[4]=0;

            for (int i=0;i<n;i++)
                IPos[i]=0;

            double t_begin, t_end, rtol, atol, h;

            t_begin = 0.0;
            t_end   = endTime;
            rtol    = 1e-10;
            atol    = 1e-12;
            h       = 1e-9;
            
            limex_(
                &n,
                &integrand_adaptor,
                NULL,
                &t_begin,
                &t_end,
                y.fortran_vec(),
                ys.fortran_vec(),
                &rtol,
                &atol,
                &h,
                Iopt,
                Ropt,
                IPos,
                IFail
            );

 
            if (IFail[0]) 
            {
                cout<<"IFAIL="<<IFail[0]<<" "<<IFail[1]<<" "<<IFail[2]<<endl;
                throw std::string("LIMEX failed!");
            }
            else
            {
                return y;
            }

        }

};

template<class ODEIntegrand>
ODEIntegrand* LIMEX<ODEIntegrand>::i_=NULL;

template<class ODEIntegrand>
void LIMEX<ODEIntegrand>::integrand_adaptor
(
    int *n,
    int *nz,
    double *t,
    double *y,
    double *f,
    double *b,
    int *ir,
    int *ic,
    int *FcInfo
)
{
    if (i_!=NULL)
    {

        // copy into ColumnVector
        ColumnVector vy(*n, 0.0);
        for (int i=0; i<*n; i++) vy(i)=y[i];
        
        // call integrand
        Matrix B;
        ColumnVector rhs;
        i_->integrand
            (
                vy, *t,
                B, rhs
            );
        
        // copy back
        for (int i=0; i<*n; i++) f[i]=rhs(i);

        if (B.rows()==0 && B.cols()==0)
        {
            int nzz=0;
            // unit matrix
            for (int i=0;i<*n;i++)
            {
                nzz++;
                b[nzz-1]=1.0;
                ir[nzz-1]=i+1;
                ic[nzz-1]=i+1;
            }
            *nz=*n;
        }
        else
        {
            int nzz=0;
            for (int i=0;i<B.rows();i++)
                for (int j=0;j<B.columns();j++)
                {
                    if (fabs(B(i,j))>1e-15)
                    {
                        nzz++;
                        b[nzz-1]=B(i,j);
                        ir[nzz-1]=i+1;
                        ic[nzz-1]=j+1;
                    }
                }
            *nz=nzz;
        }

    }
    else
    {
        std::cout<<"Error in LIMEX::integrand_adaptor():"<<std::endl
            <<"Pointer to ODEIntegrator object not set!"<<std::endl;
    }

}




#endif
