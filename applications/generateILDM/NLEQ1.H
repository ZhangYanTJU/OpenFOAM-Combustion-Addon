#ifndef NLEQ1_H
#define NLEQ1_H

typedef void (*NLEQ_FCN)
  (
   int *n,
   double *x,
   double *f,
   int *ifail
   );



typedef void (*NLEQ_JAC)
  (
   int*,
   int*,
   double*,
   double*,
   int*
  );


extern "C" void nleq1_
(
 int*,
 NLEQ_FCN,
 NLEQ_JAC,
 double*,
 double*,
 double*,
 int*,
 int*,
 int*,
 int*,
 int*,
 double*
);


template<class NLEQFunc>
class NLEQ1
{

    static NLEQFunc* f_;

    static void objective_adaptor
    (
        int *n,
        double *x,
        double *f,
        int *ifail
    );

    static void jacobian_adaptor
    (
        int* n,
        int* ldjac,
        double* x,
        double* jac,
        int* ifail
    );

public:

    ColumnVector solve(NLEQFunc* f, const ColumnVector& xstart, bool jacobian=true)
        {
            f_=f;
                        
            ColumnVector x=xstart;
            int n=xstart.length();

            ColumnVector xscal(n, 0.0);
            double rtol;
            int iopt[50], ierr=0;
            int liwk, lrwk;
            liwk=n+50;
            lrwk=(n+n+13)*n+61;
            int iwk[liwk];
            double rwk[lrwk];
            
            for (int i=0;i<50;i++)
                iopt[i]=0;
            for (int i=0;i<liwk;i++)
                iwk[i]=0;
            for (int i=0;i<lrwk;i++)
                rwk[i]=0.0;
          
            if (jacobian)
                iopt[2]=1; // user supplied Jacobian
            else
                iopt[2]=0;

            iopt[3]=0; // dense jacobian matrix
            
            //output
            
            iopt[10]=1;
            iopt[12]=1;
            iopt[14]=0;
            iopt[18]=0;
                
            //output units
            /*
                iopt[11]=9;
                iopt[13]=9;
                iopt[15]=2;
                iopt[19]=9;
                */
            
            iopt[30]=4; // extremly nonlinear
            iopt[31]=1; // allow broyden updates
            iopt[32]=1; //ordinary newton iteration, no damping
            //iopt[33]=1; //simplified newton iteration
            iopt[37]=1;
            
            iwk[30]=500; // max number of iterations
            
            //rwk[21]=1e-10; // min damping factor
            //rwk[19]=1.1;
            
            rtol=1e-5;
            if (jacobian)
                nleq1_
                    (
                        &n,
                        &objective_adaptor,
                        &jacobian_adaptor,
                        x.fortran_vec(),
                        xscal.fortran_vec(),
                        &rtol,
                        iopt,
                        &ierr,
                        &liwk,
                        iwk,
                        &lrwk,
                        rwk
                    );
            else
                nleq1_
                    (
                        &n,
                        &objective_adaptor,
                        NULL,
                        x.fortran_vec(),
                        xscal.fortran_vec(),
                        &rtol,
                        iopt,
                        &ierr,
                        &liwk,
                        iwk,
                        &lrwk,
                        rwk
                    );

            
            cout<<"RTOL="<<rtol<<std::endl;            
            
            if (ierr>0) 
            {
                std::ostringstream errmsg;
                errmsg<<"NLEQ1 return error code "<<ierr;
                throw errmsg.str();
            }
            else
            {
                return x;
            }
        }
    
};

template<class NLEQFunc>
NLEQFunc* NLEQ1<NLEQFunc>::f_;

template<class NLEQFunc>
void NLEQ1<NLEQFunc>::objective_adaptor
(
    int *n,
    double *x,
    double *f,
    int *ifail
)
{
    if (f_!=NULL)
    {
        // copy into appropriate data structures
        ColumnVector Y(*n, 0.0);
        for (int i=0;i<*n;i++)
            Y(i)=x[i]; 
 
        // compute system
        ColumnVector rhs;
        try
        {
            rhs=f_->objective(Y);
        }
        catch (int errcode)
        {
            *ifail=errcode;
        }

        // copy back into fortran types
        for (int i=0;i<rhs.length();i++)
            f[i]=rhs(i);
    }
    else
    {
        std::cout<<"Error in NLEQ1::objective_adaptor():"<<std::endl
            <<"Pointer to NLEQFunc object not set!"<<std::endl;
        *ifail=-1;
    }

}

template<class NLEQFunc>
void NLEQ1<NLEQFunc>::jacobian_adaptor
(
    int* n,
    int* ldjac,
    double* x,
    double* jac,
    int* ifail
)
{

    if (f_!=NULL)
    {
        ColumnVector Y(*n, 0.0);
        for (int i=0;i<*n;i++)
            Y(i)=x[i];
            
        try
        {
            Matrix G=f_->jacobian(Y);
                
            if ((*ldjac!=G.rows())&&(*n!=G.cols()))
            {
                std::cout<<"Error in NLEQ1::jacobian_adaptor():"<<std::endl
                    <<"Jacobian size does not match!"<<std::endl;
                *ifail=-2;
            } 
            else
            {            
                for (int i=0;i<G.rows();i++)
                    for (int j=0;j<G.cols();j++)
                        jac[j*(*ldjac)+i]=G(i,j);
            }
        }
        catch (int errcode)
        {
            *ifail=errcode;
        }
    }
    else
    {
        std::cout<<"Error in NLEQ1::jacobian_adaptor():"<<std::endl
            <<"Pointer to NLEQFunc object not set!"<<std::endl;
        *ifail=-1;
    }
     
}
 


#endif
