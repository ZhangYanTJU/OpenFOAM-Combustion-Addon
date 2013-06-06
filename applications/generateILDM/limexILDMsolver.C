#include "ILDM.H"
#include "mixtureFraction.H"
#include "limexILDMsolver.H"
#include "addToRunTimeSelectionTable.H"



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



namespace Foam
{

defineTypeNameAndDebug(limexILDMsolver, 0);
addToRunTimeSelectionTable(ILDMsolver, limexILDMsolver, dictionary);


limexILDMsolver::limexILDMsolver(ILDM& ildm)
    : ILDMsolver(ildm)
{
}

limexILDMsolver::~limexILDMsolver()
{
}

    
Cantera::compositionMap limexILDMsolver::solve
(
    double pv,
    const Cantera::compositionMap& initial
)
{
    pvval_=pv;
    csol=this;

    int n=ildm_.gas().nSpecies();
    double y[n], ys[n];
    toArray(ildm_.gas(), initial, y);
    memset(ys, 0, sizeof(ys));

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
    t_end   = 1.0;
    rtol    = 1e-9;
    atol    = 1e-9;
    h       = 1e-9;

    limex_(
        &n,
        &objectiveF_LIMEXadaptor,
        NULL,
        &t_begin,
        &t_end,
        y,
        ys,
        &rtol,
        &atol,
        &h,
        Iopt,
        Ropt,
        IPos,
        IFail
    );

    cout<<"IFAIL="<<IFail[0]<<" "<<IFail[1]<<" "<<IFail[2]<<endl;
 
    if (IFail[0]) 
    {
        throw ILDMsolverError("LIMEX failed!");
    }
    else
    {
        return toMap(ildm_.gas(), y);
    }
    
}

ColumnVector limexILDMsolver::objectiveF(const ColumnVector& Y, Matrix& B)
{
    //std::cout<<"objectiveF"<<std::endl;

    Cantera::compositionMap mY=
        toMap(ildm_.gas(), Y);
    
    ColumnVector DYDt;
    Matrix J=
        F
        (
            ildm_.gas(), 
            mY,
            ildm_.h(), 
            ildm_.p(),
            &DYDt
        );

    //std::cout<<"DYDt: "<<std::endl<<DYDt<<std::endl;
    //std::cout<<"Jacobian: "<<std::endl<<J<<std::endl;

    Matrix Zdach;

    double ratio;
    ColumnVector result(ildm_.gas().nSpecies(), 0.0);
    result.insert
        (
            ildm_.P() * Y - ildm_.tau(),
            0
        );
    result.insert
        (
            calcSchur(ildm_.gas(), J, DYDt, ildm_.nSlow(), ratio, &Zdach), 
            ildm_.nSlow()
        );

    B=Matrix(ildm_.gas().nSpecies(), ildm_.gas().nSpecies(), 0.0);
    B.insert(Zdach, ildm_.nSlow(), 0);

    //std::cout<<"nleq_F="<<result<<std::endl;

    //std::cout<<"end objectiveF"<<std::endl;
    return result;

}


limexILDMsolver* limexILDMsolver::csol=NULL;

void limexILDMsolver::objectiveF_LIMEXadaptor
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
    if (csol!=NULL)
    {
        // copy into appropriate data structures
        ColumnVector Y(*n, 0.0);
        for (int i=0;i<*n;i++)
            Y(i)=y[i]; 
 
        // computed system
        Matrix B;
        ColumnVector rhs=csol->objectiveF(Y, B);

        // copy back into fortran types
        for (int i=0;i<*n;i++)
            f[i]=rhs(i);

        int nzz=0;
        for (int i=0;i<B.rows();i++)
            for (int j=0;j<B.columns();j++)
            {
                if (fabs(B(i,j))>SMALL)
                {
                    nzz++;
                    b[nzz-1]=B(i,j);
                    ir[nzz-1]=i+1;
                    ic[nzz-1]=j+1;
                }
            }
        *nz=nzz;

    }
    else
    {
        {
            using namespace Foam;
            FatalErrorIn
                (
                    "limexILDMsolver::objectiveF_adaptor()"
                )
                << "Invalid function pointer"
                    << exit(FatalError);
        }
    }
}


}
