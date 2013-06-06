
#include "ChemicalSystem.H"
#include "octave/dbleSVD.h"

#include "mixtureFraction.H"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include <cantera/Cantera.h>
#include <cantera/IdealGasMix.h>
#include <cantera/onedim.h>
#include <cantera/equilibrium.h>
#include <cantera/transport.h>

#define MINTEMP 290.0

void fromMolesToMolefractions(Cantera::compositionMap& mix)
{
  double n_ges=0;
  for (Cantera::compositionMap::const_iterator it=mix.begin();
       it!=mix.end(); it++)
    n_ges+=it->second;
  for (Cantera::compositionMap::iterator it=mix.begin();
       it!=mix.end(); it++)
    it->second/=n_ges;
}

void fromMolefractionsToMassfractions
(const Cantera::IdealGasMix& gas, Cantera::compositionMap& mix)
{

  double denom=0;

  for (Cantera::compositionMap::const_iterator it=mix.begin();
       it!=mix.end(); it++)
    denom += it->second * gas.molarMass(gas.speciesIndex(it->first));

  for (Cantera::compositionMap::iterator it=mix.begin();
       it!=mix.end(); it++)
    it->second *= gas.molarMass(gas.speciesIndex(it->first))/denom;
}

/*
Foam::chemicalSystemState ChemicalSystem::unburntState
(
    const Foam::HashTable<Foam::label,Foam::word>& desc,
    double z, double p, double h
) const
{
    setGasTo(z, p, h);
    return FGM1D::convert(*this, desc);
}
*/

ColumnVector ChemicalSystem::unburntComposition
(
    double z, double p, double h, bool convert
) const
{
    setGasTo(z, p, h, convert);
    ColumnVector Y(gas_.nSpecies(), 0.0);
    gas_.getMassFractions(const_cast<double *>(Y.data()));
    return Y;
}


ColumnVector ChemicalSystem::equilibriumComposition
(
    double z, double p, double h, bool convert
) const
{
    setGasTo(z, p, h, convert);
    equilibrate(gas_, "HP");
    ColumnVector Y(gas_.nSpecies(), 0.0);
    gas_.getMassFractions(const_cast<double *>(Y.data()));
    return Y;
}


Matrix nullspace(Matrix A)
{
  SVD svd(A);
  Matrix U=svd.left_singular_matrix();
  DiagMatrix S=svd.singular_values();
  Matrix V=svd.right_singular_matrix();

  int rows=A.rows();
  int cols=A.cols();

  int S_nr=S.rows();
  int S_nc=S.cols();

  ColumnVector s;
  if (S_nr == 1 || S_nc == 1)
    s = S.row(1);
  else
    s = S.diag();

  double tol = A.cols() * s (1) * 1e-3;

  int rank=0;
  for (int i=0;i<s.length();i++)
    if (s(i)>tol) rank+=1;

  Matrix retval;
  if (rank < cols)
    retval = V.extract(0,rank,V.rows()-1,cols-1);
  else
    retval = Matrix(cols, cols, 0.0);

  return retval;
}




ChemicalSystem::ChemicalSystem
(
    const std::string& fname,
    const std::string& id,
    const Cantera::compositionMap& fuel,
    const Cantera::compositionMap& oxi,
    const Cantera::compositionMap& prod,
    const std::string& ctype,
    const std::string& mfatom,
    double stoichOF
)
    : gas_(fname, id),
      tr_(newTransportMgr("Mix", &gas_)),
      oxi_(oxi),
      fuel_(fuel),
      prod_(prod),
      mixFracAtom_(mfatom),
      stoichOF_(stoichOF)
{
    if (ctype=="moleNumbers")
    {
        fromMolesToMolefractions(fuel_);
        fromMolesToMolefractions(oxi_);
        fromMolesToMolefractions(prod_);
        fromMolefractionsToMassfractions(gas_, fuel_);
        fromMolefractionsToMassfractions(gas_, oxi_);
        fromMolefractionsToMassfractions(gas_, prod_);
    } else if (ctype=="moleFractions")
    {
        fromMolefractionsToMassfractions(gas_, fuel_);
        fromMolefractionsToMassfractions(gas_, oxi_);
        fromMolefractionsToMassfractions(gas_, prod_);
    } else if (ctype=="massFractions")
    {
        // already as we want it
    }

    cBr_=Z(gas().elementIndex(mixFracAtom_), gas(), fuel_);
    cProd_=Z(gas().elementIndex(mixFracAtom_), gas(), prod_);

}


#define VECTOR_DATA_PTR(x) (x.data())



ColumnVector ChemicalSystem::DYDt
(
    const double* Y,
    double h,
    double p,
    double* rho_copy
) const
{
    //std::cout<<"DYDt"<<std::endl;

    int nsp=gas_.nSpecies();

    gas_.setState_TPY(2300.0, p, Y);
    gas_.setState_HP(h, p);
    doublereal c[nsp];
    gas_.getConcentrations(c);

    doublereal DcDt[nsp];        
    gas_.getNetProductionRates(DcDt); // [kmol/m^3/s]

    double rho=gas_.density();
    if (rho_copy) *rho_copy=rho;

    double DrhoDt=0.0;
    for (int i=0; i<nsp; i++) 
        DrhoDt+=gas_.molarMass(i)*DcDt[i];

 
    ColumnVector res(nsp, 0.0);
    for (int i=0; i<nsp; i++) 
        res(i)=gas_.molarMass(i)*(rho*DcDt[i]-c[i]*DrhoDt)/(rho*rho);

    return res;
}

Matrix ChemicalSystem::F
(
    const double* Yorg,
    double h,
    double p,
    ColumnVector* DYDt_copy,
    double* rho_copy
) const
{
    int nsp=gas_.nSpecies();

    ColumnVector dydt_org=DYDt(Yorg, h, p, rho_copy);
    if (DYDt_copy) *DYDt_copy=dydt_org;

    Matrix res(nsp, nsp, 0.0);
    for (int i=0; i<nsp; i++)
    {
        double Y[nsp];
        memcpy(Y, Yorg, sizeof(Y));
        Y[i]=Y[i]+max(1e-5*Y[i], 1e-8);

        ColumnVector DdydtDy=(DYDt(Y, h, p)-dydt_org)/(Y[i]-Yorg[i]);
        res.insert(DdydtDy, 0, i);
    }

    return res;
}



extern "C" void dtrsyl_
(
 char *,
 char *,
 int *,
 int *, int *,
 double *, int *,
 double *, int *,
 double *, int *,
 double *,
 int *
);



Matrix givens(double a, double b, int n=2, int i=0)
{
  double c,s;
  double tau;
  if (b==0)
    {
      c=1;
      s=0;
    }
  else if (fabs(b)>fabs(a))
    {
      tau=-a/b; s=1/(sqrt(1+tau*tau)); c=s*tau;
    }
  else
    {
      tau=-b/a; c=1/(sqrt(1+tau*tau)); s=c*tau;
    }
  Matrix ret(n,n,0.0);

  if (n>2)
    for (int j=0;j<n;j++)
      ret(j,j)=1.0;

  ret(i,i)=c;
  ret(i,i+1)=s;
  ret(i+1,i)=-s;
  ret(i+1,i+1)=c;

  return ret;
}







void interchange(Matrix& T, Matrix& Q, int k)
{
  int n=T.rows();
  Matrix gmat=givens(T(k,k+1), T(k+1,k+1)-T(k,k)/*, n, k*/);
  
  T.insert(gmat.transpose()*T.extract(k,k,k+1,n-1), k,k);
  T.insert(T.extract(0,k,k+1,k+1)*gmat, 0,k);
  Q.insert(Q.extract(0,k,n-1,k+1)*gmat, 0,k);
  
  /*
  T=gmat.transpose()*T*gmat;
  Q=Q*gmat;
  */
}





// Q - Schur vectors
// T - Triangular matrix
void sortSchur(Matrix& T, Matrix& Q, int num)
{
  for (int i=0; i<T.rows()-1; i++)
    {
      for (int j=0; j < T.rows()-1-i; j++)
	{
	  if ( T(j,j) < T(j+1,j+1) )
	    {
	      interchange(T,Q,j);
	    }
	}
    }

}


void ChemicalSystem::subSpace
(
    const Matrix& J,
    const ColumnVector& DYDt,
    int nSlow,
    Matrix* Zdach_fast,
    Matrix* Zslow,
    Matrix* Zdach_slow,
    Matrix* Zfast,
    Matrix* D
) const
{

  int nFast=gas_.nSpecies()-nSlow;

  // compute schur decomposition
  SCHUR ss(J, "u");

  Matrix T=ss.schur_matrix();
  Matrix Q=ss.unitary_matrix();
    
  sortSchur(T, Q, nSlow);
  //std::cout<<"Q:"<<std::endl<<Q<<std::endl;
  if (D) *D=T;

  Matrix Qt=Q.transpose();

  if (Zdach_fast) *Zdach_fast = Qt.extract_n(nSlow, 0, nFast, nSlow+nFast);
  if (Zslow) *Zslow = Q.extract_n(0, 0, nSlow+nFast, nSlow);

  if (Zdach_slow!=NULL || Zfast!=NULL)
  {
      Matrix Ns=
          T.extract_n
          (
              0, 0, 
              nSlow, nSlow
          );

      Matrix Nf=
          T.extract_n
          (
              nSlow, nSlow, 
              nFast, nFast
          );

      Matrix Nfs=
          T.extract_n
          (
              0, nSlow, 
              nSlow, nFast
          );

      Matrix X=-Nfs;

      char n='N'; int p=-1;
      int a_nr=Ns.rows();
      int b_nr=Nf.rows();
      double scale=1.0;
      int info;

      dtrsyl_
          (
              &n, &n, &p,
              &a_nr, &b_nr, 
              Ns.fortran_vec(), &a_nr, 
              Nf.fortran_vec(), &b_nr,
              X.fortran_vec(), &a_nr,
              &scale, &info
          );
      
      if (info<0) 
      {
          cout<<"Error in DTRSYL"<<endl;
          exit(-2);
      }
      if(info>0) cout<<"DTRSYL INFO="<<info<<endl;
    
      
      Matrix Y(nFast, gas_.nSpecies(), 0.0);
      for (int i=0;i<Y.rows(); i++)
          Y(i,i)=1.0;
      Y.insert(-X, 0, nSlow);
      //std::cout<<Y*Qt<<std::endl;
      
      Matrix Zdach=Y*Qt;
      //std::cout<<Zdach.rows()<<" "<<Zdach.cols()<<std::endl;

      if (Zdach_slow) *Zdach_slow = Zdach.extract_n(0, 0, nSlow, nSlow+nFast);

      if (Zfast)
      {
          Y.insert(X, 0, nSlow);
          Zdach=Q*Y;
          *Zfast=Zdach.extract_n(0, nSlow, nSlow+nFast, nFast);
      }
  }

  /*    
  // find orthogonal complement of P
  Matrix mP=P(15);
  Matrix orth=nullspace(mP);

  Matrix Zds=Zdach.extract_n(0,0,nSlow,nState);
  Matrix Zdf=Zdach.extract_n(nSlow,0,nFast,nState);

  Matrix Zds_star=(Zds*mP.transpose()).inverse()*Zds;
  Matrix xx=Zdf*orth;
  cout<<xx.rows()<<" "<<xx.cols()<<endl;
  Matrix Zdf_star=(Zdf*orth).inverse()*Zdf;

  Zdach.insert(Zds_star, 0, 0);
  Zdach.insert(Zdf_star, nSlow, 0);
  */

  /*    
  Y.insert(X, 0, nSlow);

  Matrix Z=Q*Y;
      */
  /*
  if (dFdtCopy!=NULL)
    *dFdtCopy=dFdt;
  if (ZinvCopy!=NULL)
    *ZinvCopy=Zdach;
  if (Zcopy!=NULL)
    *Zcopy=Z;
  if (Jcopy!=NULL)
    *Jcopy=J;
      */
  
  //std::cout<<Zdach*DYDt<<std::endl;

  /*
  double resi=0.0;
  for (int i=0;i<sol.length();i++)
  {
      std::cout<<"res["<<i<<"]="<<sol(i)<<std::endl;
      resi+=sol(i)*sol(i);
  }
  std::cout<<"Res="<<std::sqrt(resi)<<std::endl;
      */

  //std::cout<<"end calcSchur"<<std::endl;
}



ColumnVector ChemicalSystem::calcSchur
(
    const Matrix& J,
    const ColumnVector& DYDt,
    int nSlow,
    Matrix* Zdach_copy
) const
{

  int nFast=gas_.nSpecies()-nSlow;

  // compute schur decomposition
  SCHUR ss(J, "u");

  Matrix T=ss.schur_matrix();
  Matrix Q=ss.unitary_matrix();
    
  sortSchur(T, Q, nSlow);
  //std::cout<<"Q:"<<std::endl<<Q<<std::endl;

  Matrix Qt=Q.transpose();
  /*      
  Matrix Ns=
    T.extract_n
    (
     0, 0, 
     nSlow, nSlow
     );

  Matrix Nf=
    T.extract_n
    (
     nSlow, nSlow, 
     nFast, nFast
     );

  Matrix Nfs=
    T.extract_n
    (
     0, nSlow, 
     nSlow, nFast
     );

  Matrix X=-Nfs;

  char n='N'; int p=-1;
  int a_nr=Ns.rows();
  int b_nr=Nf.rows();
  double scale=1.0;
  int info;
  dtrsyl_
    (
     &n, &n, &p,
     &a_nr, &b_nr, 
     Ns.fortran_vec(), &a_nr, 
     Nf.fortran_vec(), &b_nr,
     X.fortran_vec(), &a_nr,
     &scale, &info
     );

  if (info<0) 
    {
      cout<<"Error in DTRSYL"<<endl;
      exit(-2);
    }
  if(info>0) cout<<"DTRSYL INFO="<<info<<endl;
    
      
  Matrix Y(nFast, gas.nSpecies(), 0.0);
  for (int i=0;i<Y.rows(); i++)
    Y(i,i)=1.0;
  Y.insert(-X, 0, nSlow);
  //std::cout<<Y*Qt<<std::endl;
      
  Matrix Zdach=Y*Qt;
  std::cout<<Zdach.rows()<<" "<<Zdach.cols()<<std::endl;
  */    

  Matrix Zdach=Qt.extract_n(nSlow, 0, nFast, nSlow+nFast);

  if (Zdach_copy!=NULL) *Zdach_copy=Zdach;

  /*    
  // find orthogonal complement of P
  Matrix mP=P(15);
  Matrix orth=nullspace(mP);

  Matrix Zds=Zdach.extract_n(0,0,nSlow,nState);
  Matrix Zdf=Zdach.extract_n(nSlow,0,nFast,nState);

  Matrix Zds_star=(Zds*mP.transpose()).inverse()*Zds;
  Matrix xx=Zdf*orth;
  cout<<xx.rows()<<" "<<xx.cols()<<endl;
  Matrix Zdf_star=(Zdf*orth).inverse()*Zdf;

  Zdach.insert(Zds_star, 0, 0);
  Zdach.insert(Zdf_star, nSlow, 0);
  */

  /*    
  Y.insert(X, 0, nSlow);

  Matrix Z=Q*Y;
      */
  /*
  if (dFdtCopy!=NULL)
    *dFdtCopy=dFdt;
  if (ZinvCopy!=NULL)
    *ZinvCopy=Zdach;
  if (Zcopy!=NULL)
    *Zcopy=Z;
  if (Jcopy!=NULL)
    *Jcopy=J;
      */
  
  //std::cout<<Zdach*DYDt<<std::endl;

  ColumnVector sol=Zdach*DYDt;


  return sol;
}




struct enthalpyByProdParams
{
    const ChemicalSystem* chem;
    double z, p, h;
};

double enthalpyByYprod(double Yprod, void* params)
{
    enthalpyByProdParams* p=(enthalpyByProdParams*) params;
    return p->chem->setGasToMixFracWithProducts(p->z, p->p, Yprod, MINTEMP) - p->h;
}




double ChemicalSystem::setGasToMixFracWithProducts(double z, double p, double Yprod, double temp) const
{
    //cout<<"z="<<z()<<" cProd="<<cProd<<" cBr="<<cBr<<endl;
    double Yl=1.0 - Yprod*(cBr_ - cProd_)/cBr_ - z;
    double Ybr=1.0 - Yl - Yprod;

    //cout<<"Ybr="<<Ybr<<" Yprod="<<Yprod<<" Yl="<<Yl<<endl;

    compositionMap Y=
        Ybr * fuel_
        +
        Yprod * prod_
        +
        Yl * oxi_;
    
    gas().setState_TPY(temp, p, Y);    
    //cout<<"RES=("<<gas_->enthalpy_mass()<<"-"<<h()<<")="<<gas_->enthalpy_mass() - h()<<endl;
    return gas().enthalpy_mass();
}

void ChemicalSystem::setGasToTemp(double z, double p, double T) const
{
    Cantera::compositionMap Y=
        z * fuel_
        +
        (1.0-z) * oxi_;

    gas().setState_TPY(T, p, Y);
}


void ChemicalSystem::setGasTo(const double* Y, double p, double h) const
{
    gas().setState_TPY(300.0, p, Y);
cout<<gas().cp_mass()<<std::endl;
    gas().setState_HP(h, p);
}

void ChemicalSystem::setGasTo(double z, double p, double h, bool convert) const
{

    if (!convert)
    {
        setGasToTemp(z, p, 1.0);
        if (h<=gas().enthalpy_mass())
        {
            using namespace Foam;
            FatalErrorIn
                (
                    "ChemicalSystem::setGasTo()"
                )   << "Specified enthalpy h=" << h 
                    << " is below the absolute minimum enthalpy hzero=" << gas().enthalpy_mass() << endl
                    << exit(FatalError);
        }
        gas().setState_HP(h, p);
    }
    else
    {
        if ( h > hmin(z, p) )
        {
            //cout<<z<<" "<<p<<" "<<h<<endl;
            setGasToTemp(z, p, 1.0);
            gas().setState_HP(h, p);
        }
        else if ( h>=hzero(z, p) )
        {
            // convert part of unburnt mix into products to decrease enthalpy at constant temperature
        
            int status;
            int iter = 0, max_iter = 100;
            const gsl_root_fsolver_type *T;
            gsl_root_fsolver *s;
            double r = 0;
        
            double x_lo = 0.0, 
                x_hi = z/(1.0-((cBr_-cProd_)/cBr_));
        
            gsl_function F;
        
            //cout<<"setting function"<<endl;
            enthalpyByProdParams ps;
            ps.chem=this;
            ps.p=p;
            ps.z=z;
            ps.h=h;
            F.function = &enthalpyByYprod;
            F.params = (void *) &ps;

            T = gsl_root_fsolver_brent;
            s = gsl_root_fsolver_alloc (T);
            //cout<<"setup gsl"<<endl;
            gsl_root_fsolver_set (s, &F, x_lo, x_hi);
           /* 
                printf ("using %s method\n", 
                gsl_root_fsolver_name (s));
            
                printf ("%5s [%9s, %9s] %9s %10s %9s\n",
                "iter", "lower", "upper", "root", 
                "err", "err(est)");
            */   
            do
            {
                iter++;
                status = gsl_root_fsolver_iterate (s);
                r = gsl_root_fsolver_root (s);
                x_lo = gsl_root_fsolver_x_lower (s);
                x_hi = gsl_root_fsolver_x_upper (s);
                status = gsl_root_test_interval (x_lo, x_hi,
                0, 0.001);
              /*  
                    if (status == GSL_SUCCESS)
                    printf ("Converged:\n");
                
                    printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
                    iter, x_lo, x_hi,
                    r, 
                    x_hi - x_lo);
                */    
            }
            while (status == GSL_CONTINUE && iter < max_iter);
        
            gsl_root_fsolver_free (s);
        
            setGasToMixFracWithProducts(z, p, r, MINTEMP);
        
        } else
        {
            using namespace Foam;
            FatalErrorIn
                (
                    "ChemicalSystem::setGasTo()"
                )   << "Specified enthalpy h=" << h 
                    << " is below the absolute minimum enthalpy hzero=" << hzero(z, p) << endl
                    << exit(FatalError);
        }
    }
}

double ChemicalSystem::hmin(double z, double p) const
{
    setGasToTemp(z, p, MINTEMP);
    double h=gas().enthalpy_mass();
    std::cout<<"hmin="<<h<<std::endl;
    return h;
}

double ChemicalSystem::hzero(double z, double p) const
{
    setGasToTemp(z, p, MINTEMP);
    double h=setGasToMixFracWithProducts(z, p, z/(1.0-((cBr_-cProd_)/cBr_)), MINTEMP);
    std::cout<<"hzero="<<h<<std::endl;
    return h;
}
