
#include "REDIM.H"
#include "LIMEX.H"
//#include "DASSL2.H"
//#include "EULEX.H"
#include "Spline.H"

#include "octave/dRowVector.h"

#include "addToRunTimeSelectionTable.H"

double norm(const ColumnVector& c)
{
    double norm=0.0;
    for (int i=0; i<c.length(); i++)
        norm+=c(i)*c(i);
    return sqrt(norm);
}


/**
    Implements implicit smoother for RHS of partial differential equations
    to overcome time step size restrictions in explicit integration methods.

    F. W. Wubs, "Stabilization of Explicit Methods for Hyperbolic Partial
    Differential Equations", International Journal for Numerical Methods in
    Fluids, Vol. 6, pp 641-657 (1986)
    */

class RHSsmoother
{
    Matrix Ci;

public:

    RHSsmoother
    (
        int n,
        double deltaT,
        double minDeltaX,
        double sbC=1.0
    )
        {
            
            double mu=1.1* 0.25 * deltaT*deltaT / (sbC*sbC * minDeltaX*minDeltaX);
            cout<<"dT="<<deltaT<<" dx="<<minDeltaX<<" mu="<<mu<<endl;

            Matrix C(n, n, 0.0);
            for (int i=1; i<n-1; i++)
            {
                C(i, i+1) = -mu;
                C(i, i  ) = 1.0 + 2.0*mu;
                C(i, i-1) = -mu;
            }

            C(  0,   0) = 0.5;
            C(  0,   1) = 0.5;
            C(n-1, n-1) = 0.5;
            C(n-1, n-2) = 0.5;

            Ci = C.inverse();

            //cout<<"Ci="<<endl<<Ci<<endl;
        }

    inline ColumnVector operator()(const ColumnVector& RHS) const
        {
            return Ci * RHS;
        }
};
    

namespace Foam
{

defineTypeNameAndDebug(REDIM, 0);
addToRunTimeSelectionTable(LDM, REDIM, dictionary);


void REDIM::calculate()
{

    auto_ptr<std::vector<ColumnVector> > redim
        (
            new std::vector<ColumnVector> 
        );

    double mindx=1e3;
    for (const_iterator it=(++begin()); it!=(--end()); it++)
    {
        redim->push_back(it->second);

        // evaluate CFL criterion
        const_iterator prev=it, next=it;
        prev--;
        next++;
        mindx = min(mindx, 0.5*(next->first - prev->first));
    }

    //RHSsmoother smooth(nsp_, h_, mindx);

    double ptime=0;
    try
    {
        ofstream res("residual.t");
        for (int iter=0; iter<iterations_; iter++)
        {
            ptime+=h_;

            redim=splineInterpolate(*redim);

            /*
            if (iter%100==0)
            {
                ostringstream n;
                n << "rd_"<<iter<<".out";
                ofstream f(n.str().c_str());
                for (int k=0;k<redim->size();k++)
                    f<<(*redim)[k].transpose()<<endl;
            }
                */

            std::vector<Matrix> Psi_theta, Psi_theta2;
            splineDerivatives(*redim, Psi_theta, Psi_theta2);
            
            double residual=0.0;
            for (int i=0; i<redim->size(); i++)
            {
                Matrix Proj(nsp_, nsp_, 0.0);
                for (int k=0; k<Proj.rows(); k++)
                    Proj(k,k) = 1.0;
                
                Proj -= Psi_theta[i] * Psi_theta[i].pseudo_inverse();

                double rho;
                ColumnVector F;
                Matrix F_psi=sys().F((*redim)[i].data(), h_, p_, &F, &rho);

                ColumnVector pt=Psi_theta[i].column(0);
                ColumnVector G = (constant_/rho)* Psi_theta2[i].column(0)/(pt.transpose()*pt * double(nsp_));

                Matrix L(nsp_, nsp_, 0.0);
                for (int k=0; k<L.rows(); k++)
                    L(k,k) = 1.0;
                
                L -= h_ * Proj * F_psi;

                ColumnVector delta = L.inverse() * h_ * Proj * ( F - G );
                residual += norm(delta/h_);
                (*redim)[i] += delta;
            }

            res<<ptime<<" "<<residual<<endl;
            cout<<"REDIM-Iter #"<<iter<<": T="<<ptime<<"\t res="<<residual<<endl;
            if ( residual < 1) break;

            /*
            for (int i=0; i<redim->size(); i++)
            {
                cout<< (P_* (*redim)[i]).transpose() <<endl;
            }
                */
        }
        
    }
    catch (CanteraError)
    {
            showErrors(cout);
    }

    int j=0;
    for (iterator it=(++begin()); it!=(--end()); it++)
        it->second=(*redim)[j++];

}


REDIM::REDIM
(
    const ChemicalSystem& s, 
    const std::string& pv,
    const std::string& solverName,
    double z, double p, double h,
    int resolution,
    double c
)
    : ILDM
      (
          s, pv, 
          solverName, 
          z, p, h, 
          resolution
      ),
      nsp_(gas().nSpecies()),
      constant_(c),
      h_(1e-6),
      iterations_(300)
{
}



REDIM::REDIM
(
    const ChemicalSystem& s, 
    double z,
    double p,
    double h,
    dictionary& dict
)
    : ILDM
      (
          s, 
          z, p, h, 
          dict
      ),
      nsp_(gas().nSpecies()),
      constant_(readScalar(dict.lookup("diffusionParameter"))),
      h_(readScalar(dict.lookup("timestep"))),
      iterations_(readLabel(dict.lookup("iterations")))
{
}


REDIM::~REDIM()
{
}


std::auto_ptr<Spline> REDIM::spline
(
    const std::vector<ColumnVector>& Y
) const
{
    std::map<double, const double*> spp;

    const ColumnVector& Yu=begin()->second;
    spp[Yp(Yu)]=Yu.data();

    const ColumnVector& Ye=(--end())->second;
    spp[Yp(Ye)]=Ye.data();

    for (std::vector<ColumnVector>::const_iterator it=Y.begin(); 
         it!=Y.end(); it++)
    {
        spp[Yp(*it)]=it->data();
    }

    return auto_ptr<Spline>(new Spline(nsp_, spp));
}

auto_ptr<std::vector<ColumnVector> > REDIM::splineInterpolate
(
    const std::vector<ColumnVector>& Y
) const
{

    auto_ptr<Spline> sp=spline(Y);

    auto_ptr<std::vector<ColumnVector> > yi_ptr
        (
            new std::vector<ColumnVector>
        );
    std::vector<ColumnVector>& yi = *yi_ptr;

    for (const_iterator it=(++begin()); it!=(--end()); it++)
    {
        yi.push_back(sp->interpolate(it->first));
    }

    return yi_ptr;
}


void REDIM::splineDerivatives
(
    const std::vector<ColumnVector>& Y,
    std::vector<Matrix>& DPsiDtheta,
    std::vector<Matrix>& DPsiDtheta2
) const
{
    DPsiDtheta.clear();
    DPsiDtheta2.clear();

    auto_ptr<Spline> sp=spline(Y);

    for (int i=0; i<Y.size(); i++)
    {
        DPsiDtheta.push_back(Matrix(sp->firstDerivative(sp->x(i))));
        DPsiDtheta2.push_back(Matrix(sp->secondDerivative(sp->x(i))));
    }
    
}

}
