
#include "ILDM.H"
#include "mixtureFraction.H"
#include "cantera/equilibrium.h"

#include "ILDMsolver.H"

#include <fstream>
#include <sstream>

#include "EULEX.H"
#include "NLEQ1.H"

#include "addToRunTimeSelectionTable.H"


namespace Foam
{

defineTypeNameAndDebug(ILDM, 0);
addToRunTimeSelectionTable(LDM, ILDM, dictionary);


class ILDMpredictorError
: public std::string
{
public:
    ILDMpredictorError(const std::string& msg)
        : std::string(msg)
        {
        }
};


ILDM::ILDM
(
    const ChemicalSystem& s,
    const std::string& pv,
    const std::string& solverName,
    double z, double p, double h,
    int resolution
)
    : LDM(s, pv, z, p, h, resolution),
      solverName_(solverName)
{
}

ILDM::ILDM
(
    const ChemicalSystem& s, 
    double z,
    double p,
    double h,
    dictionary& dict
)
    : LDM(s, z, p, h, dict),
      solverName_(word(dict.lookup("ILDMsolver")))
{
}




// integrand of predictor
ColumnVector ILDM::integrand
(
    const ColumnVector& Y,
    double
) const
{
    ColumnVector DYDt;

    Matrix J=sys().F(Y.data(), h(), p(), &DYDt);

    Matrix Zslow;
    sys().subSpace(J, DYDt, nSlow(), NULL, &Zslow);

    Matrix ZP = ( P() * Zslow );
    Matrix iZP = ZP.inverse();
    if (ZP==iZP) 
    {
        throw std::string("Predictor: inversion failed");
    }

    ColumnVector d(nSlow(), 0.0);
    d(nSlow()-1)=-1;

    return Zslow*iZP*d;
}


ColumnVector ILDM::predict
(
    const ColumnVector& Ystart,
    double step
)
{
    try
    {
        EULEX<ILDM> integrator;
        
        return integrator.integrate(this, Ystart, step);
    }
    catch (std::string msg)
    {
        throw ILDMpredictorError(msg);
    }
}

bool ILDM::onILDM
(
    const ColumnVector& Y
)
{
    ColumnVector DYDt;
    Matrix T;
    
    Matrix J=sys().F(Y.data(), h(), p(), &DYDt);
    sys().subSpace(J, DYDt, nSlow(), NULL, NULL, NULL, NULL, &T);
 
    return 
        (T(nSlow_-1, nSlow_-1) / T(nSlow_, nSlow_)) 
        < 0.2;
}


void ILDM::calculate()
{

    Cantera::compositionMap Yu=sys_.Y(z_);

    gas().setState_TPY(300.0, p_, Yu);

    try 
    {
        equilibrate(gas(), "HP");
    }
    catch (CanteraError) 
    {
        showErrors(cout);
        {
            using namespace Foam;
            FatalErrorIn
            (
                "ILDM::fill()"
            )
                << "Could not solve for equilibrium composition."
                   << "Debug output from cantera displayed above"
                      << exit(FatalError);
        }
    }
    std::cout<<"adiabatic flame T="<<gas().temperature()<<std::endl;

    ColumnVector Y(gas().nSpecies(), 0.0);
    gas().getMassFractions(Y.fortran_vec()); // Y: equilibrium composition

#warning !!FIX: Assumes maximum value of progress variable at equilibrium
    double stepsize = Yp(Y) / double(resolution_);

    tau_ = P_ * Y;

    double pvval=Y(gas().speciesIndex(pv_));
    append(pvval, Y);


    // == STEP 1: calculate ILDM from equilibrium point ==

    Foam::autoPtr<Foam::ILDMsolver> solver=
        Foam::ILDMsolver::New(*this, solverName_);

    pvval-=stepsize;
    for (; pvval>0; pvval-=stepsize)
    {

        std::cout<<"Solving for PV="<<pvval<<std::endl;
        tau_(nSlow_-1)=pvval;

        try
        {


            Y=solver().solve
                (
                    pvval,
                    predict(Y, stepsize)
                );

            if ( ! onILDM(Y) ) break;

            append(pvval, Y);

        }

        catch (Foam::ILDMsolverError e)
        {
            std::cout<<"No solution possible: STOP"<<std::endl;
            std::cout<<"Message: "<<e.errorMessage()<<std::endl;
            break;
        }

        catch (ILDMpredictorError e)
        {
            std::cout<<"No prediction possible: STOP"<<std::endl;
            std::cout<<"Message: "<<e<<std::endl;
            break;
        }

    }

    // == STEP 2: extend from ILDM boundary to unburnt point ==
    Ybound=begin()->second;
    Y=Ybound;

    Yunburnt=toVector(gas(), Yu);
    Porth=nullspace(P_).transpose();

    pvval=begin()->first;
    double pbound=pvval;

    for (pvval-=stepsize; pvval>0; pvval-=stepsize)
    {

        std::cout<<"Solving for PV="<<pvval<<std::endl;
        tau_(nSlow_-1)=pvval;

        try
        {
            // linear interpolation
            append(pvval, ((pbound-pvval)*Yunburnt+pvval*Ybound)/pbound);
        }
        catch (std::string errmsg)
        {
            cout<<"ERROR: STOP"<<endl;
            cout<<errmsg<<endl;
            break;
        }
        catch (CanteraError)
        {
            showErrors(cout);
            break;
        }
    }
}


Matrix ILDM::jacobian(const ColumnVector& Y) const
{
    return Matrix();
}


void ILDM::outputEigenvalueSpectrum(std::ofstream& f)
{
    for (const_iterator it=begin();
         it!=end(); it++)
    {
        cout<<"Analyzing PV="<<it->first<<endl;

        try
        {
            const ColumnVector& Y=it->second;
            
            ColumnVector Ydot;
            Matrix D;

            Matrix J=sys().F(Y.data(), h_, p_, &Ydot);
            sys().subSpace(J, Ydot, nSlow_, NULL, NULL, NULL, NULL, &D);
         
            f<<it->first;
            for (int m=0;m<nSlow_+5;m++)
                f<<" "<<D(m,m);
            f<<std::endl;
        }
        catch (CanteraError)
        {
            showErrors(cout);
        }
    }
}

    /*
void ILDM::addToTable
(
    chemistryTable& tab,
    const LDMconservedVars& cv
)
{
    // get table index to current parameter set          
    int idx[tab.dimensionality()];
    for (int dim=0; dim<tab.dimensionality()-1; dim++)
    {
        idx[dim]=
            tab.indexAt
            (
                dim, 
                cv.getByName
                (
                    tab.nameOfPV(dim)
                )
            );
    }

    Spline ldm_interpolator(*this);

    // add 1D LDM to table
    label pvdim=tab.dimensionality()-1;
    for (
        idx[pvdim]=0;
        idx[pvdim]<tab.nElements(pvdim);
        idx[pvdim]++
    )
    {

        scalar pv= tab.valueAt(pvdim, idx[pvdim]) * pveq();
        Info<<"adding to table: pv="<<pv<<endl;
        
        tab.access(idx)=
            extendToTableEntry
            (
                tab,
                ldm_interpolator.interpolate(pv)
            );

    }

}
        */

    /*
void ILDM::appendTableDescriptions
(
    ListDescription& pvdesc,
    ListDescription& contentdesc
)
{
    pvdesc.appendEntry(pv_);

    contentdesc.appendEntry("Ydot"&pv_);

    LDM::appendTableDescriptions(pvdesc, contentdesc);
}
        */

    /*
void ILDM::setupTable
(
    ListDescription& pvdesc,
    ListDescription& contentdesc,
    label* num,
    double* start,
    double* delta
)
{
    label i=pvdesc[pv_];
    num[i]=resolution_;
    start[i]=0.0;
    delta[i]=1.0/scalar(resolution_-1);
}
        */

    /*
chemicalSystemState ILDM::extendToTableEntry
(
    const chemistryTable& tab,
    const ColumnVector& Y
)
{

    const ListDescription& desc=tab.contentDescription();              

    chemicalSystemState s=LDM::extendToTableEntry(tab, Y);

    int nsp=gas().nSpecies();

    doublereal DcDt[nsp], c[nsp];
    gas().getNetProductionRates(DcDt);
    gas().getConcentrations(c);

    double rho=gas().density();

    double DrhoDt=0.0;
    for (int i=0; i<nsp; i++) 
        DrhoDt+=gas().molarMass(i)*DcDt[i];

    s[desc["Ydot"&pv_]]=gas().molarMass(pvidx_)*
        (rho*DcDt[pvidx_]-c[pvidx_]*DrhoDt)/(rho*rho);

    return s;
}
        */

}
