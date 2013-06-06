#include "ILDM.H"
#include "mixtureFraction.H"
#include "nleq1ILDMsolver.H"
#include "addToRunTimeSelectionTable.H"
#include "NLEQ1.H"

namespace Foam
{

defineTypeNameAndDebug(nleq1ILDMsolver, 0);
addToRunTimeSelectionTable(ILDMsolver, nleq1ILDMsolver, dictionary);


nleq1ILDMsolver::nleq1ILDMsolver(ILDM& ildm)
    : ILDMsolver(ildm)
{
}

nleq1ILDMsolver::~nleq1ILDMsolver()
{
}

    
ColumnVector nleq1ILDMsolver::solve
(
    double pv,
    const ColumnVector& initial
)
{
    pvval_=pv;

    NLEQ1<nleq1ILDMsolver> solver;

    try
    {
        return solver.solve(this, initial);
    }
    catch (std::string errmsg)
    {
        throw ILDMsolverError(errmsg);
    }
    catch (CanteraError)
    {
        showErrors(cout);
        throw ILDMsolverError("Error in Cantera");
    }
    
}

ColumnVector nleq1ILDMsolver::objective(const ColumnVector& Y)
{
    ColumnVector DYDt;
    Matrix J=
        ildm_.sys().F
        (
            Y.data(),
            ildm_.h(), 
            ildm_.p(),
            &DYDt
        );

    ColumnVector result(ildm_.gas().nSpecies(), 0.0);
    result.insert
        (
            ildm_.P() * Y - ildm_.tau(),
            0
        );
    result.insert
        (
            ildm_.sys().calcSchur(J, DYDt, ildm_.nSlow()), 
            ildm_.nSlow()
        );

    return result;

}

Matrix nleq1ILDMsolver::jacobian(const ColumnVector& Y)
{
    ColumnVector DYDt;
    Matrix J=
        ildm_.sys().F
        (
            Y.data(),
            ildm_.h(), 
            ildm_.p(),
            &DYDt
        );
    Matrix Zdach;
    ildm_.sys().calcSchur(J, DYDt, ildm_.nSlow(), &Zdach);

    Matrix Jac(ildm_.gas().nSpecies(), ildm_.gas().nSpecies(), 0.0);
    Jac.insert(ildm_.P(), 0, 0);
    Jac.insert(Zdach*J, ildm_.nSlow(), 0);
 
    return Jac;

}



}
