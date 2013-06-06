
#include "ILDM.H"
#include "mixtureFraction.H"
#include "octaveILDMsolver.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(octaveILDMsolver, 0);
addToRunTimeSelectionTable(ILDMsolver, octaveILDMsolver, dictionary);


octaveILDMsolver::octaveILDMsolver(ILDM& ildm)
    : ILDMsolver(ildm)
{
}

octaveILDMsolver::~octaveILDMsolver()
{
}

    
Cantera::compositionMap octaveILDMsolver::solve
(
    double pv,
    const Cantera::compositionMap& initial
)
{
    pvval_=pv;
    csol=this;

    NLFunc func( &objectiveF_adaptor );
    NLEqn solver
        (
            toVector(ildm_.gas(), initial), 
            func
        );
    
    int state=0;
    ColumnVector cY=solver.solve(state);

    if (state==1)
    {
        return toMap(ildm_.gas(), cY);
    } else
    {
        throw ILDMsolverError(solver.error_message());
    }
    
}

ColumnVector octaveILDMsolver::objectiveF(const ColumnVector& Y)
{
    std::cout<<"objectiveF"<<std::endl;

    Cantera::compositionMap mY=
        toMap(ildm_.gas(), Y, ildm_.pv());
    
    mY[ildm_.pv()]=pvval_;
    
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

    std::cout<<"DYDt: "<<std::endl<<DYDt<<std::endl;

    //std::cout<<"Jacobian: "<<std::endl<<J<<std::endl;

    double ratio;
    ColumnVector result(ildm_.gas().nSpecies(), 0.0);
    result.insert
        (
            ildm_.P() * Y - ildm_.tau(),
            0
        );
    result.insert
        (
            calcSchur(ildm_.gas(), J, DYDt, ildm_.nSlow(), ratio), 
            ildm_.nSlow()
        );

    std::cout<<"nleq_F="<<result<<std::endl;

    std::cout<<"end objectiveF"<<std::endl;
    return result;

}


octaveILDMsolver* octaveILDMsolver::csol=NULL;

ColumnVector octaveILDMsolver::objectiveF_adaptor(const ColumnVector& Y)
{
    if (csol!=NULL)
        return csol->objectiveF(Y);
    else
    {
        {
            using namespace Foam;
            FatalErrorIn
                (
                    "octaveILDMsolver::objectiveF_adaptor()"
                )
                << "Invalid function pointer"
                    << exit(FatalError);
        }
        return ColumnVector();
    }
}


}
