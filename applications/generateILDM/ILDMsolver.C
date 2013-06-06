
#include "ILDMsolver.H"
#include "ILDM.H"
#include "autoPtr.H"

namespace Foam
{

    
ILDMsolverError::ILDMsolverError(const std::string& m)
    : errorMessage_(m)
{
}

defineTypeNameAndDebug(ILDMsolver, 0);
defineRunTimeSelectionTable(ILDMsolver, dictionary);

autoPtr<ILDMsolver> ILDMsolver::New
(
    ILDM& ildm,
    const word& solverName
)
{
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(solverName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "ILDMsolver::New()"
        )   << "Unknown ILDMsolver type " << solverName
            << endl << endl
            << "Valid ILDMsolver types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return Foam::autoPtr<ILDMsolver>(cstrIter()(ildm));
}

ILDMsolver::ILDMsolver(ILDM& ildm)
    : ildm_(ildm)
{
}

ILDMsolver::~ILDMsolver()
{
}

}
