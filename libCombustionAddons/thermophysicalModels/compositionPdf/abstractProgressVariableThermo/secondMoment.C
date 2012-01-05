#include "secondMoment.H"
#include "abstractProgressVariableThermo.H"

#include "FSDWSecondMoment.H"
#include "SIJPDFthermoStateFinder.H"
#include "SIJPDFvarZthermo.H"
#include "SIJPDFvarZthermoStateFinder.H"

#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(secondMomentSolver, 0);
defineRunTimeSelectionTable(secondMomentSolver, dictionary);

typedef FSDWSecondMoment<SIJPDFthermo<SIJPDFthermoIndexDriver> > FSDWconstZSecondMoment;
defineTemplateTypeNameAndDebugWithName(FSDWconstZSecondMoment, "FSD_W<SIJPDFconstZthermo>", 0);
addToRunTimeSelectionTable(secondMomentSolver, FSDWconstZSecondMoment, dictionary);

typedef FSDWSecondMoment<SIJPDFvarZthermo> FSDWvarZSecondMoment;
defineTemplateTypeNameAndDebugWithName(FSDWvarZSecondMoment, "FSD_W<SIJPDFvarZthermo>", 0);
addToRunTimeSelectionTable(secondMomentSolver, FSDWvarZSecondMoment, dictionary);


    autoPtr<secondMomentSolver> secondMomentSolver::New
    (
        const fvMesh& mesh,
        const abstractProgressVariableThermo& thermo,
        const word& variableName,
        const dictionary& dict,
        multivariateSurfaceInterpolationScheme<scalar>::fieldTable& mvtab
    )
    {
        word type(dict.lookup("type"));
        
        Info << "Selecting second moment determination method " << type 
            << " for variable "<< variableName << endl;
        
        dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(type);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorIn
                (
                    "secondMoment::New()"
                )   << "Unknown second moment determination method " << type
                    << endl << endl
                    << "Valid methods are :" << endl
                    << dictionaryConstructorTablePtr_->toc()
                    << exit(FatalError);
        }
        
        return autoPtr<secondMomentSolver>(cstrIter()(mesh, thermo, variableName, dict, mvtab));
        
    }


    secondMomentSolver::secondMomentSolver
    (
        const fvMesh& mesh,
        const abstractProgressVariableThermo& thermo,
        const word& variableName,
        const dictionary& dict,
        multivariateSurfaceInterpolationScheme<scalar>::fieldTable& mvtab
    )
        :
          autoPtr<volScalarField>
        (
          new volScalarField
          (
            IOobject
            (
                variableName & "_MeanSqr",
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
          )
        ),
	mesh_(mesh),
        thermo_(thermo),
        varName_(variableName)
    {
    }

    void secondMomentSolver::registerFields(compressible::LESModel& model)
    {
      //model.registerScalarField((*this)());
    }
}
