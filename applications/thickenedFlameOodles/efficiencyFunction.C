#include "efficiencyFunction.H"

namespace Foam
{

    defineTypeNameAndDebug(efficiencyFunction, 0);
    defineRunTimeSelectionTable(efficiencyFunction, dictionary);

    autoPtr<efficiencyFunction> efficiencyFunction::New
    (
        const volVectorField& U,
        const scalar& TF
    )
    {
        word typeName;

        // Enclose the creation of the dictionary to ensure it is
        // deleted before the model is created otherwise the dictionary
        // is entered in the database twice
        {
            IOdictionary propertiesDict
                (
                    IOobject
                    (
                        "efficiencyFunctionProperties",
                        U.time().constant(),
                        U.db(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                );
            
            propertiesDict.lookup("efficiencyFunction") 
                >> typeName;
        }

        dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(typeName);
        
        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorIn
                (
                    "efficiencyFunction::New()"
                )   << "Unknown efficiency function " << typeName
                    << endl << endl
                    << "Valid efficiency functions are :" << endl
                    << dictionaryConstructorTablePtr_->toc()
                    << exit(FatalError);
        }

        return autoPtr<efficiencyFunction>(cstrIter()(U, TF));
        
    }



    efficiencyFunction::efficiencyFunction
    (
        const volVectorField& U,
        const scalar& TF
    )
        : IOdictionary
    (
        IOobject
        (
            "efficiencyFunctionProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
          U_(U),
          TF_(TF),
          efficiency_
        (
            IOobject
            (
                "efficiency",
                U.time().timeName(),
                U.mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            U.mesh(),
            dimensionedScalar("", dimless, 1.0)
        ),
          delta_(LESdelta::New("delta", U.mesh(), *this))
    {
    }

    efficiencyFunction::~efficiencyFunction()
    {
    }


}
