#include "sourceTerm.H"

namespace Foam
{
    
    defineTypeNameAndDebug(sourceTerm, 0);
    defineRunTimeSelectionTable(sourceTerm, dictionary);


    autoPtr<sourceTerm> sourceTerm::New
    (
        //const volScalarField& b
        const hCombustionThermo& thermo
    )
    {
        word sourceTermTypeName;

        // Enclose the creation of the dictionary to ensure it is
        // deleted before the model is created otherwise the dictionary
        // is entered in the database twice
        {
            IOdictionary sourceTermPropertiesDict
                (
                    IOobject
                    (
                        "sourceTermProperties",
                        thermo.T().time().constant(),
                        thermo.T().db(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                );
            
            sourceTermPropertiesDict.lookup("sourceTermModel") 
                >> sourceTermTypeName;
        }

        dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(sourceTermTypeName);
        
        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorIn
                (
                    "sourceTerm::New()"
                )   << "Unknown LESmodel type " << sourceTermTypeName
                    << endl << endl
                    << "Valid sourceTerm types are :" << endl
                    << dictionaryConstructorTablePtr_->toc()
                    << exit(FatalError);
        }

        return autoPtr<sourceTerm>(cstrIter()(thermo));
        
    }



    sourceTerm::sourceTerm
    (
        const word& type,
        //const volScalarField& b
        const hCombustionThermo& thermo
    )
        : IOdictionary
    (
        IOobject
        (
            "sourceTermProperties",
            thermo.T().time().constant(),
            thermo.T().db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
          coeffsDict_(
              subDict(
                  type + "Coeffs"
              )
          ),
          thermo_(thermo),
          omega_
        (
            IOobject
            (
                "omega",
                thermo.T().time().timeName(),
                thermo.T().mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            thermo.T().mesh(),
            dimensionedScalar("", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
        )
    {
    }

    sourceTerm::~sourceTerm()
    {
    }

}
