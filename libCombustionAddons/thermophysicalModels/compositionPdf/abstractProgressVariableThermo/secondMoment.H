#ifndef SECONDMOMENT_H
#define SECONDMOMENT_H

#include "fvCFD.H"
#include "chemistryTable.H"
#include "compressible/LES/LESModel/LESModel.H"
#include "incompressible/LES/LESModel/LESModel.H"


namespace Foam
{

class abstractProgressVariableThermo;

class secondMomentSolver
: public autoPtr<volScalarField>

{

    const fvMesh& mesh_;
    const abstractProgressVariableThermo& thermo_;
    word varName_;

public:
    TypeName("secondMomentSolver");

    inline const abstractProgressVariableThermo& thermo() const
        { return thermo_; }

    inline const word& varName() const
        { return varName_; }

    // Declare run-time constructor selection table
    
    declareRunTimeSelectionTable
    (
        autoPtr,
        secondMomentSolver,
        dictionary,
        (
            const fvMesh& mesh,
            const abstractProgressVariableThermo& thermo,
            const word& variableName,
            const dictionary& dict,
            multivariateSurfaceInterpolationScheme<scalar>::fieldTable& mvtab
        ),
        (mesh, thermo, variableName, dict, mvtab)
    );

    static autoPtr<secondMomentSolver> New
    (
        const fvMesh& mesh,
        const abstractProgressVariableThermo& thermo,
        const word& variableName,
        const dictionary& dict,
        multivariateSurfaceInterpolationScheme<scalar>::fieldTable& mvtab
    );

    secondMomentSolver
    (
        const fvMesh& mesh,
        const abstractProgressVariableThermo& thermo,
        const word& variableName,
        const dictionary& dict,
        multivariateSurfaceInterpolationScheme<scalar>::fieldTable& mvtab
    );

    virtual void solve
    (
        const volScalarField& rho,
        const surfaceScalarField& phi,
        const volScalarField& firstMoment,
        const volScalarField& meanWdotC,
        fv::convectionScheme<scalar>& mvConvection,
        const compressible::LESModel& model,
        const chemistryTable& preIntegratedTable
    ) =0;

    virtual void registerFields(compressible::LESModel& model);
};

}

#endif
