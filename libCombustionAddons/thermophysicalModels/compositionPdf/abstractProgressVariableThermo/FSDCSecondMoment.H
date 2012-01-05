#ifndef FSDCSECONDMOMENT_H
#define FSDCSECONDMOMENT_H

#include "secondMoment.H"


namespace Foam
{

class FSDCSecondMoment
: public secondMomentSolver
{

    scalar ReT_;
    dimensionedScalar sL_;
    dimensionedScalar deltaL_;
    scalar uSl_, cms_, beta_;

    volScalarField Sigma_;

public:
    TypeName("FSD_C");

    FSDCSecondMoment
    (
        const fvMesh& mesh,
        const abstractProgressVariableThermo& thermo,
        const word& variableName,
        const dictionary& dict
    );

    virtual void updateSigma
    (
        const volScalarField& rho,
        const surfaceScalarField& phi,
        const volScalarField& firstMoment,
        const volScalarField& meanWdotC,
        fv::convectionScheme<scalar>& mvConvection,
        const compressible::LESModel& model,
        const chemistryTable& preIntegratedTable
    );

};

}

#endif
