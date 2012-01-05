#ifndef FSDSECONDMOMENT_H
#define FSDSECONDMOMENT_H

#include "secondMoment.H"


namespace Foam
{

template<class Thermo>
class FSDSecondMoment
: public secondMomentSolver
{

    //scalar C_;

protected:
    volScalarField Sigma_;

public:
    //TypeName("FSD");

    FSDSecondMoment
    (
        const fvMesh& mesh,
        const abstractProgressVariableThermo& thermo,
        const word& variableName,
        const dictionary& dict,
        multivariateSurfaceInterpolationScheme<scalar>::fieldTable&
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
    ) =0;

    virtual void solve
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

#include "FSDSecondMoment.C"

#endif
