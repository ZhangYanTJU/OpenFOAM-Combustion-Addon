#ifndef FSDWSECONDMOMENT_H
#define FSDWSECONDMOMENT_H

#include "FSDSecondMoment.H"


namespace Foam
{

template<class Thermo>
class FSDWSecondMoment
: public FSDSecondMoment<Thermo>
{

    dimensionedScalar uPrimeCoef_, XiShapeCoef_, XiCoef_;
    dimensionedScalar sL_;

public:
    TypeName("FSD_W");

    FSDWSecondMoment
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
    );

};

}

#include "FSDWSecondMoment.C"

#endif
