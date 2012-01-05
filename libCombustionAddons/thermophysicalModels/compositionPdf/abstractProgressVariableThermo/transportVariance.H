#ifndef TRANSPORTVARIANCE_H
#define TRANSPORTVARIANCE_H

#include "secondMoment.H"


namespace Foam
{

class transportVariance
: public secondMomentSolver
{

    scalar C_D, Sc_t;
    volScalarField variance_;

public:
    TypeName("transportVariance");

    transportVariance
    (
        const fvMesh& mesh,
        const abstractProgressVariableThermo& thermo,
        const word& variableName,
        const dictionary& dict,
        multivariateSurfaceInterpolationScheme<scalar>::fieldTable&
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
    );

    void registerFields(compressible::LESModel& model);
};

}

#endif
