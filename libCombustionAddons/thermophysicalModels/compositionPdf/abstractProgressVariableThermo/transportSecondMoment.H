#ifndef TRANSPORTSECONDMOMENT_H
#define TRANSPORTSECONDMOMENT_H

#include "secondMoment.H"


namespace Foam
{

class transportSecondMoment
: public secondMomentSolver
{

    scalar C_D, Sc_t;
//volScalarField SGS_dissipation, resolved_dissipation;

public:
    TypeName("transport");

    transportSecondMoment
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

};

}

#endif
