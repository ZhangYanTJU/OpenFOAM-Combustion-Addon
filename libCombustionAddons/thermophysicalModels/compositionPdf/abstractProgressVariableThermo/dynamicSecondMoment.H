#ifndef DYNAMICSECONDMOMENT_H
#define DYNAMICSECONDMOMENT_H

#include "secondMoment.H"
#include "LESfilter.H"


namespace Foam
{

class dynamicSecondMoment
: public secondMomentSolver
{

    autoPtr<LESfilter> gridFilterPtr_;
    autoPtr<LESfilter> testFilterPtr_;

public:
    TypeName("dynamic");

    dynamicSecondMoment
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
