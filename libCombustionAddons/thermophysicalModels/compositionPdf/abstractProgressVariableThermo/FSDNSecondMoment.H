#ifndef FSDNSECONDMOMENT_H
#define FSDNSECONDMOMENT_H

#include "secondMoment.H"


namespace Foam
{

class FSDNSecondMoment
: public secondMomentSolver
{

    scalar Ka_;
    scalar theta_, f_;
    scalar D_;
    dimensionedScalar eta_;

    volScalarField Sigma_;

public:
    TypeName("FSD_N");

    FSDNSecondMoment
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
