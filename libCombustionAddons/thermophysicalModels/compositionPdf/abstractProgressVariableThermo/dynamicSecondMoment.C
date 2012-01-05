#include "dynamicSecondMoment.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(dynamicSecondMoment, 0);
addToRunTimeSelectionTable(secondMomentSolver, dynamicSecondMoment, dictionary);


dynamicSecondMoment::dynamicSecondMoment
(
    const fvMesh& mesh,
    const abstractProgressVariableThermo& thermo,
    const word& variableName,
    const dictionary& dict,
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable& mvtab
)
    : secondMomentSolver(mesh, thermo, variableName, dict, mvtab),
      gridFilterPtr_(LESfilter::New(mesh, dict, "gridFilter")),
      testFilterPtr_(LESfilter::New(mesh, dict, "testFilter"))
{
}


void dynamicSecondMoment::solve
(
    const volScalarField& rho,
    const surfaceScalarField& phi,
    const volScalarField& firstMoment,
    const volScalarField& meanWdotC,
    fv::convectionScheme<scalar>& mvConvection,
    const compressible::LESModel& model,
    const chemistryTable& preIntegratedTable
)
{
    const LESfilter& gridFilter=gridFilterPtr_();
    const LESfilter& testFilter=testFilterPtr_();
    
    volScalarField pBarSqr=sqr(firstMoment);
    volScalarField pBarBarSqr=sqr(gridFilter(firstMoment));
    volScalarField pBarHatSqr=sqr(testFilter(firstMoment));
    volScalarField pBarBarHatSqr=sqr(testFilter(gridFilter(firstMoment)));
    
    volScalarField k=
        (testFilter(pBarSqr)-pBarHatSqr)
        /
        stabilise
        (
            testFilter(pBarBarSqr)-pBarBarHatSqr, 
            dimensionedScalar("", dimless, SMALL)
        );
    
    (*this)()=(pBarSqr+k*(gridFilter(pBarSqr)-pBarBarSqr));
}

}
