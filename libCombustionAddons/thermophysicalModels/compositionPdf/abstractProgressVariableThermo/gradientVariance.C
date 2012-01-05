#include "gradientVariance.H"
#include "SIJPDFthermoStateFinder.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(gradientVariance, 0);
addToRunTimeSelectionTable(secondMomentSolver, gradientVariance, dictionary);


gradientVariance::gradientVariance
(
    const fvMesh& mesh,
    const abstractProgressVariableThermo& thermo,
    const word& variableName,
    const dictionary& dict,
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable& mvtab
)
    : secondMomentSolver(mesh, thermo, variableName, dict, mvtab),
      C_D(readScalar(dict.lookup("C_D"))),
      Sc_t(readScalar(dict.lookup("Sc_t"))),
      variance_
      (
       IOobject
       (
        variableName&"_Var",
                        mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
      )
{
}

void  gradientVariance::registerFields(compressible::LESModel& model)
{
 //model.registerScalarField(variance_);
}

void gradientVariance::solve
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
     variance_=C_D*sqr(model.delta())*magSqr(fvc::grad(firstMoment));
    
     variance_=min(variance_, firstMoment*(1.0-firstMoment));
     variance_.max(0.0);

     (*this)()=variance_+sqr(firstMoment);
}

}
