#include "transportSecondMoment.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(transportSecondMoment, 0);
addToRunTimeSelectionTable(secondMomentSolver, transportSecondMoment, dictionary);


transportSecondMoment::transportSecondMoment
(
    const fvMesh& mesh,
    const abstractProgressVariableThermo& thermo,
    const word& variableName,
    const dictionary& dict,
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable& mvtab
)
    : secondMomentSolver(mesh, thermo, variableName, dict, mvtab),
      C_D(readScalar(dict.lookup("C_D"))),
      Sc_t(readScalar(dict.lookup("Sc_t"))) /*,
SGS_dissipation
(
 IOobject
(
variableName & "_SGSDISS",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
dimensionedScalar("", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
),
resolved_dissipation
(
 IOobject
(
variableName & "_RESDISS",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
dimensionedScalar("", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
)*/
{
  mvtab.add( (*this)() );
}


void transportSecondMoment::solve
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
/*
resolved_dissipation=-2.0*model.alpha()*magSqr(fvc::grad(firstMoment));
SGS_dissipation=-2.0*C_D*(model.alphaEff()/ (sqr(model.delta())*Sc_t) )
                *( (*this)() - sqr(firstMoment) );
*/
    Foam::solve
        (
            fvm::ddt(rho, (*this)() )
            + mvConvection.fvmDiv(phi, (*this)() )
            - 2.0*rho*meanWdotC
            ==
            //model.divRhoSgsFlux( (*this)() )
	    fvm::laplacian(model.alphaEff(), (*this)(), 
			   "laplacian(Deff,"+(*this)().name()+")")
            -2.0*
            (
                model.alpha()*magSqr(fvc::grad(firstMoment))
                +
                C_D*(model.muSgs()/ (sqr(model.delta())*Sc_t) )
                *( (*this)() - sqr(firstMoment) )
            )
        );
}

}
