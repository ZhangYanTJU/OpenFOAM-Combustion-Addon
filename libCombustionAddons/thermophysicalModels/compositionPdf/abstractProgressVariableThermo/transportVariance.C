#include "transportVariance.H"
#include "SIJPDFthermoStateFinder.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(transportVariance, 0);
addToRunTimeSelectionTable(secondMomentSolver, transportVariance, dictionary);


transportVariance::transportVariance
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
  mvtab.add( variance_ );
}

void  transportVariance::registerFields(compressible::LESModel& model)
{
  //model.registerScalarField(variance_);
}

void transportVariance::solve
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
    const volScalarField& meanWdot=
      (static_cast<const SIJPDFthermo<SIJPDFthermoIndexDriver>&>(thermo())).meanWdot(varName());
    Foam::solve
        (
            fvm::ddt(rho, variance_ )
            + mvConvection.fvmDiv(phi, variance_ )
            - 2.0*rho*(meanWdotC - meanWdot*firstMoment)
            ==
            //model.divRhoSgsFlux( variance_ )
	    fvm::laplacian(model.alphaEff(), 
			   variance_, 
			   "laplacian(Deff,"+variance_.name()+")")
            -fvm::Sp(2.0*C_D*(model.muSgs()/ (sqr(model.delta())*Sc_t) ), variance_)
	    //#warning  sign of model.sgsFlux is wrong! (inverted here)
            //+2.0*(rho*model.sgsFlux(firstMoment)()&fvc::grad(firstMoment)) 
#warning CHECK SIGN!
            -2.0*(model.alphaEff()*fvc::grad(firstMoment)&fvc::grad(firstMoment)) 
        );

     variance_=min(variance_, firstMoment*(1.0-firstMoment));
     variance_.max(0.0);

     (*this)()=variance_+sqr(firstMoment);
}

}
