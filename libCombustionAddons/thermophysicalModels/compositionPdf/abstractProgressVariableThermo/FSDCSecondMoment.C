#include "FSDCSecondMoment.H"
#include "addToRunTimeSelectionTable.H"
#include "SIJPDFthermo.H"

namespace Foam
{

defineTypeNameAndDebug(FSDCSecondMoment, 0);
addToRunTimeSelectionTable(secondMomentSolver, FSDCSecondMoment, dictionary);


FSDCSecondMoment::FSDCSecondMoment
(
    const fvMesh& mesh,
    const abstractProgressVariableThermo& thermo,
    const word& variableName,
    const dictionary& dict
)
    : secondMomentSolver(mesh, thermo, variableName, dict),
      ReT_(readScalar(dict.lookup("ReT"))),
      sL_(dict.lookup("sL")),
      deltaL_(dict.lookup("flameThickness")),
      cms_(0.28),
      beta_(1.0),
      Sigma_
    (
        IOobject
        (
            "Sigma",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless/dimLength, 0.0),
        zeroGradientFvPatchScalarField::typeName
    )
{
    if (dict.found("cms")) cms_=readScalar(dict.lookup("cms"));
    if (dict.found("beta")) beta_=readScalar(dict.lookup("beta"));

    Info<<"FSDSecondMoment: using"
        <<" ReT="<<ReT_
        <<" sL="<<sL_.value()
        <<" cms="<<cms_
        <<" beta="<<beta_<<endl;
}


#define BOUND 1e-5


void FSDCSecondMoment::updateSigma
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
    const SIJPDFthermo& pthermo=
        static_cast<const SIJPDFthermo&>(thermo());

    scalar nval=pthermo.constNormValues()[pthermo.variable()[varName()]];

    volScalarField uSl=sqrt(model.k())/sL_;
    scalar alpha=beta_*2.*Foam::log(2.)/(3.*cms_*(Foam::sqrt(ReT_)-1.));
    volScalarField Gamma=0.75*Foam::exp(-1.2/Foam::pow(uSl, 0.3))
        *Foam::pow(model.delta()/deltaL_, 2./3.);

    Sigma_ =
        // == FSDC ==
        mag(fvc::grad(firstMoment/nval))
        *
        (1.+alpha*Gamma*uSl);
}

}
