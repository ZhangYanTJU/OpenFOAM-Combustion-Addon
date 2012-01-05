#include "FSDNSecondMoment.H"
#include "addToRunTimeSelectionTable.H"
#include "SIJPDFthermo.H"

namespace Foam
{

defineTypeNameAndDebug(FSDNSecondMoment, 0);
addToRunTimeSelectionTable(secondMomentSolver, FSDNSecondMoment, dictionary);


FSDNSecondMoment::FSDNSecondMoment
(
    const fvMesh& mesh,
    const abstractProgressVariableThermo& thermo,
    const word& variableName,
    const dictionary& dict
)
    : secondMomentSolver(mesh, thermo, variableName, dict),
      Ka_(readScalar(dict.lookup("KarlovitzNumber"))),
      theta_(2.5),
      f_(0.26),
      D_(2.+(1./3.)*Foam::erf(2.*Ka_)),
      eta_("eta", dimLength, 0.0),
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
    if (dict.found("theta")) theta_=readScalar(dict.lookup("theta"));
    if (dict.found("f")) f_=readScalar(dict.lookup("f"));
    if (dict.found("eta"))
        eta_=dimensionedScalar(dict.lookup("eta"));
    else
    {
        dimensionedScalar deltaL(dict.lookup("flameThickness"));
        eta_=
        ( 
            0.345*Foam::pow(Ka_, -2)*Foam::exp(-Ka_)
            + 
            6.41*Foam::pow(Ka_, -0.5)*(1.-Foam::exp(-Ka_))
        ) * deltaL;
    }

    Info<<"FSDNSecondMoment: using"
        <<" D="<<D_
        <<" eta_i="<<eta_.value()
        <<" theta="<<theta_
        <<" f="<<f_<<endl;
}


#define BOUND 1e-5


void FSDNSecondMoment::updateSigma
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

    Sigma_ =
        // == FSDN ==
        mag(fvc::grad(firstMoment/nval))
        *
        (
            exp(-theta_*model.delta()/eta_)
            +
            (1.-f_*exp(-theta_*model.delta()/eta_)) * pow(model.delta()/eta_, D_ - 2.)
        );
}

}
