#include "colinWrinkling.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(colinWrinkling, 0);
addToRunTimeSelectionTable(efficiencyFunction, colinWrinkling, dictionary);


colinWrinkling::colinWrinkling
(
    const volVectorField& U,
    const scalar& TF
)
    : efficiencyFunction(U, TF),
      uPrime_
    (
        IOobject
        (
            "uPrime",
            U_.time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("", dimLength/dimTime, 0.0)
    ),
      alpha_(readScalar(lookup("alpha"))),
      deltal0_(lookup("deltal0")),
      sl0_(lookup("sl0"))
{
}

colinWrinkling::~colinWrinkling()
{
}

void colinWrinkling::correct()
{
    uPrime_ = 2.0 * mag(pow(delta(), 3) * fvc::laplacian(fvc::curl(U_)));

    volScalarField RL = 2.0*delta() / deltal0_;
    volScalarField RV = uPrime_ / sl0_;

    efficiency_ = 
        (1.0 + alpha_ * 0.75*exp(-1.2/pow(RV, 0.3))*pow(RL, 2./3.) * RV)
        /
        (1.0 + alpha_ * 0.75*exp(-1.2/pow(RV, 0.3))*pow(RL/TF_, 2./3.) * RV);
}

}
