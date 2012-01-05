#include "FSDWSecondMoment.H"
#include "addToRunTimeSelectionTable.H"
#include "SIJPDFthermo.H"

namespace Foam
{

    //defineTypeNameAndDebug(FSDWSecondMoment, 0);
    //addToRunTimeSelectionTable(secondMomentSolver, FSDWSecondMoment, dictionary);

template<class Thermo>
FSDWSecondMoment<Thermo>::FSDWSecondMoment
(
    const fvMesh& mesh,
    const abstractProgressVariableThermo& thermo,
    const word& variableName,
    const dictionary& dict,
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable& mvtab
)
    : FSDSecondMoment<Thermo>(mesh, thermo, variableName, dict, mvtab),
      uPrimeCoef_(dict.lookup("uPrimeCoef")),
      XiShapeCoef_(dict.lookup("XiShapeCoef")), 
      XiCoef_(dict.lookup("XiCoef")),
      sL_(dict.lookup("sL"))
{
}


#define BOUND 1e-5


template<class Thermo>
void FSDWSecondMoment<Thermo>::updateSigma
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
    const Thermo& pthermo=
        static_cast<const Thermo&>(this->thermo());

    label vI=pthermo.variable()[this->varName()];
    //scalar nval=pthermo.constNormValues()[pthermo.variable()[this->varName()]];
 
    volScalarField up = uPrimeCoef_*sqrt(max
     (
      (2.0/3.0)*model.k(), 
      dimensionedScalar("", dimVelocity*dimVelocity, SMALL))
     );
    volScalarField epsilon = pow(uPrimeCoef_, 3)*model.epsilon();
    volScalarField tauetarad=model.mu()/max(rho*epsilon, dimensionedScalar("", rho.dimensions()*epsilon.dimensions(), SMALL));
    volScalarField tauEta = sqrt(max
     (
      tauetarad, 
      dimensionedScalar("", tauetarad.dimensions(), SMALL)
     )
     );
    volScalarField Reta = up/
    (
        sqrt(epsilon*tauEta) + dimensionedScalar("1e-8", up.dimensions(), 1e-8)
    );

    tmp<volScalarField> nM1=pthermo.normalizedFirstMoment(vI);
    this->Sigma_ = mag(fvc::grad(nM1))
            *
            (
             scalar(1) 
              +
             ( scalar(1) + (2.*XiShapeCoef_)*(nM1()) - 0.5 )
            *XiCoef_*sqrt(up/sL_)*Reta);

 
}

}
