Info<< nl << "Reading thermophysicalProperties" << endl;
autoPtr<hCombustionThermo> thermoPtr
(
    hCombustionThermo::New(mesh)
);
hCombustionThermo& thermo = thermoPtr();

dimensionedScalar stoicRatio(thermo.lookup("stoichiometricAirFuelMassRatio"));

basicMultiComponentMixture& composition = thermo.composition();
volScalarField& p = thermo.p();
const volScalarField& psi = thermo.psi();
volScalarField& h = thermo.h();
const volScalarField& T = thermo.T();

// reference to regress variable
volScalarField& b = composition.Y("b");

volScalarField omega
(
    IOobject
    (
        "omega_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("", dimensionSet(1, -3, -1, 0, 0), 0.0)
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


#include "compressibleCreatePhi.H"

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
 compressible::turbulenceModel::New
 (
  rho,
  U,
  phi,
  thermo
 )
);


Info<< "Creating field DpDt\n" << endl;
volScalarField DpDt = 
    fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

if (composition.contains("ft"))
  {
    fields.add(composition.Y("ft"));
  }


fields.add(b);
fields.add(h);

scalar beta, Cmag;
{
  IOdictionary edcProperties
    (
     IOobject
     (
      "edcProperties",
      runTime.constant(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
     )
    );

  beta=readScalar(edcProperties.lookup("beta"));
  Cmag=readScalar(edcProperties.lookup("Cmag"));
}
