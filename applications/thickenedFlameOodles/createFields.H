Info<< nl << "Reading thermophysicalProperties" << endl;
autoPtr<hCombustionThermo> pthermo
(
    hCombustionThermo::New(mesh)
);
hCombustionThermo& thermo=pthermo();

basicMultiComponentMixture& composition = thermo.composition();
volScalarField& p = thermo.p();
const volScalarField& psi = thermo.psi();
volScalarField& h = thermo.h();
const volScalarField& T = thermo.T();

// reference to regress variable
volScalarField& b = composition.Y("b");

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

autoPtr<sourceTerm> chemistry
(
    sourceTerm::New(thermo)
);

autoPtr<efficiencyFunction> efficiency
(
    efficiencyFunction::New(U, TF)
);


#include "compressibleCreatePhi.H"


Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::LESModel> turbulence
(
    compressible::LESModel::New
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

fields.add(b);
fields.add(h);

