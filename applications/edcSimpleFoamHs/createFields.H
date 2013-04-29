Info<< nl << "Reading thermophysicalProperties" << endl;
autoPtr<hsCombustionThermo> thermoPtr
(
	hsCombustionThermo::NewType(mesh, "hsPsiMixtureThermo<reactingMixture<gasThermoPhysics>>")
);
hsCombustionThermo& thermo = thermoPtr();

PtrList<gasThermoPhysics> specieThermo
(
	dynamic_cast<const reactingMixture<gasThermoPhysics>&>(thermo).speciesData()
);

basicMultiComponentMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

word inertSpecie(thermo.lookup("inertSpecie"));

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


volScalarField& p = thermo.p();
const volScalarField& psi = thermo.psi();
volScalarField& hs = thermo.hs();
const volScalarField& T = thermo.T();

#include "compressibleCreatePhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);

dimensionedScalar rhoMax
(
    mesh.solutionDict().subDict("SIMPLE").lookup("rhoMax")
);

dimensionedScalar rhoMin
(
    mesh.solutionDict().subDict("SIMPLE").lookup("rhoMin")
);

dimensionedScalar pMin
(
    mesh.solutionDict().subDict("SIMPLE").lookup("pMin")
);


Info << "Creating turbulence model.\n" << nl;
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

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(hs);

DimensionedField<scalar, volMesh> chemistrySh
(
    IOobject
    (
        "chemistry::Sh",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("chemistrySh", dimEnergy/dimTime/dimVolume, 0.0)
);

scalar beta=0.0, Cmag=1.0, s;
word fuel, ox;
HashTable<scalar, word> prod;
#include "readEdcProperties.H"

