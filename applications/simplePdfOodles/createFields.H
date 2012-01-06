    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<abstractProgressVariableThermo> pthermo
    (
        abstractProgressVariableThermo::New(mesh)
    );
    abstractProgressVariableThermo& thermo=pthermo();

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

    volScalarField& p = thermo.p();
    const volScalarField& psi = thermo.psi();
    const volScalarField& T = thermo.T();


    Info<< "\nReading field U\n" << endl;
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

#   include "compressibleCreatePhi.H"


    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::LESModel> turbulence
    (
        compressible::LESModel::New(rho, U, phi, thermo)
    );

    /*
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

    if (composition.contains("ft"))
    {
        fields.add(composition.Y("ft"));
    }

    fields.add(b);
    fields.add(h);
    fields.add(hu);
    */
