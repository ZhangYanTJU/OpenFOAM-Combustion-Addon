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

/*
    volScalarField Sigma
    (
        IOobject
        (
            "Sigma",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh, 
	dimensionedScalar("", dimless/dimLength, 0.0)
    );
*/

    volScalarField& p = thermo.p();

#ifdef PDYN
    dimensionedScalar pRef(thermo.lookup("referencePressure"));
    volScalarField pdyn
    (
        IOobject
        (
            "pdyn",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    p = pdyn + pRef;
#endif

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

    thermo.registerScalarFields(turbulence());

    Info<< "Reading combustion properties\n" << endl;

    IOdictionary combustionProperties
    (
        IOobject
        (
            "combustionProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
 
    PtrList<ignitionSite> ignitionSites
	(
	 combustionProperties.lookup("ignitionSites"),
         ignitionSite::iNew(mesh)
        );

    forAll(ignitionSites, I)
     Info<<"Ignition #"<<I<<" at "<<ignitionSites[I].location()<<endl;

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
