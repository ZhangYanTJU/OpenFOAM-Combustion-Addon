rho = thermo.rho();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();

#ifdef SIMPLEC
volScalarField p0(p);

volScalarField AU(UEqn().A());
volScalarField AtU(AU - UEqn().H1());
U = UEqn().H()/AU;
#else
volScalarField rAU = 1.0/UEqn().A();
U = rAU*UEqn().H();
#endif

UEqn.clear();

if (transonic)
{
#ifndef SIMPLEC
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)
       *(
            (fvc::interpolate(U) & mesh.Sf())
        )
    );
#endif

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {

#ifdef SIMPLEC
        surfaceScalarField phid
        (
            "phid",
            fvc::interpolate(psi*U) & mesh.Sf()
        );

        surfaceScalarField phic
        (
            "phic",
            fvc::interpolate(rho/AtU - rho/AU)*fvc::snGrad(p)*mesh.magSf()
          + phid*(fvc::interpolate(p) - fvc::interpolate(p, "UD"))
        );
#endif

        fvScalarMatrix pEqn
        (
            fvm::div(phid, p)
#ifdef SIMPLEC
          + fvc::div(phic)
          - fvm::Sp(fvc::div(phid), p)
          + fvc::div(phid)*p
          - fvm::laplacian(rho/AtU, p)
#else
          - fvm::laplacian(rho*rAU, p)
#endif
        );

#ifndef SIMPLEC
        // Relax the pressure equation to ensure diagonal-dominance
        pEqn.relax(mesh.relaxationFactor("pEqn"));
#endif

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
#ifdef SIMPLEC
            phi == phic + pEqn.flux();
#else
            phi == pEqn.flux();
#endif
        }
    }
}
else
{
    phi =
        fvc::interpolate(rho)
       *(
            (fvc::interpolate(U) & mesh.Sf())
        );

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvc::div(phi)
            -fvm::laplacian(rho/AU, p)
        );

        //pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            phi += pEqn.flux();
        }
    }
}

#include "compressibleContinuityErrs.H"

p.relax();

p.max(pMin);

#ifdef SIMPLEC
U -= (fvc::grad(p0)*(1.0/AU - 1.0/AtU) + fvc::grad(p)/AtU);
#else
U -= (1.0/AU)*fvc::grad(p);
#endif
U.correctBoundaryConditions();

// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
//if (closedVolume)
//{
//    p += (initialMass - fvc::domainIntegrate(psi*p))
//        /fvc::domainIntegrate(psi);
//}

rho = thermo.rho();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
if (!transonic) rho.relax();
Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value() << endl;
