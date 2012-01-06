rho = thermo.rho();

volScalarField rUA = 1.0/UEqn.A();
U = rUA*UEqn.H();

if (transonic)
{
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(thermo.psi())
       *(
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rUA, rho, U, phi)
        )
    );

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
#ifdef PDYN
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, pdyn)
          + fvc::ddt(psi)*pRef
          + fvm::div(phid, pdyn)
          - fvm::laplacian(rho*rUA, pdyn)
        );
#else
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvm::div(phid, p)
          - fvm::laplacian(rho*rUA, p)
        );
#endif

        pEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            phi == pEqn.flux();
        }
    }
}
else
{
    phi = 
        fvc::interpolate(rho)*
        (
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rUA, rho, U, phi)
        );

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
#ifdef PDYN
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, pdyn)
          + fvc::ddt(psi)*pRef
          + fvc::div(phi)
          - fvm::laplacian(rho*rUA, pdyn)
        );
#else
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvc::div(phi)
          - fvm::laplacian(rho*rUA, p)
        );
#endif
        pEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            phi += pEqn.flux();
        }
    }
}

#ifdef PDYN
p == pdyn + pRef;
#endif

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

#ifdef PDYN
U -= rUA*fvc::grad(pdyn);
#else
U -= rUA*fvc::grad(p);
#endif

U.correctBoundaryConditions();

