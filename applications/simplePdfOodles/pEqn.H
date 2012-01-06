rho = thermo.rho();

volScalarField rUA = 1.0/UEqn.A();
U = rUA*UEqn.H();

if (transonic)
{
    surfaceScalarField phid =
    (
        (fvc::interpolate(rho*U) & mesh.Sf())
      + fvc::ddtPhiCorr(rUA, rho, U, phi)
    )/fvc::interpolate(p);

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvm::div(phid, p, "div(phid,p)")
          - fvm::laplacian(rho*rUA, p)
        );

        pEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            phi == pEqn.flux();
        }
    }
}
else
{
    phi = (fvc::interpolate(rho*U) & mesh.Sf())
        + fvc::ddtPhiCorr(rUA, rho, U, phi);

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvc::div(phi)
          - fvm::laplacian(rho*rUA, p)
        );

        pEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            phi += pEqn.flux();
        }
    }
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

p.relax();

U -= rUA*fvc::grad(p);
U.correctBoundaryConditions();
