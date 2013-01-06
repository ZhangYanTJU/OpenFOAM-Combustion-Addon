    tmp<fvVectorMatrix> UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
     ==
        rho*g
    );

if (oCorr == nOuterCorr-1)
{
    UEqn().relax(1);
}
else
{
    UEqn().relax();
}

if (momentumPredictor)
{
    if (oCorr == nOuterCorr-1)
    {
        solve(UEqn() == -fvc::grad(p), mesh.solver("UFinal"));
    }
    else
    {
        solve(UEqn() == -fvc::grad(p));
    }
}
else
{
  volScalarField rAU = 1.0/UEqn().A();
    U = rAU*(UEqn().H() - fvc::grad(p));
    U.correctBoundaryConditions();
}
