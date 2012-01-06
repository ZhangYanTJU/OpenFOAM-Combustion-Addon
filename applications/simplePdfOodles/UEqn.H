    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoBeff(U)
     ==
        rho*g
    );

    UEqn.relax();

    if (momentumPredictor)
    {
        solve(UEqn == -fvc::grad(p));
    }
