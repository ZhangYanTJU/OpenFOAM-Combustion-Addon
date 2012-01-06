    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoBeff(U)
     ==
        rho*g
    );

    if (momentumPredictor)
    {
#ifdef PDYN
        solve(UEqn == -fvc::grad(pdyn));
#else
        solve(UEqn == -fvc::grad(p));
#endif
    }
