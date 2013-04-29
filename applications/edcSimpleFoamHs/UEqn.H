    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
     ==
        rho*g
    );

    UEqn().relax();

    solve(UEqn() == -fvc::grad(p));
