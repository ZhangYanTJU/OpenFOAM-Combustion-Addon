{
    fvScalarMatrix hsEqn
    (
        mvConvection->fvmDiv(phi, hs)
      - fvm::Sp(fvc::div(phi), hs)
      - fvm::laplacian(turbulence->alphaEff(), hs)
//      - fvm::laplacian(turbulence->muEff(), hs)  // unit lewis no.
      - chemistrySh
     ==
//        DpDt
/*
        fvc::div(phi/fvc::interpolate(rho), rho/psi, "div(U,p)")
      - (rho/psi)*fvc::div(phi/fvc::interpolate(rho))
*/
      - fvc::div(phi, 0.5*magSqr(U), "div(phi,K)")
    );

    hsEqn.relax();
    hsEqn.solve();

//    if (runTime.outputTime()) hs.write();

    thermo.correct();

    Info<< "T gas min/max   = " << min(T).value() << ", "
        << max(T).value() << endl;
}
