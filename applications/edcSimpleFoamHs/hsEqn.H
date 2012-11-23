{
    fvScalarMatrix hsEqn
    (
        mvConvection->fvmDiv(phi, hs)
      - fvm::laplacian(turbulence->alphaEff(), hs)
//      - fvm::laplacian(turbulence->muEff(), hs)  // unit lewis no.
     ==
//        DpDt
       chemistrySh
    );

    hsEqn.relax();
    hsEqn.solve();

    if (runTime.outputTime()) hs.write();

    thermo.correct();

    Info<< "T gas min/max   = " << min(T).value() << ", "
        << max(T).value() << endl;
}
