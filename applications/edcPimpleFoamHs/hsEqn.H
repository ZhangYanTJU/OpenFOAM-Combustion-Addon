{
    fvScalarMatrix hsEqn
    (
        fvm::ddt(rho, hs)
      + mvConvection->fvmDiv(phi, hs)
      - fvm::laplacian(turbulence->alphaEff(), hs)
//      - fvm::laplacian(turbulence->muEff(), hs)  // unit lewis no.
      - chemistrySh
     ==
        DpDt
    );

    if (oCorr == nOuterCorr-1)
    {
        hsEqn.relax();
        hsEqn.solve(mesh.solver("hsFinal"));
    }
    else
    {
        hsEqn.relax();
        hsEqn.solve();
    }


//    if (runTime.outputTime()) hs.write();

    thermo.correct();

    Info<< "T gas min/max   = " << min(T).value() << ", "
        << max(T).value() << endl;
}
