{
    fvScalarMatrix hEqn
    (
        fvm::ddt(rho, h)
      + mvConvection->fvmDiv(phi, h)
      - fvm::laplacian(turbulence->alphaEff(), h)
      //- fvm::laplacian(turbulence->alpha(), h)
/*     ==
        DpDt
*/
    );

    hEqn.relax();
    hEqn.solve();

    thermo.correct();
}
