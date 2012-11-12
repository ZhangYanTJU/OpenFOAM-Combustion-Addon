if (composition.contains("ft"))
{
    volScalarField& ft = composition.Y("ft");

    fvScalarMatrix ftEqn
    (
        fvm::ddt(rho, ft)
      + mvConvection->fvmDiv(phi, ft)
      - fvm::laplacian(turbulence->muEff(), ft)
    );

    ftEqn.relax();
    ftEqn.solve();

    ft.max(0.0);
    ft.min(1.0);
}
