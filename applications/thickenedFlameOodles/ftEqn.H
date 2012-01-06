if (composition.contains("ft"))
{
    volScalarField& ft = composition.Y("ft");

    solve
    (
        fvm::ddt(rho, ft)
      + mvConvection->fvmDiv(phi, ft)
      - fvm::laplacian(turbulence->muEff(), ft)
    );
}
