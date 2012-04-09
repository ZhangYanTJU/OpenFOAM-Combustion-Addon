volScalarField Ymin = 0.0*b;

if (composition.contains("ft"))
{
  volScalarField& ft = composition.Y("ft");

  volScalarField Yfu = b*ft + 
    (1.0 - b)*max(ft - (1.0 - ft)/stoicRatio.value(), 0.0);

  volScalarField Yox = 1.0 - ft - (ft - Yfu)*stoicRatio.value();

  volScalarField Ypr = 1.0 - Yfu - Yox;

  Ymin=Yfu;
  Ymin=min(Ymin, Yox/stoicRatio.value());
  //Ymin=min(Ymin, beta*Ypr/(1.+stoicRatio.value()));
}

omega=Cmag * rho * (turbulence->epsilon()/turbulence->k()) * Ymin;

fvScalarMatrix bEqn
(
 fvm::ddt(rho, b)
 + mvConvection->fvmDiv(phi, b)
 ==
 fvm::laplacian(turbulence->muEff(), b)
 - omega
);

bEqn.relax();
bEqn.solve();
Info<< "min(b) = " << min(b).value() << endl;

b.min(1.0);
b.max(0.0);
