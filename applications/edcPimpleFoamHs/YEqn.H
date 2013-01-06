tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);


{
  chemistrySh=dimensionedScalar("chemistrySh", dimEnergy/dimTime/dimVolume, 0.0);

  const volScalarField& Yfu=composition.Y(fuel);
  const volScalarField& Yox=composition.Y(ox);

  tmp<volScalarField> k=turbulence->k();
  volScalarField omega_fu=
      Cmag * rho * (turbulence->epsilon() / (k() + dimensionedScalar("", k().dimensions(), SMALL))) *
      min(Yfu, Yox/s);

  label inertIndex = -1;
  volScalarField Yt("Yt", 0.0*Y[0]);

  for (label i = 0; i < Y.size(); i++)
  {
      scalar f=0.0;

      if (Y[i].name() != inertSpecie)
      {

          if (Y[i].name() == fuel)
            {
              f=-1.0;
            }
          else if (Y[i].name() == ox)
            {
              f=-s;
            }
          else if ( prod.found(Y[i].name()) )
            {
              f=(1.+s)*prod[Y[i].name()];
            }

          forAll(chemistrySh, cellI)
          {
              scalar hi = specieThermo[i].Hc();
              chemistrySh[cellI] += hi*omega_fu[cellI]*f;
          }

          volScalarField& Yi = Y[i];

          fvScalarMatrix YEqn
          (
   	      fvm::ddt(rho, Yi)
            + mvConvection->fvmDiv(phi, Yi)
            - fvm::laplacian(turbulence->muEff(), Yi)
           ==
              omega_fu*f
          );

	  if (oCorr == nOuterCorr-1)
	    {
	      YEqn.relax(mesh.relaxationFactor("Yi"));
	      YEqn.solve(mesh.solver("YiFinal"));
	    }
	  else
	    {
	      YEqn.relax(mesh.relaxationFactor("Yi"));
	      YEqn.solve(mesh.solver("Yi"));
	    }

          Info<<Y[i].name()<<" min/max : "<<min(Yi).value()<<" / "<<max(Yi).value()<<endl;
          Yi.max(0.0);
          Yi.min(1.0);
          Yt += Yi;
      }
      else
      {
        inertIndex = i;
      }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Info<<Y[inertIndex].name()<<" min/max : "<<min(Y[inertIndex]).value()<<" / "<<max(Y[inertIndex]).value()<<endl;
    Y[inertIndex].max(0.0);
    Y[inertIndex].min(1.0);

}
