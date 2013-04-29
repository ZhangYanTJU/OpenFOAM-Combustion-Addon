rho = thermo.rho();
//rho = max(rho, rhoMin);
//rho = min(rho, rhoMax);
//rho.relax();

volScalarField rAU = 1.0/UEqn().A();
U = rAU*UEqn().H();

if (nCorr <= 1)
{
    UEqn.clear();
}

scalar tmppResidual=0;

if (transonic)
{
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)
       *(
            (fvc::interpolate(U) & mesh.Sf())
	    + fvc::ddtPhiCorr(rAU, rho, U, phi)
        )
    );

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {

        fvScalarMatrix pEqn
        (
	    fvm::ddt(psi, p)
          + fvm::div(phid, p)
          - fvm::laplacian(rho*rAU, p)
        );

	pEqn.relax(mesh.relaxationFactor("pEqn"));
        //pEqn.setReference(pRefCell, pRefValue);

        if
        (
            oCorr == nOuterCorr-1
            && corr == nCorr-1
            && nonOrth == nNonOrthCorr
        )
        {
	  tmppResidual=pEqn.solve(mesh.solver("pFinal")).initialResidual();
        }
        else
        {
	  tmppResidual=pEqn.solve().initialResidual();
        }


        if (nonOrth == nNonOrthCorr)
        {
            phi == pEqn.flux();
	}
        if (nonOrth == 0)
        {
	    pResidual=tmppResidual;
        }
    }
}
else
{
    phi =
        fvc::interpolate(rho)
       *(
            (fvc::interpolate(U) & mesh.Sf())
        );

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
             fvm::ddt(psi, p)
          +  fvc::div(phi)
          -  fvm::laplacian(rho*rAU, p)
        );

        pEqn.setReference(pRefCell, pRefValue);

        if
        (
            oCorr == nOuterCorr-1
         && corr == nCorr-1
         && nonOrth == nNonOrthCorr
        )
        {
	  tmppResidual=pEqn.solve(mesh.solver("pFinal")).initialResidual();
        }
        else
        {
	  tmppResidual=pEqn.solve().initialResidual();
        }

        if (nonOrth == nNonOrthCorr)
        {
            phi += pEqn.flux();
	}
        if (nonOrth == 0)
        {
	    pResidual=tmppResidual;
        }
    }
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

//if (oCorr != nOuterCorr-1)
{
    // Explicitly relax pressure for momentum corrector
    p.relax();

    rho = thermo.rho();
    rho = max(rho, rhoMin);
    rho = min(rho, rhoMax);
    rho.relax();
    Info<< "rho max/min : " << max(rho).value()
        << " " << min(rho).value() << endl;
}

U -= rAU*fvc::grad(p);
U.correctBoundaryConditions();

DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

bound(p, pMin);

Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value() << endl;
