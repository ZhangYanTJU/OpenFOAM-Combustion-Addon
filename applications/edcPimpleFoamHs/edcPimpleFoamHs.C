/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    reactingFoam

Description
    Solver for combustion with chemical reactions.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "hCombustionThermo.H"
#include "turbulenceModel.H"
#include "psiChemistryModel.H"
#include "chemistrySolver.H"
#include "multivariateScheme.H"
#include "thermoPhysicsTypes.H"
#include "bound.H"

#define SIMPLEC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "readGravitationalAcceleration.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readEdcProperties.H"
        #include "readTimeControls.H"
        #include "readPIMPLEControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (nOuterCorr != 1)
        {
            p.storePrevIter();
            rho.storePrevIter();
        }

        #include "rhoEqn.H"

	scalar pResidual=GREAT;

        // --- Pressure-velocity PIMPLE corrector loop
        int oCorr=0;
	for (; oCorr<nOuterCorr; oCorr++)
        {
#           include "UEqn.H"
#           include "YEqn.H"
#           include "hsEqn.H"

            // --- PISO loop
            for (int corr=0; corr<nCorr; corr++)
            {
#              include "pEqn.H"
	    }

	    if (pResidual<minpResidual) // Pressure converged
	      {
		if (adjustDeltaTByConvergence && 
		    (oCorr < 4*nOuterCorr/5))
		  {
		    runTime.setDeltaT(runTime.deltaT().value()/0.8);
		  }
		Info<<"SIMPLE loop converged after "<<oCorr<<" iterations"<<endl;
		break;
	      }

	    turbulence->correct();
        }

	if (pResidual>=minpResidual) // Pressure not converged
	  {
	    if (adjustDeltaTByConvergence && 
		(oCorr >= nOuterCorr))
	      {
		runTime.setDeltaT(runTime.deltaT().value()*0.8);
	      }
	  }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
