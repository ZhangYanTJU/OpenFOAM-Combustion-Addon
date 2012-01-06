/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    presumedPdfOodles

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "abstractProgressVariableThermo.H"
#include "LESModel.H"
#include "IFstream.H"
#include "OFstream.H"
#include "ignition/ignitionSite.H"

//#define PDYN

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMeshNoClear.H"
#   include "readGravitationalAcceleration.H"
#   include "createFields.H"
#   include "readPISOControls.H"
#   include "initContinuityErrs.H"
#   include "readTimeControls.H"
#   include "compressibleCourantNo.H"
#   include "setInitialDeltaT.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readTimeControls.H"
#       include "readPISOControls.H"
#       include "compressibleCourantNo.H"
#       include "setDeltaT.H"

	runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "rhoEqn.H"

        turbulence->correct();

#       include "UEqn.H"
/*
            scalar Ka=0.07788;
            scalar delta_l=0.00033;
            dimensionedScalar D("D", dimless, 2.+(1./3.)*Foam::erf(2.*Ka));
            dimensionedScalar eta("eta", dimLength, 
                (0.345*Foam::pow(Ka, -2)*Foam::exp(-Ka)
                +6.41*Foam::pow(Ka,-0.5)*(1.-Foam::exp(-Ka)))*delta_l);
            dimensionedScalar theta("theta", dimless, 2.5);
            dimensionedScalar f("f", dimless, 0.26);
*/
        // --- PISO loop
        for (int corr=1; corr<=nCorr; corr++)
        {

            bool patched=false;
            forAll(ignitionSites, I)
            {
                //Pout<<"patching ignition site #"<<I<<endl;
                patched = 
                    ignitionSites[I].patchProgressVariable
                    (
                        runTime,
                        thermo.firstMoment
                        (
                            ignitionSites[I].progressVariableID()
                        )
                    ) || patched;
            }

            if (patched) thermo.boundProgressVariables();


            thermo.solveProgressVariables(rho, phi, turbulence());
            thermo.correct();
/*
            const volScalarField& c = thermo->firstMoment(0);
            Sigma=mag(fvc::grad(c))*(exp(-theta*turbulence->delta()/eta)
             +(1.-f*exp(-theta*turbulence->delta()/eta))*pow(turbulence->delta()/eta, D-2.));
*/
            Info 
              <<  "Tmin=" << min(T).value() 
              << " Tmax=" << max(T).value() <<endl;

#           include "pEqn.H"
        }

        runTime.write();

        rho = thermo.rho();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
