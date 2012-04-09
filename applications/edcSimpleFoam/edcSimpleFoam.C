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

Description
    the transport velocity is the same as that for the Xi equation.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "hCombustionThermo.H"
#include "turbulenceModel.H"
#include "ignition.H"
#include "IFstream.H"
#include "OFstream.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMeshNoClear.H"
#   include "readGravitationalAcceleration.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"
//#   include "readIgnitionProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readSIMPLEControls.H"

        p.storePrevIter();
        rho.storePrevIter();

        turbulence->correct();

        {

#           include "UEqn.H"

            tmp<fv::convectionScheme<scalar> > mvConvection
                (
                    fv::convectionScheme<scalar>::New
                    (
                        mesh,
                        fields,
                        phi,
                        mesh.divScheme("div(phi,b_ft_h)")
                    )
                );

#           include "ftEqn.H"

#           include "bEqn.H"

#           include "hEqn.H"

            thermo.correct();

#           include "pEqn.H"
        }


        rho = thermo.rho();
        rho.relax();
        Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value() << endl;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
