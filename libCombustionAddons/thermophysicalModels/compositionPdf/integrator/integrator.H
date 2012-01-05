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

Class
    integrator

Description

SourceFiles
    integrator.C

\*---------------------------------------------------------------------------*/

#ifndef integrator_H
#define integrator_H

#include "integrableFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class T>
class integrator
{
    
    static integrationParameter X__;
    static integrationParameter constants__;
    static const integrableFunction<T>* F__;

    static int Integrand
    (
        const int *ndim, const double xx[],
        const int *ncomp, double ff[], void *userdata
    );

    static double octaveIntegrand(double);

    static double gsl_integrand(double x, void * params);
public:
  
  /**
  integrate over unit hypercube
  */  
    static T integrate
    (
        const integrableFunction<T>&,
        const integrationParameter&,
        const integrationParameter&
    );

  /**
  integrate over unit hypercube
  */  
    static T integrate
    (
        const sumOfIntegrableFunctions<T>&,
        const integrationParameter&,
        const integrationParameter&
    );

};


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "integrator.C"

//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
