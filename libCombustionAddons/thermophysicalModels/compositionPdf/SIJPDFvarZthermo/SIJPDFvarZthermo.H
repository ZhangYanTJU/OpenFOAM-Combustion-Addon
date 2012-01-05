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
    SIJPDFvarZthermo

Description
    Thermophysical properties integrated over statistical independent 
    joint PDFs. The value range of the properties varies with mixture
    fraction.

SourceFiles
    SIJPDFvarZthermo.C

\*---------------------------------------------------------------------------*/

#ifndef SIJPDFvarZthermo_H
#define SIJPDFvarZthermo_H

#include "SIJPDFthermo.H"
//#include "SIJPDFvarZthermoStateFinder.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    class SIJPDFvarZthermoIndexDriver;

// Class forward declarations

/*---------------------------------------------------------------------------*\
                           Class SIJPDFvarZthermo Declaration
\*---------------------------------------------------------------------------*/

class SIJPDFvarZthermo
:
    public SIJPDFthermo<SIJPDFvarZthermoIndexDriver>
{

    friend class SIJPDFvarZthermoIndexDriver;
    // Private data

    // Private Member Functions

    //- Disallow default bitwise copy construct
    SIJPDFvarZthermo(const SIJPDFvarZthermo&);
    
    //- Disallow default bitwise assignment
    void operator=(const SIJPDFvarZthermo&);
    
protected:
    word mixFracName_;
    label mixFracI_;
    HashPtrTable<integrableGraph> normValues_;
    
    HashPtrTable<MultidimensionalLookupTable<scalarField> > derivedNormFactors_;
    
    virtual void performIntegration
    (
        //int[],
        List<label>&,
        chemicalSystemState&,
        integrationParameter&,
        integrationParameter&
    );

    virtual void initialize();

    bool readNormFactorsFromCache();

public:
    //- Runtime type information
    TypeName("SIJPDFvarZthermo");


    // Constructors

        //- Construct from components
        SIJPDFvarZthermo(const fvMesh&);



    // Destructor

        ~SIJPDFvarZthermo();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

    // Friend Functions

    // Friend Operators

    // IOstream Operators
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
