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

Class
    combustionProgressFunctionObject

Description
    Abstract base-class for Time/database function objects.

SourceFiles
    combustionProgressFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef combustionProgressFunctionObject_H
#define combustionProgressFunctionObject_H

#include "functionObject.H"
#include "fvCFD.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class combustionProgressFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class combustionProgressFunctionObject
: public functionObject
{

    word name_;
    const Time& time_;
    const fvMesh& mesh_;

    const volScalarField& T_;
    dimensionedScalar Tthreshold_;
    dimensionedScalar cubeLength_;

    autoPtr<OFstream> f_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        combustionProgressFunctionObject(const combustionProgressFunctionObject&);

        //- Disallow default bitwise assignment
        void operator=(const combustionProgressFunctionObject&);


public:
    TypeName("combustionProgress");

    // Constructors
    
    //- Construct null
    //combustionProgressFunctionObject();
    
    combustionProgressFunctionObject
    (
        const word& name,
        const Time& t, 
        const dictionary& dict
    );

    // Destructor    
    virtual ~combustionProgressFunctionObject();
    
    
    // Member Functions
    
    //- start is called at the start of the time-loop
    virtual bool start();
    
    //- execute is called at each ++ or += of the time-loop
    virtual bool execute();
    
    //- Read and set the function object if it's data has changed
    virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
