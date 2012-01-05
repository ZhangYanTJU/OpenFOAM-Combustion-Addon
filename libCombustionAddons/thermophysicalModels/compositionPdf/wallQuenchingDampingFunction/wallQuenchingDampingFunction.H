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
    wallQuenchingDampingFunction

Description

Author

SourceFiles
    wallQuenchingDampingFunction.C

\*---------------------------------------------------------------------------*/

#ifndef wallQuenchingDampingFunction_H
#define wallQuenchingDampingFunction_H

#include "fvCFD.H"
#include "cellDistFuncs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class wallQuenchingDampingFunction Declaration
\*---------------------------------------------------------------------------*/

class wallQuenchingDampingFunction
 : public cellDistFuncs
{


    // Private Member Functions
    wallQuenchingDampingFunction(const wallQuenchingDampingFunction&);
    void operator=(const wallQuenchingDampingFunction&);

protected:
    //- Number of unset cells and faces.
    label nUnset_;
    volScalarField y_;

    void correctWallDist(const wordList& patches);

public:
    TypeName("wallQuenchingDampingFunction");


    inline const volScalarField& wallY() const { return y_; }
    // Constructors

    //- Construct from Istream
    wallQuenchingDampingFunction(const fvMesh& mesh, const dictionary&);


   // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        wallQuenchingDampingFunction,
        dictionary,
        (
            const fvMesh& mesh,
            const dictionary& d
        ),
        (mesh, d)
    );

    // Selectors

    //- Select null constructed
    static autoPtr<wallQuenchingDampingFunction> New
    (
            const fvMesh& mesh,
            const dictionary& d,
            word keyword="type"
    );

    virtual ~wallQuenchingDampingFunction();

    virtual tmp<volScalarField> operator()(const volScalarField& fm) const =0;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
