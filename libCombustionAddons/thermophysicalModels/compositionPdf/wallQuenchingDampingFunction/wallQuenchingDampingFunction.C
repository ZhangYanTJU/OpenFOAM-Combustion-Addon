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

\*----------------------------------------------------------------------------*/

#include "wallQuenchingDampingFunction.H"
#include "wallDist.H"
#include "patchWave.H"
#include "fvMesh.H"
#include "wallPolyPatch.H"
#include "fvPatchField.H"
#include "Field.H"
#include "emptyFvPatchFields.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(wallQuenchingDampingFunction, 0);
defineRunTimeSelectionTable(wallQuenchingDampingFunction, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void wallQuenchingDampingFunction::correctWallDist(const wordList& patches)
{
    // Get patchids of walls
    labelHashSet wallPatchIDs;
    if (patches.size()!=0)
    {
        Info<<"Quenching enabled on patches "<<patches<<endl;
        wallPatchIDs=getPatchIDs(patches); // get only selected patches
    }
    else
    {
        Info<<"Quenching enabled on all wall patches"<<endl;
        wallPatchIDs=getPatchIDs(wordList(1, wallPolyPatch::typeName)); // get all walls
    }
    
    // Calculate distance starting from wallPatch faces.
    patchWave wave(cellDistFuncs::mesh(), wallPatchIDs, true);
    
    // Transfer cell values from wave into *this 
    y_.transfer(wave.distance());
    
    // Transfer values on patches into boundaryField of *this
    forAll(y_.boundaryField(), patchI)
    {
        if 
            (
                y_.boundaryField()[patchI].type() 
                != emptyFvPatchScalarField::typeName
            )
        {
            scalarField& waveFld = wave.patchDistance()[patchI];
            
            y_.boundaryField()[patchI].transfer(waveFld);
        }
    }
    // Transfer number of unset values
    nUnset_ = wave.nUnset();

}

wallQuenchingDampingFunction::wallQuenchingDampingFunction
(
    const fvMesh& mesh,
    const dictionary& d
)
:
    cellDistFuncs(mesh),
    nUnset_(0),
    y_
    (
        IOobject
        (
            "y",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("y", dimLength, GREAT)
    )
{
    wordList patches;
    if (d.found("patches"))
        patches=wordList(d.lookup("patches"));
    correctWallDist(patches);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<wallQuenchingDampingFunction> wallQuenchingDampingFunction::New
(
    const fvMesh& mesh,
    const dictionary& d,
    word keyword
)
{
    word tn(d.lookup(keyword));

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(tn);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "wallQuenchingDampingFunction::New()"
        )   << "Unknown damping function type " << tn
            << endl << endl
            << "Valid damping function types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    Info<<"Loading wallQuenchingDampingFunction "<<tn<<endl;

    return autoPtr<wallQuenchingDampingFunction>(cstrIter()(mesh, d));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

wallQuenchingDampingFunction::~wallQuenchingDampingFunction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

}

// ************************************************************************* //
