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
    timeVaryingUniformFixedValueFvPatchField

Description

SourceFiles
    timeVaryingUniformFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef freeFlameInletFvPatchVectorField_H
#define freeFlameInletFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "cuttingPlane.H"
#include "patchToPatchInterpolation.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class timeVaryingUniformFixedValueFvPatch Declaration
\*---------------------------------------------------------------------------*/

class freeFlameInletFvPatchVectorField
:
    public fixedValueFvPatchField<vector>
{
    // Private data

    scalar distanceFromPatch_;
    vector patchNormal_;
    autoPtr<cuttingPlane> sourcePlane_;
    autoPtr<primitiveFacePatch> planePatch_;
    autoPtr<PatchToPatchInterpolation<primitiveFacePatch, primitivePatch> > interPatch_;

    void createInterpolationObjects(const scalar&);

public:

    //- Runtime type information
    TypeName("freeFlameInlet");


    // Constructors

        //- Construct from patch and internal field
        freeFlameInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        freeFlameInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given freeFlameInletFvPatchVectorField
        //  onto a new patch
        freeFlameInletFvPatchVectorField
        (
            const freeFlameInletFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        freeFlameInletFvPatchVectorField
        (
            const freeFlameInletFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new freeFlameInletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        freeFlameInletFvPatchVectorField
        (
            const freeFlameInletFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchField<vector> > clone(const DimensionedField<vector, volMesh>& iF) const
        {
            return tmp<fvPatchField<vector> >
            (
                new freeFlameInletFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
