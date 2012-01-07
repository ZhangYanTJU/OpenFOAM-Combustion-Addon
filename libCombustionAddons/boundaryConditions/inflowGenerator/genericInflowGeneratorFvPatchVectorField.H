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
    genericInflowGeneratorFvPatchVectorField

Description

SourceFiles
    genericInflowGeneratorFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef genericInflowGeneratorFvPatchVectorField_H
#define genericInflowGeneratorFvPatchVectorField_H

#include "Random.H"
#include "fixedValueFvPatchFields.H"
#include "homogeneousTurbulence.H"
#include "hat.H"
#include "SLList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class genericInflowGeneratorFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class genericInflowGeneratorFvPatchVectorField
:
    public fixedValueFvPatchField<vector>
{
    
    // Private data

protected:
    
    Random ranGen_;
    vector fluctuationScale_;
    Field<vector> referenceField_;
    Field<symmTensor> RField_;
    label curTimeIndex_;

public:

    //- Runtime type information
    TypeName("genericInflowGenerator");


    // Constructors

        //- Construct from patch and internal field
        genericInflowGeneratorFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        genericInflowGeneratorFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given genericInflowGeneratorFvPatchVectorField
        //  onto a new patch
        genericInflowGeneratorFvPatchVectorField
        (
            const genericInflowGeneratorFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        genericInflowGeneratorFvPatchVectorField
        (
            const genericInflowGeneratorFvPatchVectorField&
        );

        //- Construct as copy setting internal field reference
        genericInflowGeneratorFvPatchVectorField
        (
            const genericInflowGeneratorFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );



    // Member functions

        // Access

            //- Return the fluctuation scale
            const vector& fluctuationScale() const
            {
                return fluctuationScale_;
            }

            //- Return reference to the fluctuation scale to allow adjustment
            vector& fluctuationScale()
            {
                return fluctuationScale_;
            }

            //- Return the reference field
            const Field<vector>& referenceField() const
            {
                return referenceField_;
            }

            //- Return reference to the reference field to allow adjustment
            Field<vector>& referenceField()
            {
                return referenceField_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<vector>&,
                const labelList&
            );


        // Evaluation functions

            virtual void doUpdate() =0;

            virtual void updateCoeffs();



        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
