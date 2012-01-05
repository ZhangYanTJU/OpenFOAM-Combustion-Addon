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
    oscillatingFixedValueFvPatchField

Description

SourceFiles
    oscillatingFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef ORACLESinflowFvPatchVectorField_H
#define ORACLESinflowFvPatchVectorField_H

#include "Random.H"
#include "inflowGeneratorFvPatchVectorField.H"
#include "graph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class ORACLESinflowFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class ORACLESinflowFvPatchVectorField
:
    public inflowGeneratorFvPatchVectorField<homogeneousTurbulence>
{
    // Private data

    struct xy
    {
        scalar x_, y_;
        
        xy()
            {}
        
        // Friend Operators
        
        friend bool operator==(const xy& a, const xy& b)
            {
                return equal(a.x_, b.x_) && equal(a.y_, b.y_);
            }
        
        friend bool operator!=(const xy& a, const xy& b)
            {
                return !(a == b);
            }
        
        friend Istream& operator>>(Istream& is, xy& xyd)
            {
                is >> xyd.x_ >> xyd.y_;
                return is;
            }
        
        friend Ostream& operator<<(Ostream& os, const xy& xyd)
            {
                os << xyd.x_ << ' ' << xyd.y_;
                return os;
            }
    };
    

    //- Reference value
    scalarField y_;
    scalarField UbyU0_;
    scalar U0_;

        //- Amplitude
        scalar amplitude_;

        //- Frequency
        scalar frequency_;

        //- Current time index
        label curTimeIndex_;


    // Private member functions

        //- Return current scale
        scalar currentScale() const;

    tmp<vectorField> refFieldFromProfile();


public:

    //- Runtime type information
    TypeName("ORACLESinflow");


    // Constructors

        //- Construct from patch and internal field
        ORACLESinflowFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        ORACLESinflowFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given ORACLESinflowFvPatchVectorField
        //  onto a new patch
        ORACLESinflowFvPatchVectorField
        (
            const ORACLESinflowFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        ORACLESinflowFvPatchVectorField
        (
            const ORACLESinflowFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new ORACLESinflowFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        ORACLESinflowFvPatchVectorField
        (
            const ORACLESinflowFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new ORACLESinflowFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return amplitude
            scalar amplitude() const
            {
                return amplitude_;
            }

            scalar& amplitude()
            {
                return amplitude_;
            }

            //- Return frequency
            scalar frequency() const
            {
                return frequency_;
            }

            scalar& frequency()
            {
                return frequency_;
            }



        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
