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

\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"
#include "ORACLESinflowFvPatchVectorField.H"
#include "mathematicalConstants.H"
#include "volMesh.H"
#include "interpolateXY.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

scalar ORACLESinflowFvPatchVectorField::currentScale() const
{
    return
        1.0
      + amplitude_*
        sin(2*mathematicalConstant::pi*frequency_*this->db().time().value());
}

tmp<vectorField> ORACLESinflowFvPatchVectorField::refFieldFromProfile()
{
    tmp<vectorField> tref(new vectorField(this->size()));
    vectorField& ref=tref();

    ref=
        vector(U0_,0,0)
        *
        interpolateXY
        (
            patch().Cf().component(vector::Y),
            y_, UbyU0_
        );

    return tref;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ORACLESinflowFvPatchVectorField::ORACLESinflowFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    inflowGeneratorFvPatchVectorField<homogeneousTurbulence>(p, iF),
    y_(0),
    UbyU0_(0),
    U0_(0.0),
    amplitude_(0.0),
    frequency_(0.0),
    curTimeIndex_(-1)
{}


ORACLESinflowFvPatchVectorField::ORACLESinflowFvPatchVectorField
(
    const ORACLESinflowFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inflowGeneratorFvPatchVectorField<homogeneousTurbulence>(ptf, p, iF, mapper),
    y_(ptf.y_),
    UbyU0_(ptf.UbyU0_),
    U0_(ptf.U0_),
    amplitude_(ptf.amplitude_),
    frequency_(ptf.frequency_),
    curTimeIndex_(-1)
{}


ORACLESinflowFvPatchVectorField::ORACLESinflowFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    inflowGeneratorFvPatchVectorField<homogeneousTurbulence>(p, iF, dict),
    y_(0),
    UbyU0_(0),
    U0_(readScalar(dict.lookup("U0"))),
    amplitude_(readScalar(dict.lookup("amplitude"))),
    frequency_(readScalar(dict.lookup("frequency"))),
    curTimeIndex_(-1)
{
    List<xy> xyData(dict.lookup("profile"));
    
    y_.setSize(xyData.size());
    UbyU0_.setSize(xyData.size());

    forAll (xyData, i)
    {
        y_[i] = xyData[i].x_;
        UbyU0_[i] = xyData[i].y_;
    }

    this->referenceField_=refFieldFromProfile();

}


ORACLESinflowFvPatchVectorField::ORACLESinflowFvPatchVectorField
(
    const ORACLESinflowFvPatchVectorField& ptf
)
:
    inflowGeneratorFvPatchVectorField<homogeneousTurbulence>(ptf),
    y_(ptf.y_),
    UbyU0_(ptf.UbyU0_),
    U0_(ptf.U0_),
    amplitude_(ptf.amplitude_),
    frequency_(ptf.frequency_),
    curTimeIndex_(-1)
{}


ORACLESinflowFvPatchVectorField::ORACLESinflowFvPatchVectorField
(
    const ORACLESinflowFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    inflowGeneratorFvPatchVectorField<homogeneousTurbulence>(ptf, iF),
    y_(ptf.y_),
    UbyU0_(ptf.UbyU0_),
    U0_(ptf.U0_),
    amplitude_(ptf.amplitude_),
    frequency_(ptf.frequency_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void ORACLESinflowFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        referenceField_ = refFieldFromProfile()*currentScale();

        curTimeIndex_ = this->db().time().timeIndex();
    }

    inflowGeneratorFvPatchVectorField<homogeneousTurbulence>::updateCoeffs();
}


void ORACLESinflowFvPatchVectorField::write(Ostream& os) const
{
    inflowGeneratorFvPatchVectorField<homogeneousTurbulence>::write(os);
    os.writeKeyword("U0")
        << U0_ << token::END_STATEMENT << nl;
    os.writeKeyword("profile")<<nl<<token::BEGIN_LIST<<nl;
    forAll(y_, i)
        os<<y_[i]<<token::SPACE<<UbyU0_[i]<<nl;
    os << token::END_LIST << token::END_STATEMENT << nl;
    os.writeKeyword("amplitude")
        << amplitude_ << token::END_STATEMENT << nl;
    os.writeKeyword("frequency")
        << frequency_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, ORACLESinflowFvPatchVectorField);

} // End namespace Foam

// ************************************************************************* //
