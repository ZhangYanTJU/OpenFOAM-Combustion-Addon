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

#include "genericInflowGeneratorFvPatchVectorField.H"
#include "transform.H"
#include "transformField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
defineTypeNameAndDebug(genericInflowGeneratorFvPatchVectorField, 0);

genericInflowGeneratorFvPatchVectorField::genericInflowGeneratorFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    ranGen_(label(0)),
    fluctuationScale_(pTraits<vector>::one),
    referenceField_(p.size()),
    RField_(p.size()),
    curTimeIndex_(-1)
{
}

genericInflowGeneratorFvPatchVectorField::genericInflowGeneratorFvPatchVectorField
(
    const genericInflowGeneratorFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    ranGen_(label(0)),
    fluctuationScale_(ptf.fluctuationScale_),
    referenceField_(ptf.referenceField_, mapper),
    RField_(ptf.RField_, mapper),
    curTimeIndex_(-1)
{
}

genericInflowGeneratorFvPatchVectorField::genericInflowGeneratorFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF),
    ranGen_(label(0)),
    fluctuationScale_(dict.lookup("fluctuationScale")),
    referenceField_("referenceField", dict, p.size()),
    RField_(p.size(), I),
    curTimeIndex_(-1)
{
    if (dict.found("R"))
     RField_=symmTensorField("R", dict, p.size());
    if (dict.found("value"))
    {
        fixedValueFvPatchField<vector>::operator==
        (
            Field<vector>("value", dict, p.size())
        );
    }
    else
    {
        fixedValueFvPatchField<vector>::operator==(referenceField_);
    }
}

genericInflowGeneratorFvPatchVectorField::genericInflowGeneratorFvPatchVectorField
(
    const genericInflowGeneratorFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    ranGen_(ptf.ranGen_),
    fluctuationScale_(ptf.fluctuationScale_),
    referenceField_(ptf.referenceField_),
    RField_(ptf.RField_),
    curTimeIndex_(-1)
{
}

genericInflowGeneratorFvPatchVectorField::genericInflowGeneratorFvPatchVectorField
(
    const genericInflowGeneratorFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    ranGen_(ptf.ranGen_),
    fluctuationScale_(ptf.fluctuationScale_),
    referenceField_(ptf.referenceField_),
    RField_(ptf.RField_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void genericInflowGeneratorFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    Field<vector>::autoMap(m);
    referenceField_.autoMap(m);
    RField_.autoMap(m);
}


void genericInflowGeneratorFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<vector>::rmap(ptf, addr);

    const genericInflowGeneratorFvPatchVectorField& tiptf =
        refCast<const genericInflowGeneratorFvPatchVectorField >(ptf);

    referenceField_.rmap(tiptf.referenceField_, addr);
    RField_.rmap(tiptf.RField_, addr);
}


void genericInflowGeneratorFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {

        doUpdate();

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchField<vector>::updateCoeffs();
}



void genericInflowGeneratorFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeKeyword("fluctuationScale")
        << fluctuationScale_ << token::END_STATEMENT << nl;
    referenceField_.writeEntry("referenceField", os);
    RField_.writeEntry("R", os);
    this->writeEntry("value", os);
}


} // End namespace Foam

// ************************************************************************* //
