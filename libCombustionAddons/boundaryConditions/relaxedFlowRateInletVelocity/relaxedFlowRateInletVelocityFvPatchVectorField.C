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

#include "relaxedFlowRateInletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
relaxedFlowRateInletVelocityFvPatchVectorField::
relaxedFlowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    flowRate_(0),
    phiName_("phi"),
    rhoName_("rho"),
    relaxUp_(1.0),
    relaxDown_(1.0),
    lastUpdate_(db().time().startTime().value())
{}


Foam::
relaxedFlowRateInletVelocityFvPatchVectorField::
relaxedFlowRateInletVelocityFvPatchVectorField
(
    const relaxedFlowRateInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    flowRate_(ptf.flowRate_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    relaxUp_(ptf.relaxUp_),
    relaxDown_(ptf.relaxDown_),
    lastUpdate_(ptf.lastUpdate_)
{}


Foam::
relaxedFlowRateInletVelocityFvPatchVectorField::
relaxedFlowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    flowRate_(readScalar(dict.lookup("flowRate"))),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    relaxUp_(dict.lookupOrDefault<scalar>("relaxUp", 0.001)),
    relaxDown_(dict.lookupOrDefault<scalar>("relaxDown", 1.)),
    lastUpdate_(db().time().startTime().value())
{}


Foam::
relaxedFlowRateInletVelocityFvPatchVectorField::
relaxedFlowRateInletVelocityFvPatchVectorField
(
    const relaxedFlowRateInletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    flowRate_(ptf.flowRate_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    relaxUp_(ptf.relaxUp_),
    relaxDown_(ptf.relaxDown_),
    lastUpdate_(ptf.lastUpdate_)
{}


Foam::
relaxedFlowRateInletVelocityFvPatchVectorField::
relaxedFlowRateInletVelocityFvPatchVectorField
(
    const relaxedFlowRateInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    flowRate_(ptf.flowRate_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    relaxUp_(ptf.relaxUp_),
    relaxDown_(ptf.relaxDown_),
    lastUpdate_(ptf.lastUpdate_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::relaxedFlowRateInletVelocityFvPatchVectorField::updateCoeffs()
{
  if (updated() || (db().time().value()<=lastUpdate_))
    {
        return;
    }
  lastUpdate_=db().time().value();

    // A simpler way of doing this would be nice
    scalar avgU = -flowRate_/gSum(patch().magSf());

    vectorField n = patch().nf();

    scalarField newU(n.size());

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        // Volumetric flow-rate
      newU=avgU;
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        // Mass flow-rate
        newU=avgU/rhop;
    }
    else
    {
        FatalErrorIn
        (
            "relaxedFlowRateInletVelocityFvPatchVectorField::updateCoeffs()"
        )   << "dimensions of " << phiName_ << " are incorrect" << nl
            << "    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << nl << exit(FatalError);
    }
    //Info<<newU<<((*this)&n)<<endl;
    scalarField relax( n.size() );
    relax = (newU<((*this)&n) ? relaxUp_ : relaxDown_ );
    //Info<<relax<<endl;
    operator==( relax*n*newU + (1.-relax)*(*this) );

    Info<<"Setting new velocity to "<<(gAverage(*this))<<endl;

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::relaxedFlowRateInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeKeyword("flowRate") << flowRate_ << token::END_STATEMENT << nl;
    os.writeKeyword("relaxUp") << relaxUp_ << token::END_STATEMENT << nl;
    os.writeKeyword("relaxDown") << relaxDown_ << token::END_STATEMENT << nl;
    if (phiName_ != "phi")
    {
        os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    }
    if (rhoName_ != "rho")
    {
        os.writeKeyword("rho") << rhoName_ << token::END_STATEMENT << nl;
    }
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       relaxedFlowRateInletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
