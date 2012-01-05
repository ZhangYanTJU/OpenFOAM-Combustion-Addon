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

\*---------------------------------------------------------------------------*/
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "freeFlameInlet.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

freeFlameInletFvPatchVectorField::
freeFlameInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    distanceFromPatch_(0.0),    
    patchNormal_(0,0,0)
{}


freeFlameInletFvPatchVectorField::
freeFlameInletFvPatchVectorField
(
    const freeFlameInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    distanceFromPatch_(ptf.distanceFromPatch_),
    patchNormal_(0,0,0)
{
    createInterpolationObjects(distanceFromPatch_);
}


freeFlameInletFvPatchVectorField::
freeFlameInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF),
    distanceFromPatch_(readScalar(dict.lookup("patchNormalDistance"))),
    patchNormal_(0,0,0)
{
    createInterpolationObjects(distanceFromPatch_);

    if (dict.found("value"))
    {
        fvPatchField<vector>::operator==(Field<vector>("value", dict, p.size()));
    }
    else
    {
        updateCoeffs();
    }
}


freeFlameInletFvPatchVectorField::
freeFlameInletFvPatchVectorField
(
    const freeFlameInletFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    distanceFromPatch_(ptf.distanceFromPatch_),
    patchNormal_(0,0,0)
{
    createInterpolationObjects(distanceFromPatch_);
}


freeFlameInletFvPatchVectorField::
freeFlameInletFvPatchVectorField
(
    const freeFlameInletFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    distanceFromPatch_(ptf.distanceFromPatch_),
    patchNormal_(0,0,0)
{
    createInterpolationObjects(distanceFromPatch_);
}


void freeFlameInletFvPatchVectorField::createInterpolationObjects
(
    const scalar& dist
)
{
    if (size()>0)
    {
        patchNormal_=-average(this->patch().nf());
        sourcePlane_.reset
            (
                new cuttingPlane
                (
                    this->patch().boundaryMesh().mesh(),
                    plane
                    (
                        this->patch().Cf()[0]+patchNormal_*dist,
                        patchNormal_
                    )
                )
            );
        planePatch_.reset
            (
                new primitiveFacePatch
                (                    
                    sourcePlane_().faces(),
                    sourcePlane_().points()
                )
            );
        interPatch_.reset
            (
                new PatchToPatchInterpolation<primitiveFacePatch, primitivePatch>
                (
                    planePatch_(),
                    this->patch().patch()
                )
            );
    }
}
            
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
// Update the coefficients associated with the patch field

void freeFlameInletFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (size()>0)
    {

        if (!sourcePlane_().cut())
            FatalErrorIn("freeFlameInletFvPatchVectorField::updateCoeffs()")
                << "could not intersect cutting plane with mesh" << endl
                    << " distance=" << distanceFromPatch_ << endl           
                    << " in direction " << patchNormal_ << endl
                    << abort(FatalError);
        
        const volVectorField& U =
            patch().boundaryMesh().mesh().objectRegistry::lookupObject<volVectorField>("U");
        const volScalarField& rho =
            patch().boundaryMesh().mesh().objectRegistry::lookupObject<volScalarField>("rho");
       
	vectorField Uc=interPatch_().faceInterpolate(sourcePlane_().sample(U)); 
	scalarField rhoc=interPatch_().faceInterpolate(sourcePlane_().sample(rho)); 
      
        scalar flux=gSum(rhoc*(patch().Sf() & Uc));
        scalar thisflux = gSum(rho.boundaryField()[patch().index()]*(patch().Sf() & (*this)));
        
        scalar cfac= (thisflux!=0.0)? flux/thisflux : 1.0;
        
        Pout << "Correction factor: " << cfac
            << ", Uncorrected flux: "
            << thisflux
            << ", Target flux: "
            << flux
            << endl;
        
        operator==(cfac*(*this));
    }

    fixedValueFvPatchField<vector>::updateCoeffs();
}


// Write
void freeFlameInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeKeyword("patchNormalDistance")
        << distanceFromPatch_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


makePatchTypeField
(
    fvPatchVectorField,
    freeFlameInletFvPatchVectorField
);


} // End namespace Foam

// ************************************************************************* //
