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

#include "combustionProgressFunctionObject.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(combustionProgressFunctionObject, 0);
    addToRunTimeSelectionTable(functionObject, combustionProgressFunctionObject, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
    /*
combustionProgressFunctionObject::combustionProgressFunctionObject()
{}
        */

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

combustionProgressFunctionObject::combustionProgressFunctionObject
(
    const word& n,
    const Time& t,
    const dictionary& functionDict
):
    name_(n),
    time_(t),
    mesh_
    (
        refCast<const fvMesh>
        (
            time_.lookupObject<objectRegistry>(polyMesh::defaultRegion)
        )
    ),
    T_(mesh_.lookupObject<volScalarField>("T")),
    Tthreshold_(functionDict.lookup("Tthreshold")),
    cubeLength_(functionDict.lookup("cubeLength"))
 
{
    fileName locfiledir;
    if (Pstream::parRun())
        locfiledir=time_.path()/".."/name_/time_.timeName();
    else
        locfiledir=time_.path()/name_/time_.timeName();

    mkDir(locfiledir);
    f_.reset(new OFstream(locfiledir/"combustionProgress.t"));

    //read(functionDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

combustionProgressFunctionObject::~combustionProgressFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- start is called at the start of the time-loop
bool combustionProgressFunctionObject::start()
{
    return true;
}
    
//- execute is called at each ++ or += of the time-loop
bool combustionProgressFunctionObject::execute()
{

    scalar maxT=max(T_).value();
    scalar progress=pos(T_ - Tthreshold_)().weightedAverage(mesh_.V()).value();
    vector COG=vector::zero;
    if (progress>SMALL) COG=gSum(mesh_.C().internalField()*(pos(T_ - Tthreshold_)*mesh_.V()))/
     gSum((pos(T_ - Tthreshold_)*mesh_.V()));
    scalar R=pow(3.0*pow(cubeLength_.value(), 3)*progress/(4.0*3.14), 1./3.);

    if (Pstream::master())
    {
     Info<<"Combustion Progress = "<<100.0*progress<<" (R="<<R<<")"<<endl;
     f_() << time_.timeName() << " " << R << " " << progress << " " << maxT <<" "
      <<COG.x()<<" "<<COG.y()<<" "<<COG.z()<<" "
      <<endl;
    }

    return true;
}
    
 
//- Read and set the function object if it's data has changed
bool combustionProgressFunctionObject::read(const dictionary& dict)
{
    Tthreshold_=dimensionedScalar(dict.lookup("Tthreshold"));
    cubeLength_=dimensionedScalar(dict.lookup("cubeLength"));
    Pout<<"Tthreshold="<<Tthreshold_<<", cuebLength="<<cubeLength_<<endl;
    return true;
}

}

// ************************************************************************* //
