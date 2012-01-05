/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

#include "abstractProgressVariableThermo.H"
#include "fvMesh.H"
#include "OFstream.H"

#include "SIJPDFthermo.H"
#include "SIJPDFthermoStateFinder.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

#ifdef DEBUGFILE
word pnum()
{
OStringStream ss;
ss<<Pstream::myProcNo();
Pout<<Pstream::myProcNo();
return ss.str();
}

autoPtr<OFstream> dbgfile;

#endif


/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */
            
    defineTypeNameAndDebug(abstractProgressVariableThermo, 0);
    defineRunTimeSelectionTable(abstractProgressVariableThermo, fvMesh);

    typedef SIJPDFthermo<SIJPDFthermoIndexDriver> SIJPDFconstZthermo;
    defineTemplateTypeNameAndDebugWithName(SIJPDFconstZthermo, "SIJPDFconstZthermo", 0);
    
    addToRunTimeSelectionTable
    (
        abstractProgressVariableThermo,
        SIJPDFconstZthermo,
        fvMesh
    );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    autoPtr<abstractProgressVariableThermo> 
     abstractProgressVariableThermo::New
    (
        const fvMesh& mesh
    )
    {
        word thermoTypeName;

        // Enclose the creation of the thermophysicalProperties to ensure it is
        // deleted before the turbulenceModel is created otherwise the dictionary
        // is entered in the database twice
        {
            IOdictionary thermoDict
                (
                    IOobject
                    (
                        "thermophysicalProperties",
                        mesh.time().constant(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                );

            thermoDict.lookup("thermoType") >> thermoTypeName;
        }
        
        Info<< "Selecting thermodynamics package "
            << thermoTypeName << endl;

        fvMeshConstructorTable::iterator cstrIter =
            fvMeshConstructorTablePtr_->find(thermoTypeName);
        
        if (cstrIter == fvMeshConstructorTablePtr_->end())
        {
            FatalErrorIn("abstractProgressVariableThermo::New()")
                << "Unknown abstractProgressVariableThermo type "
                    << thermoTypeName << endl << endl
                    << "Valid abstractProgressVariableThermo types are :" << endl
                    << fvMeshConstructorTablePtr_->toc()
                    << exit(FatalError);
        }

        autoPtr<abstractProgressVariableThermo> newthermo(cstrIter()(mesh));
        newthermo->initialize();
        return newthermo;
    }
        
    // Constructors
        
    //- Construct from dictionary and mesh
    abstractProgressVariableThermo::abstractProgressVariableThermo
    (
        const fvMesh& mesh
    )
        : basicThermo(mesh)
    {
#ifdef DEBUGFILE
 dbgfile.reset(new OFstream(word("scatterplot.debug.")&pnum()));
#endif
    }
        
        
    // Destructor        
    abstractProgressVariableThermo::~abstractProgressVariableThermo()
    {
    }


} // End namespace Foam

// ************************************************************************* //
