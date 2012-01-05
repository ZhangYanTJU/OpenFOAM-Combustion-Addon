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

Class

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef ABSTRACTPROGRESSVARIABLETHERMO_H
#define ABSTRACTPROGRESSVARIABLETHERMO_H

#include "basicThermo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "IFstream.H"
#include "OFstream.H"
#include "compressible/LES/LESModel/LESModel.H"
#include "compressible/RAS/RASModel/RASModel.H"
#include "secondMoment.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#define DEBUGFILE
#undef DEBUGFILE

namespace Foam
{
#ifdef DEBUGFILE
    extern autoPtr<OFstream> dbgfile;
#endif

    class abstractProgressVariableThermo
    : public basicThermo
    {


    protected:

        virtual void initialize() =0;

        PtrList<volScalarField> firstMoment_;
        //PtrList<volScalarField> secondMoment_;
        PtrList<secondMomentSolver> secondMoment_;

        autoPtr<volScalarField> h_;

    public:
        //- Runtime type information
        TypeName("progressVariableThermo");

        // Declare run-time constructor selection tables
        declareRunTimeSelectionTable
        ( 
            autoPtr,
            abstractProgressVariableThermo,
            fvMesh,
            (
                const fvMesh& mesh
            ),
            (mesh)
        );

      
        static autoPtr<abstractProgressVariableThermo> New
        (
            const fvMesh& mesh
        );
        
        // Constructors
        
        //- Construct from dictionary and mesh
        abstractProgressVariableThermo
        (
            const fvMesh&
        );
        
        
        // Destructor        
        virtual ~abstractProgressVariableThermo();


        // Member functions    

        inline volScalarField& firstMoment(label i)
            {
                return firstMoment_[i];
            }

        inline secondMomentSolver& secondMoment(label i)
            {
                return secondMoment_[i];
            }

        inline const volScalarField& firstMoment(label i) const
            {
                return firstMoment_[i];
            }

        inline const secondMomentSolver& secondMoment(label i) const
            {
                return secondMoment_[i];
            }

        inline volScalarField& T()
            {
                return T_;
            }

        inline bool solvesEnthalpy() const
            {
                return h_.valid();
            }

        virtual void registerScalarFields
        (
            compressible::LESModel&
        ) =0;

        // Solve transport equations of progress variables (LES)
        virtual void solveProgressVariables
        (
            const volScalarField& rho, 
            const surfaceScalarField& phi,
            const compressible::LESModel&
        ) =0;
       
        // Solve transport equations of progress variables (RAS)
        virtual void solveProgressVariables
        (
            const volScalarField& rho, 
            const surfaceScalarField& phi,
            const compressible::RASModel&
        ) =0;
       
        virtual void boundProgressVariables() =0;
  
    };

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
