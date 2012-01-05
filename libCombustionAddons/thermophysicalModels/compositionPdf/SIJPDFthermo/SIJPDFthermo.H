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
    SIJPDFthermo

Description
    Thermophysical properties integrated over statistical independent 
    joint PDFs

SourceFiles
    SIJPDFthermo.C

\*---------------------------------------------------------------------------*/

#ifndef SIJPDFTHERMO_H
#define SIJPDFTHERMO_H

#include "abstractProgressVariableThermo.H"
#include "autoPtr.H"
#include "IFstream.H"

#include "chemistryTable.H"
#include "pdf.H"
#include "wallQuenchingDampingFunction.H"

#include "graph.H"

#include "integrationScheduleEntry.H"

//#include "SIJPDFthermoStateFinder.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
        Class SIJPDFthermo Declaration
\*---------------------------------------------------------------------------*/





    template<class IndexDriver>
    class SIJPDFthermo
    : public abstractProgressVariableThermo
    {

        //static label elementsdone;

    protected:
        IOdictionary rawTableDict_;
        autoPtr<wallQuenchingDampingFunction> quench_;

   	chemistryTable rawTable_;

        PtrList<pdf> pdf_;
        HashTable<label> variable_;
        
        autoPtr<chemistryTable> preIntegratedTable_;

        List<scalar> constNormValues_;

        PtrList<volScalarField> meanWdot_;
        PtrList<volScalarField> meanWdotC_;

        multivariateSurfaceInterpolationScheme<scalar>::fieldTable mvfields_;


        // Private member functions
        
        sumOfIntegrableFunctions<scalar> multipliedPdfs(List<label>, label idx=0);

        void preIntegrate();

        void fillIntegrationSchedule //preIntegrate
        (
            List<label>&, 
            DynamicList<integrationScheduleEntry>& schedule,
            label varDim=0
        );

        virtual void boundProgressVariable(const word&);
        virtual void boundProgressVariable(label);
        virtual void boundProgressVariables();

        virtual void performIntegration
        (
            List<label>&,
            chemicalSystemState&,
            integrationParameter&,
            integrationParameter&
        );

        virtual HashTable<label,word> createPVDescription() const;
        virtual HashTable<label,word> createContentDescription() const;

        virtual void createPreIntegratedTable();

        virtual void initialize();


    public:

        //- Runtime type information
        TypeName("SIJPDFthermo");

        typedef IndexDriver ThermoIndexDriver;
      
        // Constructors
        
        //- Construct from dictionary and mesh
        SIJPDFthermo
        (
            const fvMesh&
        );
        
        
        // Destructor        
        virtual ~SIJPDFthermo();


        // Member functions


        // Access to thermodynamic state variables
        
        //- Enthalpy [J/kg]
        //  Non-const access allowed for transport equations
        volScalarField& h()
            {
                return h_();
            }
        
        //- Enthalpy [J/kg]
        const volScalarField& h() const
            {
                return h_();
            }
        
        // Fields derived from thermodynamic state variables
        
        //- Enthalpy for cell-set [J/kg]
        tmp<scalarField> h
        (
            const scalarField& T,
            const labelList& cells
        ) const;
        
        //- Enthalpy for patch [J/kg]
        tmp<scalarField> h
        (
            const scalarField& T,
            const label patchi
        ) const;
        
        
        //- Heat capacity at constant pressure for patch [J/kg/K]
        tmp<scalarField> Cp(const scalarField& T, const label patchi) const;
        
        //- Heat capacity at constant pressure [J/kg/K]
        tmp<volScalarField> Cp() const;
        
        /*
        //- Heat capacity at constant volume [J/kg/K]
        tmp<volScalarField> Cv() const;           
            */

        //- Density [kg/m^3]
        virtual tmp<volScalarField> rho() const
        {
	  return p()*psi();
        }

        
        scalar lookupCellByTemp
        (
            scalar,
            label,
            chemicalSystemState* state=NULL
        ) const;

        scalar lookupFaceByTemp
        (
            scalar, 
            label,
            label,
            chemicalSystemState* state=NULL
        ) const;

        chemicalSystemState lookupCell
        (
            label
        ) const;

        chemicalSystemState lookupFace
        (
            label, 
            label
        ) const;
 
            
        inline const chemistryTable& preIntegratedTable() const
            {
                return preIntegratedTable_();
            }

        tmp<volScalarField> normalizedFirstMoment(label I) const;

        inline const List<scalar>& constNormValues() const
            {
                return constNormValues_;
            }

        inline const PtrList<pdf>& presumedPdfs() const
            {
                return pdf_;
            }

        inline const HashTable<label>& variable() const
            {
                return variable_;
            }

	inline const volScalarField& meanWdot(const word& pv) const
	 {
	   return meanWdot_[variable_[pv]];
	 }
       
        //- Update properties
        virtual void correct();

        virtual void registerScalarFields
        (
            compressible::LESModel&
        );

        virtual void solveProgressVariables
        (
            const volScalarField&,
            const surfaceScalarField&,
            const compressible::LESModel&
        );
        
        virtual void solveProgressVariables
        (
            const volScalarField&,
            const surfaceScalarField&,
            const compressible::RASModel&
        );
        

        //- Read thermophysicalProperties dictionary
        virtual bool read();


        fileName cacheFile() const;
        bool readPreIntegratedTableFromCache();
        void writeCache() const;

    };

    #define forAllVariables(it)                                         \
                                                                        \
for (                                                                   \
    HashTable<label>::const_iterator it=variable_.begin();              \
    it!=variable_.end();                                                \
    it++                                                                \
)                                                                

     #define forAllVariablesIn(thermo, it)                           \
                                                                        \
for (                                                                   \
    HashTable<label>::const_iterator it=thermo.variable().begin();      \
    it!=thermo.variable().end();                                        \
    it++                                                                \
)                                                                       
    

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "SIJPDFthermo.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
