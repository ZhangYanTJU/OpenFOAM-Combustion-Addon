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

\*----------------------------------------------------------------------------*/

#include "SIJPDFvarZthermo.H"
#include "SIJPDFvarZthermoStateFinder.H"
#include "fvMesh.H"
#include "OFstream.H"

#include "beta.H"
#include "addToRunTimeSelectionTable.H"
#include "integrableFunction.H"
#include "integrator.H"

namespace Foam
{
    typedef SIJPDFthermo<SIJPDFvarZthermoIndexDriver> dummy;
    defineTemplateTypeNameAndDebug(dummy, 0);

    defineTemplateTypeNameAndDebug(SIJPDFvarZthermo, 0);
    
    addToRunTimeSelectionTable
    (
        abstractProgressVariableThermo,
        SIJPDFvarZthermo,
        fvMesh
    );
        
// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void SIJPDFvarZthermo::performIntegration
(
    //int pIcursor[],
    List<label>& pIcursor,
    chemicalSystemState& integratedState,
    integrationParameter& ix,
    integrationParameter& constants
)
{
    // fill integrated state
    sumOfIntegrableFunctions<scalar> pdfs=multipliedPdfs(pIcursor);
      
    const ListDescription& cdesc=preIntegratedTable_().contentDescription();

    // source terms for moment transport equations
    //forAll(pdf_, I)
    forAllVariables(var)
    {
          
        integratedState[2*var()]=
            integrator<scalar>::integrate
            (
                rawTable_.entry("Ydot" & pdf_[var()].parameterName())
                * pdfs,
                ix, constants
            );

        integratedState[2*var()+1]=
            integrator<scalar>::integrate
            (
                singleParameter(pdf_[var()].parameterName())
                * (*normValues_[pdf_[var()].parameterName()])
                * rawTable_.entry("Ydot" & pdf_[var()].parameterName())
                * pdfs,
                ix, constants
            );
               
        word n="grad("&pdf_[var()].parameterName()&")";
        if (cdesc.found(n))
        {
            integratedState[cdesc[n]]=
                integrator<scalar>::integrate
                (
                    inv(*normValues_[pdf_[var()].parameterName()])
                    * rawTable_.entry(n) 
                    * pdfs,
                    ix, constants
                );
        }
    }

    label Irho=2*pdf_.size();

    // mean mixture properties
    integratedState[Irho]=
        integrator<scalar>::integrate
        (
            rawTable_.entry("rho") * pdfs,
            ix, constants
        );
    integratedState[Irho+1]=
        integrator<scalar>::integrate
        (
            rawTable_.entry("T") * pdfs,
            ix, constants
        );

    integratedState[Irho+2]=
        integrator<scalar>::integrate
        (
            rawTable_.entry("lambda") * pdfs,
            ix, constants
        );
    integratedState[Irho+3]=
        integrator<scalar>::integrate
        (
            rawTable_.entry("Cp") * pdfs,
            ix, constants
        );
    integratedState[Irho+4]=
        integrator<scalar>::integrate
        (
            rawTable_.entry("mu") * pdfs,
            ix, constants
        );

    //Info<<"I="<<integratedState<<endl;

}


void SIJPDFvarZthermo::initialize()
{

    // * * * * * * * * * * * * * * * * * * * * * * * * 
    //      - normalization quantities for pdf1:
    //        <(1/NPV)>, <(1/NPV)^2>, <(NPV'')^2>
    // * * * * * * * * * * * * * * * * * * * * * * * * 

    const label n[]=
        {
            pdf_[mixFracI_].tabResolutionFM(),
            pdf_[mixFracI_].tabResolutionSM()
        };
    const scalar s[]={0.0, 0.0};
    const scalar d[]={1.0/scalar(n[0]-1), 1.0/scalar(n[1]-1)};
    
    bool cacheMatches=readNormFactorsFromCache();
    
    if (cacheMatches)
    {
        // has required resolution?
        forAllVariables(var)
            cacheMatches=
            cacheMatches
            &&
            (derivedNormFactors_[var.key()]->nElements(0) == n[0])
            &&
            (derivedNormFactors_[var.key()]->nElements(1) == n[1]);
    }  

    if (!cacheMatches)
    {
        // build table
        derivedNormFactors_.clear();
        //derivedNormFactors_.setSize( variable_.size() );
        
        Info << "Creating tables of normalization factors" << endl;
        
        forAllVariables(var)
        {
            Info<<" Creating table for " << var.key() << endl;
            
            MultidimensionalLookupTable<scalarField>* nftab=
                new MultidimensionalLookupTable<scalarField>(2, n, s, d);

            pdf& mixFracPdf=pdf_[mixFracI_];
            
            int i[2];
            for (
                i[0] = 0;
                i[0] < nftab->nElements(0);
                i[0]++
            )
                for (
                  i[1] = 0;
                  i[1] < nftab->nElements(1);
                  i[1]++
                  )
            {

                scalar mean=nftab->valueAt(0, i[0]);
                scalar variance=nftab->valueAt(1, i[1]);

               mixFracPdf.setParameters
                    (
                        mean,
                        variance*mean*(1.0-mean)+mean*mean
                    );


                integrationParameter ix;
                ix.insert(mixFracName_, 0.0);
                //ix.insert(mixFracName_&"_Var", 0.0);

                integrationParameter constants;

                scalarField I(3, 0.0);

                I[0]=
                    integrator<scalar>::integrate
                    (
                        inv( *normValues_[var.key()] )
                        * mixFracPdf,
                        ix, constants
                    );

                I[1]=
                    integrator<scalar>::integrate
                    (
                        inv( *normValues_[var.key()] ) 
                        * inv( *normValues_[var.key()] )
                        * mixFracPdf,
                        ix, constants
                    );

                scalar normcmean=
                    integrator<scalar>::integrate
                    (
                        (*normValues_[var.key()])
                        * mixFracPdf,
                        ix, constants
                    );

                I[2]=
                    integrator<scalar>::integrate
                    (
                        ( 
                            (*normValues_[var.key()]) 
                            * 
                            (*normValues_[var.key()])
                        )
                        * mixFracPdf,
                        ix, constants
                    ) - sqr(normcmean);

                nftab->access(i)=I;

                Info<<mean<<" "<<variance<<":"
                    <<I[0]<<" "<<I[1]<<" "<<I[2]<<endl;
            }

            derivedNormFactors_.insert(var.key(), nftab);            
        }
      }

    SIJPDFthermo<SIJPDFvarZthermoIndexDriver>::initialize();

}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SIJPDFvarZthermo::SIJPDFvarZthermo(const fvMesh& mesh)
:
    SIJPDFthermo<SIJPDFvarZthermoIndexDriver>(mesh),
    mixFracName_(word(lookup("mixtureFractionVariableName"))),
    mixFracI_(-1)
{

    // correct list of mixture fraction dependent variables
    mixFracI_=variable_[mixFracName_];
    variable_.erase(mixFracName_);

    if (variable_.size()!=pdf_.size()-1)
    {
      FatalErrorIn("SIJPDFvarZthermo::SIJPDFvarZthermo()")
          << " mixture fraction declared as variable (name " 
              << mixFracName_ <<")" << endl
              << " but either none PDF or multiple PDFs given"
              <<" for the mixture fraction." << endl
              <<  variable_
              << exit(FatalError);
    }

    // read graphs of normalization factors
    forAllVariables(var)
    {
        graph g(rawTableDict_.lookup("normValues_" & var.key()));

        normValues_.insert
            (
                var.key(),
                new integrableGraph(mixFracName_, g.x(), g.y())
            );
    }
}




bool SIJPDFvarZthermo::readNormFactorsFromCache()
{
  const Time& time=psi_.mesh().time();
  fileName fullPath=time.path()/time.constant()/cacheFile();
  if (exists(fullPath))
  {
    IOdictionary tableDict
    (
      IOobject
      (
        cacheFile(),
        psi_.mesh().time().constant(),
        psi_.mesh().time(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
      )
    );

    if (tableDict.found("derivedNormFactors"))
    {

        Info << "Reading tables of normalization factors from" 
            << fullPath << endl;

        tableDict.lookup("derivedNormFactors") >> derivedNormFactors_;

    } else return false;

    return true;
  }
  else return false;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SIJPDFvarZthermo::~SIJPDFvarZthermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

}

// ************************************************************************* //
