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

#include "SIJPDFthermo.H"
#include "fvMesh.H"
#include "OFstream.H"

#include "beta.H"
#include "addToRunTimeSelectionTable.H"
#include "integrableFunction.H"
#include "integrator.H"

#include "ListListOps.H"
#include "PstreamReduceOps.H"



//#define COMBUSTAPPROX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    /*
    defineTemplateTypeNameAndDebug(SIJPDFthermo, 0);
    
    addToRunTimeSelectionTable
    (
        abstractProgressVariableThermo,
        SIJPDFthermo,
        fvMesh
    );
        */

template<class IndexDriver>
scalar SIJPDFthermo<IndexDriver>::lookupCellByTemp
(
    scalar temperature, 
    label cellI,
    chemicalSystemState* state
) const
{
    scalar h;

    string error;
    chemicalSystemState s=
        preIntegratedTable().reverseLookup
        (
            preIntegratedTable().indexOfPV("h"),
            IndexDriver::cellTableIndex(this, cellI),
            preIntegratedTable().indexOf("T"),
            temperature,
            h,
            error
        );

    if (state!=NULL) *state=s;
    return h;
}

template<class IndexDriver>
scalar SIJPDFthermo<IndexDriver>::lookupFaceByTemp
(
    scalar temperature,
    label patchI,
    label faceI,
    chemicalSystemState* state
) const
{
    scalar h;
        
    string error;
    chemicalSystemState s=
        preIntegratedTable().reverseLookup
        (
            preIntegratedTable().indexOfPV("h"),
            IndexDriver::faceTableIndex(this, patchI, faceI),
            preIntegratedTable().indexOf("T"),
            temperature,
            h,
            error
        );

    if (state!=NULL) *state=s;
    return h;
}

template<class IndexDriver>
chemicalSystemState SIJPDFthermo<IndexDriver>::lookupCell
(
    label cellI
) const
{
    return preIntegratedTable().lookup(IndexDriver::cellTableIndex(this, cellI));
}

template<class IndexDriver>
chemicalSystemState SIJPDFthermo<IndexDriver>::lookupFace
(
    label patchI,
    label faceI
) const
{
    return preIntegratedTable().lookup
        (
            IndexDriver::faceTableIndex(this, patchI, faceI)
        );
}

template<class IndexDriver>
tmp<volScalarField> SIJPDFthermo<IndexDriver>::normalizedFirstMoment(label I) const
{
    tmp<volScalarField> tnormM1
        (
            new volScalarField
            (
                IOobject
                (
                    pdf_[I].parameterName()&"_norm",
                    psi_.mesh().time().timeName(),
                    psi_.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                psi_.mesh(),
                dimensionedScalar("", dimless, 0.0)
            )
        );
    volScalarField& normM1=tnormM1();

    forAll(normM1, cellI)
        normM1[cellI]=IndexDriver::cellNormM1(this, cellI, I);

    forAll(normM1.boundaryField(), pI)
    {
        forAll(normM1.boundaryField()[pI], faceI)
        {
            normM1.boundaryField()[pI][faceI]=
                IndexDriver::faceNormM1(this, pI, faceI, I);
        }
    }

    return tnormM1;
}



template<class IndexDriver>
sumOfIntegrableFunctions<scalar> SIJPDFthermo<IndexDriver>::multipliedPdfs
(
  List<label> piCursor,
  label idx
)
{
  scalar mean=preIntegratedTable_().valueAt(2*idx, piCursor[2*idx]);
  scalar var=preIntegratedTable_().valueAt(2*idx+1, piCursor[2*idx+1]);
  
  //Info<<"setting "<<mean<<" "<<var<<endl;

  pdf_[idx].setParameters
  (
   mean,
   var*mean*(1.0-mean)+mean*mean
  );
  
  if (idx < pdf_.size()-1)
      return
          pdf_[idx]
          * multipliedPdfs(piCursor, idx+1);
  else
      return
          pdf_[idx];

}


template<class IndexDriver>
void SIJPDFthermo<IndexDriver>::preIntegrate()
{
    Info<<"Building integration schedule."<<endl;

    List<List<integrationScheduleEntry> > allSchedules(Pstream::nProcs());
    List<integrationScheduleEntry>& mySchedule=allSchedules[Pstream::myProcNo()];
    DynamicList<integrationScheduleEntry> schedule;

    label elementstotal=1;
    for (label i=0;i<preIntegratedTable_().dimensionality();i++)
    {
        elementstotal*=preIntegratedTable_().nElements(i);
    }

    if (elementstotal%Pstream::nProcs()!=0)
    {
        FatalErrorIn("SIJPDFthermo::preIntegrate()")
            <<"Total number of lookup table entries ("<<elementstotal<<") "
                <<"must be divideable by the number of processors ("<<Pstream::nProcs()<<"."
                <<abort(FatalError);
    }

    List<label> sizes(Pstream::nProcs(), 0);
    label sum=0;
    forAll(sizes, I)
    {
        sizes[I]=elementstotal/Pstream::nProcs();
        sum+=sizes[I];
    }
    //sizes[0] -= sum-elementstotal;
    //Info<<sizes<<endl;

    mySchedule.setSize(sizes[Pstream::myProcNo()]);

    List<label> icursor(preIntegratedTable_().dimensionality(), 0);
    fillIntegrationSchedule(icursor, schedule);
    schedule.shrink();
        
    Info<<"Distributing integration to processors"<<endl;
    label ofs=0;
    for (label p=0; p<Pstream::myProcNo(); p++)
        ofs+=sizes[p];

    for (label k=0; k<sizes[Pstream::myProcNo()]; k++)
    {
        if (ofs+k<schedule.size())
            mySchedule[k]=schedule[ofs+k];
    }

    Pstream::gatherList(allSchedules);
    Pstream::scatterList(allSchedules);

    //Pout<<"msize="<<mySchedule.size()<<endl;
        
    Info<<"Performing integration"<<endl;
    label integrated=0;
    label lastRep=0;
    Info<<"Progress: "<<flush;
    forAll(mySchedule, I)
    {
        integrationScheduleEntry& e = mySchedule[I];

        integrated++;
        label sumi=integrated;
        reduce(sumi, sumOp<label>());
        label percent=(100*sumi/elementstotal);
        if (percent>lastRep)
        {
            lastRep=percent;
            if (percent%10==0)
                Info<<percent<<"%"<<flush;
            else
                Info<<"."<<flush;
        }
        
        performIntegration(e.piCursor, e.integratedState, e.ix, e.constants);
    }
    Info<<" (Finished)"<<endl;

    Pstream::gatherList(allSchedules);
    Pstream::scatterList(allSchedules);

    Info<<"Combining table"<<endl;

    List<integrationScheduleEntry>  completeSchedule=
        ListListOps::combine<List<integrationScheduleEntry> >
        (
            allSchedules,
            accessOp<List<integrationScheduleEntry> >()
        );
    /*
    // combine entries
    if (Pstream::master())
    {
        Pout<<"master csize="<<combinedSchedule.size()<<endl;

        forAll(allSchedules, pI)
            allSchedules[pI]=combinedSchedule;

    }
        */
    // all processors get complete schedule
    //Pstream::scatterList(allSchedules);
    //Pout<<"msize="<<completeSchedule.size()<<endl;

    Info<<"Filling table"<<endl;
    // fill table on each processor
    forAll(completeSchedule, sI)
    {
        integrationScheduleEntry& e=completeSchedule[sI];

        // save integrated state to preintegrated table
        int ic[e.piCursor.size()];
        forAll(e.piCursor, I) ic[I]=e.piCursor[I];
        preIntegratedTable_().access(ic)=e.integratedState;
    }

}

/*
template<class IndexDriver>
int SIJPDFthermo<IndexDriver>::elementsdone=0;
    */

template<class IndexDriver>
void SIJPDFthermo<IndexDriver>::fillIntegrationSchedule //preIntegrate
(
  List<label>& pIcursor,
  DynamicList<integrationScheduleEntry>& schedule,
  label varDim
)
{
  
    if (varDim<preIntegratedTable_().dimensionality())
    {        
        for
            (
                pIcursor[varDim] = 0;
                pIcursor[varDim] < preIntegratedTable_().nElements(varDim);
                pIcursor[varDim]++
            )
        {
            /*preIntegrate*/fillIntegrationSchedule(pIcursor, schedule, varDim+1);
        }
    }
    else
    {
        /*
        elementsdone++;
        label elementstotal=1;
        for (label i=0;i<preIntegratedTable_().dimensionality();i++)
        {
            elementstotal*=preIntegratedTable_().nElements(i);
        }
        if ((100*elementsdone)%elementstotal==0)
            Info<<"PROGRESS: "<<(100*elementsdone)/elementstotal<<"%"<<endl;
            */

        // perform integration
        
        // variables over which will be integrated
        integrationParameter ix;    
        forAll(pdf_, I)
        {
            ix.insert(pdf_[I].parameterName(), 0.0);
        }
        
        
        // constant parameters during integration
        integrationParameter constants;
        for (ListDescription::const_iterator it=
                 rawTable_.progressVariableDescription().begin();
             it!=rawTable_.progressVariableDescription().end(); it++)
        {
            
            if ( !ix.found(it.key()) ) // constant
            {
                label pIdim=preIntegratedTable_().indexOfPV(it.key());
                constants.insert
                    (
                        it.key(), 
                        preIntegratedTable_().valueAt(pIdim, pIcursor[pIdim])
                    );
            }
        }
        
        chemicalSystemState integratedState
            (
                preIntegratedTable_().contentDescription()
            );
        
        /*
        performIntegration(pIcursor, integratedState, ix, constants);

        // save integrated state to preintegrated table
        preIntegratedTable_().access(pIcursor)=integratedState;
            */

        integrationScheduleEntry e;
        e.piCursor=pIcursor;
        e.ix=ix;
        e.constants=constants;
        e.integratedState=integratedState;

        schedule.append(e);
    }
}

template<class IndexDriver>
void SIJPDFthermo<IndexDriver>::performIntegration
(
    List<label>& pIcursor,
    chemicalSystemState& integratedState,
    integrationParameter& ix,
    integrationParameter& constants
)
{
    // fill integrated state
    sumOfIntegrableFunctions<scalar> pdfs=multipliedPdfs(pIcursor);

    const ListDescription& cdesc=preIntegratedTable_().contentDescription();

    label Irho=2*pdf_.size();

    // mean mixture properties
    scalar meanrho=
        integrator<scalar>::integrate
        (
            rawTable_.entry("rho") * pdfs,
            ix, constants
        );

    integratedState[Irho]=meanrho;
      
    // source terms for moment transport equations
    forAll(pdf_, I)
    {
          
        integratedState[2*I]=
            integrator<scalar>::integrate
            (
                rawTable_.entry("Ydot" & pdf_[I].parameterName())
                * pdfs,
                ix, constants
            );

        integratedState[2*I+1]=
            constNormValues_[I]*
            integrator<scalar>::integrate
            (
                singleParameter(pdf_[I].parameterName())
                * rawTable_.entry("Ydot" & pdf_[I].parameterName())
                * pdfs,
                ix, constants
            );
      
        word n="grad("&pdf_[I].parameterName()&")";
        if (cdesc.found(n))
        {
            integratedState[cdesc[n]]=
                (1./constNormValues_[I]) *
                integrator<scalar>::integrate
                (
                    rawTable_.entry(n) 
                    * pdfs,
                    ix, constants
                );
        }
    }

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

}






template<class IndexDriver>
void SIJPDFthermo<IndexDriver>::boundProgressVariable(const word& variable)
{
    label i=preIntegratedTable_().indexOfPV(variable);
    if ( variable=="h" )
    {
        h_().max(preIntegratedTable_().startValue(i));
        h_().min(preIntegratedTable_().maxValue(i));
        Info<<"min(h)="<<min(h_()).value()<<" max(h)="<<max(h_()).value()<<endl;
    } else
        boundProgressVariable(preIntegratedTable_().indexOfPV(variable));
}

template<class IndexDriver>
void SIJPDFthermo<IndexDriver>::boundProgressVariable(label i)
{
    label I=i/2;
    if (i%2==0)
    {
        // First moment
        firstMoment_[I].max(0.0);
        firstMoment_[I].min(1.0);
        Info<<"min(fm"<<I<<")="<<min(firstMoment_[I]).value()
            <<" max(fm"<<I<<")="<<max(firstMoment_[I]).value()<<endl;
    } else
    {
        // Second moment
        secondMoment_[I]().volScalarField::operator=
            (
                max
		(
                    min
                    (
                        secondMoment_[I](), 
                        firstMoment_[I]
	            ), 
                    sqr(firstMoment_[I])
                )
            );
        Info<<"min(sm"<<I<<")="<<min(secondMoment_[I]()).value()
            <<" max(sm"<<I<<")="<<max(secondMoment_[I]()).value()<<endl;
    }
}


template<class IndexDriver>
void SIJPDFthermo<IndexDriver>::boundProgressVariables()
{
    if (solvesEnthalpy()) boundProgressVariable("h");
    for (label i=0; i<pdf_.size()*2; i++)
        boundProgressVariable(i);
}


/**
    
    Creates directory of coordinates (progress variables) in preintegrated chemistry
    table

    New preintegrated chemistry table should have coordinates in the following order:
    - first moment of pdf 1
    - second moment of pdf 1
    - first moment of pdf 2
    - second moment of pdf 2
    (and so on for all pdfs)
    - remaining coordinates from raw table
*/
template<class IndexDriver>
HashTable<label,word> SIJPDFthermo<IndexDriver>::createPVDescription() const
{
    // create description of progress variables
    HashTable<label,word> new_pvdesc;
    
    // insert mean and variance of each progress variable as first table coordinates
    label pvindex=0;
    forAll(pdf_, I)
    {
        if (!rawTable_.progressVariableDescription().found(pdf_[I].parameterName()))
        {
            FatalErrorIn("SIJPDFthermo<IndexDriver>::SIJPDFthermo()")
                << " Selected progress variable " << pdf_[I].parameterName()
                    << " of PDF " << pdf_[I].name()
                    << " not found in provided chemistry table." << endl
                    << " Valid progress variables are:" << endl
                    << rawTable_.progressVariableDescription().toc() << endl
                    << exit(FatalError);
        }
        new_pvdesc.insert(pdf_[I].parameterName(), pvindex++);
        new_pvdesc.insert(pdf_[I].parameterName() & "_Var",  pvindex++);
    }
    
    // transfer all remaining table coordinates to new preintegrated table
    HashTable<label,word> cpvdesc(rawTable_.progressVariableDescription());
    forAll(pdf_, I)
    {
        cpvdesc.erase(pdf_[I].parameterName());
    }
    for (HashTable<label,word>::const_iterator it=cpvdesc.begin();
         it!=cpvdesc.end(); it++)
        new_pvdesc.insert(it.key(), pvindex++);
    
    return new_pvdesc;
}

/**

    Creates directory of fields in a cell in the preintegrated chemistry table

    Order of entries should be:
    - source term of first moment eqn for pdf 1
    - source term of second moment eqn for pdf 1
    - source term of first moment eqn for pdf 2
    - source term of second moment eqn for pdf 2
    (and so on for all pdfs)
    - density of mixture
    - temperature
    - lambda
    - Cp
    - mu
    - remaining fields
*/
template<class IndexDriver>
HashTable<label,word> SIJPDFthermo<IndexDriver>::createContentDescription() const
{
    // create description of preintegrated table content
    HashTable<label,word> new_contentdesc;
    
    // insert source terms of mean and variance equation of each
    // progress variable as first table entries
    label cindex=0;
    forAll(pdf_, I)
    {
        new_contentdesc.insert("Ydot" & pdf_[I].parameterName(), cindex++);
        new_contentdesc.insert
            ("Ydot" & pdf_[I].parameterName() & "_" & pdf_[I].parameterName(),  cindex++);
    }
    
    HashTable<label,word> ccontentdesc(rawTable_.contentDescription());
    
    new_contentdesc.insert("rho", cindex++);  ccontentdesc.erase("rho");
    new_contentdesc.insert("T", cindex++);  ccontentdesc.erase("T");
    new_contentdesc.insert("lambda", cindex++);  ccontentdesc.erase("lambda");
    new_contentdesc.insert("Cp", cindex++);  ccontentdesc.erase("Cp");
    new_contentdesc.insert("mu", cindex++);  ccontentdesc.erase("mu");
    
    // transfer all remaining table entries to new preintegrated table
    for (HashTable<label,word>::const_iterator it=ccontentdesc.begin();
         it!=ccontentdesc.end(); it++)
        new_contentdesc.insert(it.key(), cindex++);
    
    return new_contentdesc;
}


/**
    Creates an (empty) table to store preintegrated chemistry source terms and
    thermophysical properties with first and second moments of PDF's as coordinates
*/
template<class IndexDriver>
void SIJPDFthermo<IndexDriver>::createPreIntegratedTable()
{

    HashTable<label,word> new_contentdesc
        =createContentDescription();
    HashTable<label,word> new_pvdesc
        =createPVDescription();

    // create preintegrated lookup table
    label num[new_pvdesc.size()];
    scalar start[new_pvdesc.size()*2];
    scalar delta[new_pvdesc.size()*2];

    // set up first and second moments
    forAll(pdf_, I)
    {
        num[2*I]=pdf_[I].tabResolutionFM();
        num[2*I+1]=pdf_[I].tabResolutionSM();
    }

    for (label i=0;i<pdf_.size()*2;i++)
    {
        start[i]=0.0;
        delta[i]=1.0/scalar(num[i]-1);
    }
    

    // set up remaining coordinates (h, ...)
    for (label i=pdf_.size()*2;i<new_pvdesc.size();i++)
    {

        // find name of current progress variable i
        word pv="";
        for (HashTable<label,word>::const_iterator it=new_pvdesc.begin();
             it!=new_pvdesc.end(); it++)
            if (it()==i) pv=it.key();            
        
        // use same resolution as in raw table
        num[i]=rawTable_.nElements(rawTable_.progressVariableDescription()[pv]);
        start[i]=rawTable_.startValue(rawTable_.progressVariableDescription()[pv]);
        delta[i]=rawTable_.deltaValue(rawTable_.progressVariableDescription()[pv]);
    }

    preIntegratedTable_.reset
        (
            new chemistryTable
            (
                new_pvdesc.size(),
                num, start, delta,
                new_contentdesc,
                new_pvdesc
            )
        );    
}



template<class IndexDriver>
void SIJPDFthermo<IndexDriver>::initialize()
{
     
    bool cacheMatches=readPreIntegratedTableFromCache();

    if (cacheMatches)
    {
        /*
        for (HashTable<label,word>::const_iterator it=new_pvdesc.begin();
             it!=new_pvdesc.end(); it++)
            cacheMatches=
                cacheMatches
                &&
                preIntegratedTable_().progressVariableDescription().found(it.key());
        
        for (HashTable<label,word>::const_iterator it=new_contentdesc.begin();
             it!=new_contentdesc.end(); it++)
            cacheMatches=
                cacheMatches
                &&
                preIntegratedTable_().contentDescription().found(it.key());
        
        for (label i=0;i<preIntegratedTable_().dimensionality();i++)
            cacheMatches=
                cacheMatches
                &&
                (preIntegratedTable_().nElements(i)==num[i]);
            */
    }
    else
    {

        Info << endl << "Creating preintegrated chemistry table" << endl;
        
        createPreIntegratedTable();          
        
        // Integrate with pdfs over raw chemistry table
        
        // index into preintegrated table
        label icursor[preIntegratedTable_().dimensionality()];
        memset(&icursor, sizeof(icursor), 0);
          
        Info << endl << "Integrating chemistry table" << endl;

        //elementsdone=0;
        preIntegrate();
        
        writeCache();
      }
    
    
    if (preIntegratedTable_().progressVariableDescription().found("h"))
    {

        Info << "Creating field h"<<endl<<endl;

        h_.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        "h",
                        psi_.mesh().time().timeName(),
                        psi_.mesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    psi_.mesh(),
                    dimensionSet(0, 2, -2, 0, 0),
                    hBoundaryTypes()
                )
            );
        
        mvfields_.add(h_());
        
        forAll(h_(), celli)
        {
            h_()[celli] = lookupCellByTemp(T_[celli], celli);
        }
        
        forAll(h_().boundaryField(), patchi)
        {
               forAll (h_().boundaryField()[patchi], facei)
               {
                   h_().boundaryField()[patchi][facei] = 
                       lookupFaceByTemp
                       (
                           T_.boundaryField()[patchi][facei], 
                           patchi, 
                           facei
                       );
               }
        }
           
        hBoundaryCorrection(h_());
        
        h_().write();
           
    }

    
    // -----------------------------------------
 
    boundProgressVariables();
    correct();
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class IndexDriver>
SIJPDFthermo<IndexDriver>::SIJPDFthermo
(
    const fvMesh& mesh
)
    : abstractProgressVariableThermo(mesh),
      rawTableDict_
      (
       IOobject
       (
	word(lookup("LDMtableName")),
	mesh.time().constant(),
	mesh,
	IOobject::MUST_READ,
	IOobject::NO_WRITE
	)
      ),      
      rawTable_(rawTableDict_.lookup("tableData")),                                    
      pdf_(lookup("compositionPDFs"), pdf::iNew()),
      constNormValues_(pdf_.size(), 1.0)
{

    if (found("wallQuenchingDampingFunction"))
    {
        quench_=wallQuenchingDampingFunction::New(mesh, subDict("wallQuenchingDampingFunction"));
    }
    
    // Create directory of variable names
    forAll(pdf_, I)
        variable_.insert(pdf_[I].parameterName(), I);

    forAllVariables(var)
    {
        // check upper and lower bounds
        scalar lower=rawTable_.startValue(rawTable_.indexOfPV(var.key()));
        scalar upper=rawTable_.maxValue(rawTable_.indexOfPV(var.key()));
        if 
            (
                (lower!=0.0)
                ||
                mag(upper-1.0)>1e-3
            )
        {
            FatalErrorIn("SIJPDFthermo<IndexDriver>::SIJPDFthermo()")
                << " Tabulation range of variable " << var.key()
                    << " is ["<<lower<<","<<upper<<"]." << endl
                    << " Variables, which will be integrated over PDFs "
                    << "must be in range [0,1]." << endl
                    << exit(FatalError);
        }

        // look for normalization values
        word keyword="constantNormValue_"+var.key();
        if (rawTableDict_.found(keyword))
            constNormValues_[var()]=readScalar(rawTableDict_.lookup(keyword));
    }

    // * * * * * * * * * * *  Initialize fields  * * * * * * * * * * *

    // First moments
    firstMoment_.setSize(pdf_.size());
    forAll(firstMoment_, I)
    {
        firstMoment_.set
            (
                I,
                new volScalarField
                (
                    IOobject
                    (
                        pdf_[I].name() & pdf_[I].parameterName() & "_Mean",
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
        mvfields_.add(firstMoment_[I]);
    }

    // Second moments
    secondMoment_.setSize(pdf_.size());
    forAll(secondMoment_, I)
    {
        secondMoment_.set
            (
                I,
                secondMomentSolver::New
                (
                    mesh,
                    *this,
                    pdf_[I].name() & pdf_[I].parameterName(),
                    subDict(pdf_[I].name() & pdf_[I].parameterName()&"secondMoment"),
                    mvfields_
                )
            );
        //mvfields_.add(secondMoment_[I]());
    }

    // source terms of first moments
    meanWdot_.setSize(pdf_.size());
    forAll(meanWdot_, I)
      meanWdot_.set
          (
           I,
           new volScalarField
           (
              IOobject
              (
                pdf_[I].name() & "meanW" & pdf_[I].parameterName() & "dot",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
              ),
              mesh,
              dimensionedScalar("meanWdot", dimless/dimTime, 0.0)
            )
          );
   
    // source terms of second moments
    meanWdotC_.setSize(pdf_.size());
    forAll(meanWdotC_, I)
      meanWdotC_.set
          (
           I,
           new volScalarField
           (
            IOobject
            (
              pdf_[I].name() & "meanW" 
              & pdf_[I].parameterName() 
              & "dot_" & pdf_[I].parameterName(),
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("meanWdotC", dimless/dimTime, 0.0)
           )
          );

    // correct() is called through initialize()
    // since preintegration is not done yet

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class IndexDriver>
SIJPDFthermo<IndexDriver>::~SIJPDFthermo()
{}



template<class IndexDriver>
tmp<scalarField> SIJPDFthermo<IndexDriver>::h
(
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> th(new scalarField(T.size()));
    scalarField& h = th();

    forAll(T, celli)
    {
        h[celli] = lookupCellByTemp(T[celli], cells[celli]);
    }

    return th;
}


template<class IndexDriver>
tmp<scalarField> SIJPDFthermo<IndexDriver>::h
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> th(new scalarField(T.size()));
    scalarField& h = th();

    forAll(T, facei)
    {
        h[facei] = lookupFaceByTemp(T[facei], patchi, facei);
    }

    return th;
}
 

template<class IndexDriver>
tmp<scalarField> SIJPDFthermo<IndexDriver>::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCp(new scalarField(T.size()));
    scalarField& cp = tCp();

    forAll(T, facei)
    {
        chemicalSystemState state;
        lookupFaceByTemp(T[facei], patchi, facei, &state);

        //REVERSELOOKUPH(.boundaryField()[patchi][facei], T[facei], state, reverseLookup);
        cp[facei] = state[preIntegratedTable_().indexOf("Cp")];;
    }

    return tCp;
}

template<class IndexDriver>
tmp<volScalarField> SIJPDFthermo<IndexDriver>::Cp() const
{
    const fvMesh& mesh = T_.mesh();

    tmp<volScalarField> tCp
    (
        new volScalarField
        (
            IOobject
            (
                "Cp",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 2, -2, -1, 0)
        )
    );

    volScalarField& cp = tCp();

    forAll(T_, celli)
    {
        chemicalSystemState state;
        lookupCellByTemp(T_[celli], celli, &state);
        cp[celli] = state[preIntegratedTable_().indexOf("Cp")];;
    }

    forAll(T_.boundaryField(), patchi)
    {
        cp.boundaryField()[patchi] = Cp(T_.boundaryField()[patchi], patchi);
    }

    return tCp;
}
/*
tmp<volScalarField> hThermo<MixtureType>::Cv() const
{
    const fvMesh& mesh = T_.mesh();

    tmp<volScalarField> tCv
    (
        new volScalarField
        (
            IOobject
            (
                "Cv",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 2, -2, -1, 0),
            T_.boundaryField().types()
        )
    );

    volScalarField& cv = tCv();

    forAll(T_, celli)
    {
        cv[celli] = this->cellMixture(celli).Cv(T_[celli]);
    }


    forAll(T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pT = T_.boundaryField()[patchi];
        fvPatchScalarField& pCv = cv.boundaryField()[patchi];

        forAll(pT, facei)
        {
            pCv[facei] = this->patchFaceMixture(patchi, facei).Cv(pT[facei]);
        }
    }

    return tCv;
}
    */  
/*

#define LOOKUP(_i_)                                                     \
    INDEX(_i_)                                                          \
        if (h_.valid())                                                 \
        {                                                               \
            label i=preIntegratedTable_().indexOfPV("h");               \
            st[i]=h_()_i_;                                              \
        }                                                               \
        return preIntegratedTable_().lookup(st);                 


template<class Finder>
chemicalSystemState SIJPDFthermo<IndexDriver>::lookupCell(label cellI) const
{
  LOOKUP( [cellI] )
}


template<class Finder>
chemicalSystemState SIJPDFthermo<IndexDriver>::lookupFace(label patchI, label faceI) const
{
  LOOKUP( .boundaryField()[patchI][faceI] )
}
    */



template<class IndexDriver>
void SIJPDFthermo<IndexDriver>::correct()
{

  psi_.oldTime(); // Important! Pressure correction won't work if omitted

#ifdef DEBUGFILE
  dbgfile()<<endl<<endl<<endl;
#endif

  label Irho=2*pdf_.size();
  
  forAll(psi_, cellI)
    {
        chemicalSystemState state=lookupCell(cellI);
        
#ifdef COMBUSTAPPROX
#warning Combustion approximation of density is activated (intended for testing purposes only)
        scalar theta=firstMoment(0)[cellI]/0.161685;
        scalar rho=1.18114*(1.0/(1.0+0.86*theta/(1-0.86)));
        forAll(pdf_, I)
        {
          //meanWdot_[I][cellI]=firstMoment(0)[cellI]*rho*(1-theta)*exp(-18.4*(1-theta)/(1-0.86*(1-theta)));
          meanWdot_[I][cellI]=(1200/0.15)*pow(theta, 2)*exp(-30*sqr(theta - 0.3));//state[ 2*I   ];
          meanWdotC_[I][cellI]=0.0;
        }
        psi_[cellI]=rho/1e5;
#else
        forAll(pdf_, I)
        {
          meanWdot_[I][cellI]=state[ 2*I   ];
          meanWdotC_[I][cellI]=state[ 2*I+1 ];
        }
#warning isobaric table for p=1bar assumed!        
        psi_[cellI]=state[Irho]/1e5;
#endif
        //Info<<state[Irho]<<" "<<psi_[cellI]<<endl;
        T_[cellI]=state[Irho+1];//*pow(p_[cellI]/1e5, 0.285);
#warning gas constant of air assumed
        alpha_[cellI]=state[Irho+2]/state[Irho+3];
        mu_[cellI]=state[Irho+4];
    }

  if (solvesEnthalpy())
  {
      forAll(psi_.boundaryField(), pI)
      {
          fvPatchScalarField& pT = T_.boundaryField()[pI]; 
          fvPatchScalarField& ph = h_().boundaryField()[pI]; 
          
          if (pT.fixesValue())
          {
              forAll(pT, facei)
              {
                  ph[facei] = lookupFaceByTemp(pT[facei], pI, facei);
              }
          }
      }
  }

  forAll(psi_.boundaryField(), pI)
  {
      forAll(psi_.boundaryField()[pI], faceI)
      {
          chemicalSystemState state=lookupFace(pI, faceI);
          
          forAll(pdf_, I)
          {
              meanWdot_[I].boundaryField()[pI][faceI]=state[ 2*I   ];
              meanWdotC_[I].boundaryField()[pI][faceI]=state[ 2*I+1 ];
          }
          
          if (!T_.boundaryField()[pI].fixesValue())
              T_.boundaryField()[pI][faceI]=state[Irho+1];//*pow(p_.boundaryField()[pI][faceI]/1e5, 0.285);;
#ifdef COMBUSTAPPROX
          psi_.boundaryField()[pI][faceI]
           =1.18114*(1.0/(1.0+0.86*(firstMoment(0).boundaryField()[pI][faceI]/0.161685)/(1-0.86)))/1e5;
#else
          psi_.boundaryField()[pI][faceI]=state[Irho]/1e5;
#endif
          alpha_.boundaryField()[pI][faceI]=state[Irho+2]/state[Irho+3];
          mu_.boundaryField()[pI][faceI]=state[Irho+4];
      }
  }
#ifdef DEBUGFILE 
 dbgfile()<<endl<<endl<<endl;
#endif
  //Info<<psi()<<endl;
}

template<class IndexDriver>
void SIJPDFthermo<IndexDriver>::registerScalarFields
(
    compressible::LESModel& model
)
{
  /*
 if (h_.valid())
 	model.registerScalarField(h_());

 forAll(firstMoment_, I)
        model.registerScalarField(firstMoment_[I]);

 forAll(secondMoment_, I)
        secondMoment_[I].registerFields(model);
  */
}


template<class IndexDriver>
void SIJPDFthermo<IndexDriver>::solveProgressVariables
(
    const volScalarField& rho,
    const surfaceScalarField& phi,
    const compressible::LESModel& model
)
{
    tmp<fv::convectionScheme<scalar> > mvConvection
        (
            fv::convectionScheme<scalar>::New
            (
                phi.mesh(),
                mvfields_,
                phi,
                phi.mesh().divScheme("div(phi,mvScalars)")
            )
        );
    
    if (solvesEnthalpy())
    {
        solve
            (
                fvm::ddt(rho, h_())
                + mvConvection->fvmDiv(phi, h_())
                == 
                fvm::laplacian(model.alphaEff(), h_())
                //model.divRhoSgsFlux(h_())
            );
        boundProgressVariable("h");
    }

  forAll(firstMoment_, I)
  {
      

      fvScalarMatrix fmEqn 
          (
              fvm::ddt(rho, firstMoment_[I])
              + mvConvection->fvmDiv(phi, firstMoment_[I])
              - rho*meanWdot_[I]
              ==
              //model.divRhoSgsFlux(firstMoment_[I])
	      fvm::laplacian(model.alphaEff(), 
			     firstMoment_[I], 
			     "laplacian(Deff,"+firstMoment_[I].name()+")")
              //fvm::laplacian(alpha()/1.45, firstMoment_[I], "laplacian(Deff,CO2)")
          );
      
      if (quench_.valid())
      {
          solve(fmEqn == quench_()(firstMoment_[I]));
      }
      else
      {
          solve(fmEqn);
      }
      
      boundProgressVariable(2*I);
  }
  
  forAll(secondMoment_, I)
  {
      secondMoment_[I].solve
          (
              rho,
              phi,
              firstMoment_[I],
              meanWdotC_[I],
              mvConvection(),
              model,
              preIntegratedTable_()
          );
      boundProgressVariable(2*I+1);
  }
}

template<class IndexDriver>
void SIJPDFthermo<IndexDriver>::solveProgressVariables
(
    const volScalarField& rho,
    const surfaceScalarField& phi,
    const compressible::RASModel& model
)
{
    tmp<fv::convectionScheme<scalar> > mvConvection
        (
            fv::convectionScheme<scalar>::New
            (
                phi.mesh(),
                mvfields_,
                phi,
                phi.mesh().divScheme("div(phi,mvScalars)")
            )
        );
    
    if (solvesEnthalpy())
    {
        solve
            (
                fvm::ddt(rho, h_())
                + mvConvection->fvmDiv(phi, h_())
                == 
                fvm::laplacian(model.alphaEff(), h_())
            );
        boundProgressVariable("h");
    }

  forAll(firstMoment_, I)
  {
      

      fvScalarMatrix fmEqn 
          (
              fvm::ddt(rho, firstMoment_[I])
              + mvConvection->fvmDiv(phi, firstMoment_[I])
              - rho*meanWdot_[I]
              ==
              fvm::laplacian(model.muEff(), firstMoment_[I])
          );
      
      if (quench_.valid())
      {
          solve(fmEqn == quench_()(firstMoment_[I]));
      }
      else
      {
          solve(fmEqn);
      }
      
      boundProgressVariable(2*I);
  }
  
  forAll(secondMoment_, I)
  {

      Foam::solve
          (
              fvm::ddt(rho, secondMoment_[I]() )
              + mvConvection->fvmDiv(phi, secondMoment_[I]() )
              - 2.0*rho*meanWdotC_[I]
              ==
              fvm::laplacian(model.muEff(), secondMoment_[I]() )
              -2.0*
              (
                  model.mu()*magSqr(fvc::grad(firstMoment_[I]))
                  +
                  2.0*rho*(model.epsilon()/model.k())
                  *( secondMoment_[I]() - sqr(firstMoment_[I]) )
              )
          );
      
      boundProgressVariable(2*I+1);
  }
}

template<class IndexDriver>
fileName SIJPDFthermo<IndexDriver>::cacheFile() const
{
  fileName f;
  forAll(pdf_, I)
  {
    if (I>0) f = f & "-";
    f = f & pdf_[I].name() & pdf_[I].parameterName();
  }
  return f;
}


template<class IndexDriver>
bool SIJPDFthermo<IndexDriver>::readPreIntegratedTableFromCache()
{
    const Time& time=psi_.mesh().time();

    fileName fullPath;
    if (Pstream::parRun())
    {
        fullPath = time.path()/".."/time.constant()/cacheFile();
    }
    else
    {
        fullPath = time.path()/time.constant()/cacheFile();;
    }

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

        Info << "Reading preintegrated table from " 
            << fullPath << endl;

        preIntegratedTable_.reset
            (
                new chemistryTable(tableDict.lookup("tableData"))
            );

        return true;
    }
    else return false;
}



template<class IndexDriver>
void SIJPDFthermo<IndexDriver>::writeCache() const
{
    if (Pstream::master())
    {
        const Time& time=psi_.mesh().time();
        
        fileName fullPath;
        if (Pstream::parRun())
        {
            fullPath = time.path()/".."/time.constant()/cacheFile();
        }
        else
        {
            fullPath = time.path()/time.constant()/cacheFile();;
        }

        IOdictionary tableDict
            (
                IOobject
                (
                    cacheFile(),
                    psi_.mesh().time().constant(),
                    psi_.mesh().time(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                )
            );
        
        OFstream tableDictFile(fullPath);
        
        tableDict.writeHeader(tableDictFile);
        
        tableDictFile
            <<"tableData"
                <<token::SPACE
                <<preIntegratedTable_()
                <<token::END_STATEMENT
                <<endl;
        /*    
            if (mixtureFractionIsVariable_)
            {
            tableDictFile
            <<"normFactors"
            <<token::SPACE
            <<normFactors_
            <<token::END_STATEMENT
            <<endl;
            } else
            */
        forAllVariables(var)
        {
            tableDictFile
                <<"constantNormValue_"<<var.key()
                    <<token::SPACE
                    <<constNormValues_[var()]
                    <<token::END_STATEMENT
                    <<endl;
        }
    }
}

template<class IndexDriver>
bool SIJPDFthermo<IndexDriver>::read()
{
  return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
