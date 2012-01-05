#include "FSDSecondMoment.H"
#include "addToRunTimeSelectionTable.H"
#include "SIJPDFthermo.H"
#include "SIJPDFthermoStateFinder.H"

namespace Foam
{

    //defineTypeNameAndDebug(FSDSecondMoment, 0);
    //addToRunTimeSelectionTable(secondMomentSolver, FSDSecondMoment, dictionary);

template<class Thermo>
FSDSecondMoment<Thermo>::FSDSecondMoment
(
    const fvMesh& mesh,
    const abstractProgressVariableThermo& thermo,
    const word& variableName,
    const dictionary& dict,
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable& mvtab
)
    : secondMomentSolver(mesh, thermo, variableName, dict, mvtab),
      //C_(readScalar(dict.lookup("C"))),
      Sigma_
    (
        IOobject
        (
            "Sigma",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless/dimLength, 0.0),
        zeroGradientFvPatchScalarField::typeName
    )
{
}


#define BOUND 1e-5


template<class Thermo>
void FSDSecondMoment<Thermo>::solve
(
    const volScalarField& rho,
    const surfaceScalarField& phi,
    const volScalarField& firstMoment,
    const volScalarField& meanWdotC,
    fv::convectionScheme<scalar>& mvConvection,
    const compressible::LESModel& model,
    const chemistryTable& preIntegratedTable
)
{

    updateSigma
        (
            rho, 
            phi, 
            firstMoment, 
            meanWdotC, 
            mvConvection, 
            model, 
            preIntegratedTable
        );

    Sigma_.correctBoundaryConditions();

    Info<<"max Sigma="<<Foam::max(Sigma_).value()<<endl;
    Info<<"total flame surface="<<fvc::domainIntegrate(Sigma_).value()<<endl;

    label pIdx=preIntegratedTable.indexOfPV(varName()&"_Var");
    label cIdx=preIntegratedTable.indexOf("grad("&varName()&")");
    string error;
    forAll( (*this)(), cellI)
    {
        scalar S=0.0;
        scalar Sigma=Sigma_[cellI];
        if (Sigma>BOUND)
        {
	//List<scalar> cursor=Thermo::ThermoIndexDriver::cellTableIndex(&thermo(), cellI);
	//Info<<cursor<<endl;
	//error="n";
            preIntegratedTable.reverseLookup
                (
                    pIdx,
                    Thermo::ThermoIndexDriver::cellTableIndex(&thermo(), cellI),
                    cIdx,
                    Sigma,
                    S,
                    error
                );
        } else S=1.0;
        S=max(1e-4, min(0.9999, S));
	//dbgfile() << " "<<firstMoment[cellI]<<" " <<S<<" "<<Sigma<<nl;
	//if (S-1.0>1e-4) Info<<"FSD #"<<cellI<<": S="<<S<<endl;
        (*this)()[cellI]=S*firstMoment[cellI]+(1.-S)*sqr(firstMoment[cellI]);
    }
    //dbgfile()<<endl<<endl;
    forAll( (*this)().boundaryField(), patchI)
        forAll( (*this)().boundaryField()[patchI], faceI)
    {
        scalar S=0.0;
        scalar Sigma=Sigma_.boundaryField()[patchI][faceI];
        if (Sigma>BOUND)
        {
            preIntegratedTable.reverseLookup
                (
                    pIdx,
                    Thermo::ThermoIndexDriver::faceTableIndex(&thermo(), patchI, faceI),
                    cIdx,
                    Sigma,
                    S,
                    error
                );
        } else S=1.0;
        S=max(1e-4, min(0.9999, S));
        (*this)().boundaryField()[patchI][faceI]=
            S*firstMoment.boundaryField()[patchI][faceI]
            +(1.-S)*sqr(firstMoment.boundaryField()[patchI][faceI]);
    }

}

}
