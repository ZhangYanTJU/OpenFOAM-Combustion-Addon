#include "SIJPDFvarZthermo.H"
#include "fvMesh.H"
#include "OFstream.H"

#include "beta.H"
#include "addToRunTimeSelectionTable.H"
#include "integrableFunction.H"
#include "integrator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    //defineTypeNameAndDebug(SIJPDFvarZthermoStateFinder, 0);

#define EPSILON 1e-3    

#define INDEX(_i_) \
                                                                        \
    const SIJPDFvarZthermo/*<SIJPDFvarZthermoStateFinder>*/& thermo_=   \
        *static_cast<const SIJPDFvarZthermo/*<SIJPDFvarZthermoStateFinder>*/ *>(t); \
                                                                        \
        List<scalar> st                                                 \
            (                                                           \
                thermo_.preIntegratedTable().progressVariableDescription().size(), \
                0.0                                                     \
            );                                                          \
                                                                        \
        scalarField norm(3, 1.0);                                       \
        List<scalar> mZ(2, 0.0);                                        \
                                                                        \
        mZ[0]=thermo_.firstMoment(thermo_.mixFracI_)_i_;                \
        mZ[1]=thermo_.secondMoment(thermo_.mixFracI_)_i_;               \
        scalar denomZ=(mZ[0]*(1.0-mZ[0]));                              \
        mZ[1]=denomZ>EPSILON?                                           \
               (mZ[1]-sqr(mZ[0]))/denomZ : 0.0;	                        \
        mZ[0]/=thermo_.constNormValues_[thermo_.mixFracI_];             \
        /*Info<<mZ[0]<<" "<<mZ[1]<<endl;*/                              \
                                                                        \
        st[2*thermo_.mixFracI_]=mZ[0];                                  \
        st[2*thermo_.mixFracI_+1]=mZ[1];                                \
                                                                        \
        forAllVariablesIn(thermo_, var)                                 \
        {                                                               \
            norm = thermo_.derivedNormFactors_[var.key()]->lookup(mZ);  \
                                                                        \
            scalar M = min                                              \
                (                                                       \
                    1.0,                                                \
                    thermo_.firstMoment(var())_i_ * norm[0]             \
                );                                                      \
            st[2*var()] = M;                                            \
                                                                        \
            if (M>EPSILON)                                              \
            {                                                           \
                scalar m1=thermo_.firstMoment(var())_i_;                \
                scalar m2=thermo_.secondMoment(var())_i_;               \
                                                                        \
                /*scalar denom= ( m1 * ( 1.0 - m1) );                   \
                    scalar S =                                          \
                    denom > EPSILON                                     \
                    ?                                                   \
                    (m2-sqr(m1)) / denom                                \
                    :                                                   \
                    0.0;*/                                              \
                scalar v=max( m2-sqr(m1), 0.0);                         \
                /*     if (norm[2]*v<EPSILON)                           \
                    Info<<norm[2]<<"*"<<v<<endl;     */                 \
                scalar S =                                              \
                    (                                                   \
                        norm[2] + v - 2.0 *                             \
                        sqrt(max(norm[2] * v, 0.0))                     \
                    ) * norm[1];                                        \
                scalar denom= ( M * ( 1.0 - M) );                       \
                if (denom>EPSILON) S/=denom; else S=0.0;                \
                /*if ((S>1)||(M>1))                                     \
                    Info<<"S="<<S<<" M="<<M<<endl;*/                    \
                st[2*var()+1] = min(1.0, S);                            \
            } else st[2*var()+1]=0.0;                                   \
        }                                                               \
                                                                        \
        if (thermo_.solvesEnthalpy())                                   \
        {                                                               \
            label i=thermo_.preIntegratedTable().indexOfPV("h");        \
            st[i]=thermo_.h()_i_;                                       \
        }                                                               \
    

scalar SIJPDFvarZthermoStateFinder::lookupCellByTemp
(
    const abstractProgressVariableThermo* t,
    scalar temperature, 
    label cellI,
    chemicalSystemState* state
)
{
    scalar h;

    INDEX([cellI]);
        
    chemicalSystemState s=
        thermo_.preIntegratedTable().reverseLookup
        (
            thermo_.preIntegratedTable().indexOfPV("h"),
            st,
            thermo_.preIntegratedTable().indexOf("T"),
            temperature,
            h
        );

    if (state!=NULL) *state=s;
    return h;
}

scalar SIJPDFvarZthermoStateFinder::lookupFaceByTemp
(
    const abstractProgressVariableThermo* t,
    scalar temperature,
    label patchI,
    label faceI,
    chemicalSystemState* state
)
{
    scalar h;

    INDEX(.boundaryField()[patchI][faceI]);
        
    chemicalSystemState s=
        thermo_.preIntegratedTable().reverseLookup
        (
            thermo_.preIntegratedTable().indexOfPV("h"),
            st,
            thermo_.preIntegratedTable().indexOf("T"),
            temperature,
            h
        );

    if (state!=NULL) *state=s;
    return h;
}

chemicalSystemState SIJPDFvarZthermoStateFinder::lookupCell
(
    const abstractProgressVariableThermo* t,
    label cellI
)
{
    INDEX( [cellI] );
    return thermo_.preIntegratedTable().lookup(st);
}

chemicalSystemState SIJPDFvarZthermoStateFinder::lookupFace
(
    const abstractProgressVariableThermo* t,
    label patchI,
    label faceI
)
{
    INDEX( .boundaryField()[patchI][faceI] )
    return thermo_.preIntegratedTable().lookup(st);
}

#undef INDEX

}
