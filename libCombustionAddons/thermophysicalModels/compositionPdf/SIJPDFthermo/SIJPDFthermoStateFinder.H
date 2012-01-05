
#ifndef SIJPDFTHERMOSTATEFINDER_H
#define SIJPDFTHERMOSTATEFINDER_H

#include "SIJPDFthermo.H"

namespace Foam
{

#define EPSILON 1e-3


#define SIJPDFTHERMO_NORM_M1(I, _i_)                            \
                                                                        \
    const SIJPDFthermo<SIJPDFthermoIndexDriver>& thermo=                \
        static_cast<const SIJPDFthermo<SIJPDFthermoIndexDriver>&>       \
    (*th);                                                              \
        scalar normM1 = thermo.firstMoment(I)_i_ /                      \
            thermo.constNormValues()[I];                                \

#ifdef DEBUGFILE
#define OUTPUTCMD dbgfile() << " "<<m1 << " "<<S<<" "<<m2;
#else
#define OUTPUTCMD ;
#endif

#define SIJPDFTHERMO_INDEX(_i_)                                         \
                                                                        \
    const SIJPDFthermo<SIJPDFthermoIndexDriver>& thermo=                \
        static_cast<const SIJPDFthermo<SIJPDFthermoIndexDriver>&>       \
    (*th);                                                              \
                                                                        \
    List<scalar> st                                                     \
    (                                                                   \
        thermo.preIntegratedTable().progressVariableDescription().size(), \
        0.0                                                             \
    );                                                                  \
                                                                        \
    forAll(thermo.presumedPdfs(), I)                                    \
    {                                                                   \
        scalar m1=thermo.firstMoment(I)_i_;                             \
        scalar m2=thermo.secondMoment(I)()_i_;                          \
        scalar m1n=thermo.constNormValues()[I];                         \
        scalar denom=m1n*m1-m1*m1;   \
        scalar nom=m2-sqr(m1); \
        st[2*I  ] = m1 / m1n;                   \
                                                                        \
        scalar S = denom>EPSILON? nom/denom : (nom<EPSILON? 0.0:1.0);          \
        /*if (S-1.0>1e-4) Info<<"SF S="<<S<<" m1="<<m1<<endl;*/                \
        st[2*I+1] = S;                                                  \
        OUTPUTCMD; \
    }                                                                   \
                                                                        \
    if (thermo.solvesEnthalpy())                                        \
    {                                                                   \
        label i=thermo.preIntegratedTable().indexOfPV("h");             \
        st[i]=thermo.h()_i_;                                            \
    }


    template<class IndexDriver>
    class SIJPDFthermo;

    class SIJPDFthermoIndexDriver
    {
    public:

        static inline scalar cellNormM1
        (
            const abstractProgressVariableThermo* th,
            label cellI,
            label pdfI
        ) 
            {
                SIJPDFTHERMO_NORM_M1( pdfI, [cellI] );
                return normM1;
            }

        static inline scalar faceNormM1
        (
            const abstractProgressVariableThermo* th,
            label patchI, 
            label faceI,
            label pdfI
        ) 
            {
                SIJPDFTHERMO_NORM_M1( pdfI, .boundaryField()[patchI][faceI] );
                return normM1;
            }

        static inline List<scalar> cellTableIndex
        (
            const abstractProgressVariableThermo* th,
            label cellI
        ) 
            {
                SIJPDFTHERMO_INDEX( [cellI] );
#ifdef DEBUGFILE
		dbgfile()<<nl;
#endif
                return st;
            }

        static inline List<scalar> faceTableIndex
        (
            const abstractProgressVariableThermo* th,
            label patchI, 
            label faceI
        ) 
            {
                SIJPDFTHERMO_INDEX( .boundaryField()[patchI][faceI] )
#ifdef DEBUGFILE
		dbgfile()<<nl;
#endif
                return st;
            }
    };
}

#endif
