#ifndef SIJPDFVARZTHERMOSTATEFINDER_H
#define SIJPDFVARZTHERMOSTATEFINDER_H

namespace Foam
{
#define EPSILON 1e-4

#define SIJPDFVARZTHERMO_NORM_M1(I, _i_)                            \
                                                                        \
    const SIJPDFvarZthermo& thermo=                                     \
        static_cast<const SIJPDFvarZthermo&>                            \
        (*th);                                                          \
                                                                        \
    List<scalar> mZ(2, 0.0);                                            \
                                                                        \
    scalar nZ=thermo.constNormValues_[thermo.mixFracI_];                \
    mZ[0]=thermo.firstMoment(thermo.mixFracI_)_i_;                      \
    mZ[1]=thermo.secondMoment(thermo.mixFracI_)()_i_;                   \
    scalar denomZ=(nZ*mZ[0]*-mZ[0]*mZ[0]);                              \
    mZ[1]=denomZ>EPSILON?                                               \
    (mZ[1]-sqr(mZ[0]))/denomZ : 0.0;                                    \
    mZ[0]/=nZ;                   \
                                                                        \
    scalarField norm = thermo.derivedNormFactors_                       \
        [thermo.presumedPdfs()[I].parameterName()]->lookup(mZ);         \
                                                                        \
    scalar normM1 = min                                                 \
        (                                                               \
            1.0,                                                        \
            thermo.firstMoment(I)_i_ * norm[0]                          \
        );                                                              \
    

#define SIJPDFVARZTHERMO_INDEX(_i_) \
                                                                        \
    const SIJPDFvarZthermo& thermo=                                     \
        *static_cast<const SIJPDFvarZthermo*>(th);                      \
                                                                        \
        List<scalar> st                                                 \
            (                                                           \
                thermo.preIntegratedTable().progressVariableDescription().size(), \
                0.0                                                     \
            );                                                          \
                                                                        \
        scalarField norm(3, 1.0);                                       \
        List<scalar> mZ(2, 0.0);                                        \
                                                                        \
        scalar nZ=thermo.constNormValues_[thermo.mixFracI_];            \
        mZ[0]=thermo.firstMoment(thermo.mixFracI_)_i_;                  \
        mZ[1]=thermo.secondMoment(thermo.mixFracI_)()_i_;               \
        scalar denomZ=(nZ*mZ[0]-mZ[0]*mZ[0]);                           \
        mZ[1]=denomZ>EPSILON?                                           \
               (mZ[1]-sqr(mZ[0]))/denomZ : 0.0;	                        \
        mZ[0]/=nZ;               \
                                                                        \
        st[2*thermo.mixFracI_]=mZ[0];                                   \
        st[2*thermo.mixFracI_+1]=mZ[1];                                 \
                                                                        \
        forAllVariablesIn(thermo, var) /*note: mixture fraction has been removed from variable list*/  \
        {                                                               \
            norm = thermo.derivedNormFactors_[var.key()]->lookup(mZ);   \
                                                                        \
            scalar M = min                                              \
                (                                                       \
                    1.0,                                                \
                    thermo.firstMoment(var())_i_ * norm[0]              \
                );                                                      \
            st[2*var()] = M;                                            \
                                                                        \
            if (M>EPSILON)                                              \
            {                                                           \
                scalar m1=thermo.firstMoment(var())_i_;                 \
                scalar m2=thermo.secondMoment(var())()_i_;              \
                                                                        \
                scalar v=max( m2-sqr(m1), 0.0);                         \
                scalar S =                                              \
                    (                                                   \
                        norm[2] + v - 2.0 *                             \
                        sqrt(max(norm[2] * v, 0.0))                     \
                    ) * norm[1];                                        \
                scalar denom= ( M - M*M );                              \
                if (denom>EPSILON) S/=denom; else S=0.0;                \
                st[2*var()+1] = min(1.0, S);                            \
            } else st[2*var()+1]=0.0;                                   \
        }                                                               \
                                                                        \
        if (thermo.solvesEnthalpy())                                    \
        {                                                               \
            label i=thermo.preIntegratedTable().indexOfPV("h");         \
            st[i]=thermo.h()_i_;                                        \
        }                                                               \
    

    class SIJPDFvarZthermo;

    class SIJPDFvarZthermoIndexDriver
    {

    public:

        static inline scalar cellNormM1
        (
            const abstractProgressVariableThermo* th,
            label cellI,
            label pdfI
        ) 
            {
                SIJPDFVARZTHERMO_NORM_M1( pdfI, [cellI] );
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
                SIJPDFVARZTHERMO_NORM_M1( pdfI, .boundaryField()[patchI][faceI] );
                return normM1;
            }

        static inline List<scalar> cellTableIndex
        (
            const abstractProgressVariableThermo* th,
            label cellI
        ) 
            {
                SIJPDFVARZTHERMO_INDEX( [cellI] );
                return st;
            }

        static inline List<scalar> faceTableIndex
        (
            const abstractProgressVariableThermo* th,
            label patchI, 
            label faceI
        ) 
            {
                SIJPDFVARZTHERMO_INDEX( .boundaryField()[patchI][faceI] )
                return st;
            }
    };

}

#endif
