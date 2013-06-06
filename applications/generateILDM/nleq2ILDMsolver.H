#ifndef nleq2ILDMsolver_H
#define nleq2ILDMsolver_H

#include "ILDMsolver.H"

namespace Foam
{

class nleq2ILDMsolver
: public ILDMsolver
{

    ColumnVector objectiveF(const ColumnVector& Y, Matrix& B);

    // static
    static nleq2ILDMsolver* csol;
    static void objectiveF_NLEQ2adaptor
    (
        long int *n,
        long int *nz,
        double *t,
        double *y,
        double *f,
        double *b,
        long int *ir,
        long int *ic,
        long int *FcInfo
    );

    double pvval_;

public:
    TypeName("nleq2");

    nleq2ILDMsolver(ILDM&);
    virtual ~nleq2ILDMsolver();

    virtual Cantera::compositionMap solve
    (
        double pv,
        const Cantera::compositionMap& inital
    );
};

}

#endif
