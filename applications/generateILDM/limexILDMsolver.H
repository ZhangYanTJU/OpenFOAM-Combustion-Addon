#ifndef limexILDMsolver_H
#define limexILDMsolver_H

#include "ILDMsolver.H"

namespace Foam
{

class limexILDMsolver
: public ILDMsolver
{

    ColumnVector objectiveF(const ColumnVector& Y, Matrix& B);

    // static
    static limexILDMsolver* csol;
    static void objectiveF_LIMEXadaptor
    (
        int *n,
        int *nz,
        double *t,
        double *y,
        double *f,
        double *b,
        int *ir,
        int *ic,
        int *FcInfo
    );

    double pvval_;

public:
    TypeName("limex");

    limexILDMsolver(ILDM&);
    virtual ~limexILDMsolver();

    virtual Cantera::compositionMap solve
    (
        double pv,
        const Cantera::compositionMap& inital
    );
};

}

#endif
