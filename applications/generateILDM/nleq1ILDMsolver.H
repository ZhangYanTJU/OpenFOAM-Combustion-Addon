#ifndef nleq1ILDMsolver_H
#define nleq1ILDMsolver_H

#include "ILDMsolver.H"

namespace Foam
{


class nleq1ILDMsolver
: public ILDMsolver
{

    double pvval_;

public:
    TypeName("nleq1");

    nleq1ILDMsolver(ILDM&);
    virtual ~nleq1ILDMsolver();

    virtual ColumnVector solve
    (
        double pv,
        const ColumnVector& inital
    );

    ColumnVector objective(const ColumnVector& Y);
    Matrix jacobian(const ColumnVector& Y);

};

}

#endif
