#ifndef octaveILDMsolver_H
#define octaveILDMsolver_H

#include "ILDMsolver.H"

namespace Foam
{

class octaveILDMsolver
: public ILDMsolver
{

    ColumnVector objectiveF(const ColumnVector& Y);

    // static
    static octaveILDMsolver* csol;
    static ColumnVector objectiveF_adaptor(const ColumnVector& Y);

    double pvval_;

public:
    TypeName("octave");

    octaveILDMsolver(ILDM&);
    virtual ~octaveILDMsolver();

    virtual Cantera::compositionMap solve
    (
        double pv,
        const Cantera::compositionMap& inital
    );
};

}

#endif
