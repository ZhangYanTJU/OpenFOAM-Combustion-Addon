#ifndef INTERPOLATEDLDM_H
#define INTERPOLATEDLDM_H

#include "LDM.H"

namespace Foam
{

class interpolatedLDM
: public LDM
{

public:
    TypeName("interpolatedLDM");

    struct FinalState
    {
        LDMconservedVars cv;
        ColumnVector Y;

        FinalState(const ChemicalSystem& sys, LDMconservedVars fcv, bool unburnt=false);
    };

    interpolatedLDM
    (
        const LDM& lowest,
        LDMconservedVars cv,
        const FinalState& finalState
    );

    virtual ~interpolatedLDM();

    void calculate();
};

}

#endif
