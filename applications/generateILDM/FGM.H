#ifndef FGM_H
#define FGM_H

#include "LDM.H"
#include "Switch.H"

namespace Foam
{

class FGM
: public LDM
{

    label initialResolution_;
    scalar domainLength_;
    Switch onlyCache_;

public:
    TypeName("FGM");

    FGM
    (
        const ChemicalSystem& s, 
        const std::string& pv,
        const std::string& solverName,
        LDMconservedVars cv,
        int resolution=20
    );

    FGM
    (
        const ChemicalSystem& s, 
        LDMconservedVars cv,
        dictionary& dict
    );

    FGM
    (
        const LDM& lowest,
        LDMconservedVars cv,
        const FinalState& finalState
    );


    void calculate();

    virtual void appendTableDescriptions
    (
        ListDescription& pvdesc,
        ListDescription& contentdesc
    );

    chemicalSystemState extendToTableEntry
    (
        const chemistryTable& tab,
        const ColumnVector& Y
    );
};

}


#endif
