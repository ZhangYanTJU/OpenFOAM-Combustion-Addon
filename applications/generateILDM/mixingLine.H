#ifndef MIXINGLINE_H
#define MIXINGLINE_H

#include "LDM.H"

namespace Foam
{

class mixingLine
: public LDM
{

    bool convert_;

public:
    TypeName("mixingLine");

    mixingLine
    (
        const ChemicalSystem& s, 
        const std::string& pvname,
        LDMconservedVars cv,
        bool convert=false
    );

    mixingLine
    (
        const ChemicalSystem& s, 
        LDMconservedVars cv,
        dictionary& dict,
        bool convert=false
    );

    mixingLine(const mixingLine&);

    virtual ~mixingLine();

    void calculate();
};

}

#endif
