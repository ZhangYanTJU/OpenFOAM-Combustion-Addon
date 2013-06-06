#ifndef enthalpyGenerationLoop_H
#define enthalpyGenerationLoop_H

#include "generationLoop.H"

namespace Foam
{

class enthalpyGenerationLoop
: public generationLoop
{

public:
    TypeName("h");

    enthalpyGenerationLoop(const enthalpyGenerationLoop&);
    enthalpyGenerationLoop(const ChemicalSystem&, dictionary& dict);

    autoPtr<LDM::FinalState> finalLowerState(const LDMconservedVars& cv) const;
 
    virtual double getBestConvergenceValue() const;
    virtual void setParameter(LDMconservedVars& cv, double value) const;
    virtual double getParameter(LDMconservedVars& cv) const;
    virtual autoPtr<generationLoop> clone() const;
    

};

}

#endif
