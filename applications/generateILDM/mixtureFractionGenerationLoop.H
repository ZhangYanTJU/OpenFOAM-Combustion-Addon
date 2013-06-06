#ifndef mixtureFractionGenerationLoop_H
#define mixtureFractionGenerationLoop_H

#include "generationLoop.H"

namespace Foam
{

class mixtureFractionGenerationLoop
: public generationLoop
{

public:
    TypeName("z");

    mixtureFractionGenerationLoop(const mixtureFractionGenerationLoop&);
    mixtureFractionGenerationLoop(const ChemicalSystem&, dictionary& dict);

    autoPtr<LDM::FinalState> finalLowerState(const LDMconservedVars& cv) const;
    autoPtr<LDM::FinalState> finalUpperState(const LDMconservedVars& cv) const;

    virtual double getBestConvergenceValue() const;
    virtual void setParameter(LDMconservedVars& cv, double value) const;
    virtual double getParameter(LDMconservedVars& cv) const;
    virtual autoPtr<generationLoop> clone() const;
    
    void writeAdditionalInformationToTableFile(Ostream& f);

};

}

#endif
