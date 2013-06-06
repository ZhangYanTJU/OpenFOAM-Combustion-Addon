#include "enthalpyGenerationLoop.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

    defineTypeNameAndDebug(enthalpyGenerationLoop, 0);
    addToRunTimeSelectionTable(generationLoop, enthalpyGenerationLoop, dictionary);

    enthalpyGenerationLoop::enthalpyGenerationLoop(const enthalpyGenerationLoop& o)
        : generationLoop(*this)
    {
    }

    enthalpyGenerationLoop::enthalpyGenerationLoop
    (
        const ChemicalSystem& sys,
        dictionary& dict
    )
        : generationLoop(sys, dict)
    {
    }

    double enthalpyGenerationLoop::getBestConvergenceValue() const
    {
        return 0.0;
    }

    void enthalpyGenerationLoop::setParameter(LDMconservedVars& p, double value) const
    {
        p.seth(value);
    }

    double enthalpyGenerationLoop::getParameter(LDMconservedVars& p) const
    {
        return p.h();
    }

    autoPtr<LDM::FinalState> enthalpyGenerationLoop::finalLowerState
     (const LDMconservedVars& cv) const
    {
        LDMconservedVars fcv(cv);
        fcv.seth(sys().hzero(cv.z(), cv.p()));
        return autoPtr<LDM::FinalState>(new LDM::FinalState(sys(), fcv, true));
    }    

    autoPtr<generationLoop> enthalpyGenerationLoop::clone() const
    {
        return autoPtr<generationLoop>(new enthalpyGenerationLoop(*this));
    }

};
