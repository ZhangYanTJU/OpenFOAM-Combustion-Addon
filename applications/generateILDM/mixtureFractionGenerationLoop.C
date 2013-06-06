#include "mixtureFractionGenerationLoop.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

    defineTypeNameAndDebug(mixtureFractionGenerationLoop, 0);
    addToRunTimeSelectionTable(generationLoop, mixtureFractionGenerationLoop, dictionary);

    mixtureFractionGenerationLoop::mixtureFractionGenerationLoop(const mixtureFractionGenerationLoop& o)
        : generationLoop(*this)
    {
    }

    mixtureFractionGenerationLoop::mixtureFractionGenerationLoop
    (
        const ChemicalSystem& sys,
        dictionary& dict
    )
        : generationLoop(sys, dict)
    {
        if (mag(start_)>SMALL)
            FatalErrorIn("mixtureFractionGenerationLoop::mixtureFractionGenerationLoop")
                << "mixture fraction loop has to start at 0"
                    << endl << abort(FatalError);
    }

    autoPtr<LDM::FinalState> mixtureFractionGenerationLoop::finalLowerState
     (const LDMconservedVars& cv) const
    {
        LDMconservedVars fcv(cv);
        fcv.setz(0.0);
        return autoPtr<LDM::FinalState>(new LDM::FinalState(sys(), fcv));
    }

    autoPtr<LDM::FinalState> mixtureFractionGenerationLoop::finalUpperState
     (const LDMconservedVars& cv) const
    {
        LDMconservedVars fcv(cv);
        fcv.setz(1.0);
        return autoPtr<LDM::FinalState>(new LDM::FinalState(sys(), fcv));
    }


    double mixtureFractionGenerationLoop::getBestConvergenceValue() const
    {
        return sys().stoichZ();
    }

    void mixtureFractionGenerationLoop::setParameter(LDMconservedVars& p, double value) const
    {
        p.setz(value);
    }

    double mixtureFractionGenerationLoop::getParameter(LDMconservedVars& p) const
    {
        return p.z();
    }
    

    autoPtr<generationLoop> mixtureFractionGenerationLoop::clone() const
    {
        return autoPtr<generationLoop>(new mixtureFractionGenerationLoop(*this));
    }

    void mixtureFractionGenerationLoop::writeAdditionalInformationToTableFile(Ostream& f)
    {
    }

};
