#ifndef GENERATIONLOOP_H
#define GENERATIONLOOP_H


#include <map>
#include <vector>
#include <string>

#include "scalar.H"
#include "label.H"
#include "word.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "Istream.H"
#include "HashTable.H"
#include "chemistryTable.H"

#include "ChemicalSystem.H"
#include "LDM.H"
#include "interpolatedLDM.H"

namespace Foam
{


class generationLoop
{

protected:
    const ChemicalSystem& sys_;

    scalar start_;
    scalar end_;
    label resolution_;

    typedef std::vector<double> valueList;
    valueList upwardValueList_;
    valueList downwardValueList_;

    typedef std::vector<generationLoop*> loopList;
    static loopList allLoops_;

    loopList::iterator myIndex();
    void registerMe();
    void unregisterMe();
    
  public:

    inline const ChemicalSystem& sys() const { return sys_; }

    //- Runtime type information
    TypeName("generationLoop");
    
    // Declare run-time constructor selection table
    
    declareRunTimeSelectionTable
    (
        autoPtr,
        generationLoop,
        dictionary,
        (const ChemicalSystem& sys, dictionary& dict),
        (sys, dict)
    );


    // Selectors

    //- Return a reference to the selected surfaceWriter
    static autoPtr<generationLoop> New(const ChemicalSystem& sys, dictionary& dict);


    enum loopDirection {up,down};

    class iNew
    {
        const ChemicalSystem& sys_;

      public:

        iNew(const ChemicalSystem& sys)
            : sys_(sys)
        {}

        autoPtr<generationLoop> operator()(Istream& is) const
        {
            dictionary d(is);
            return generationLoop::New(sys_, d);
        }
    };




    generationLoop(const generationLoop&);
    generationLoop(const ChemicalSystem&, dictionary& dict);
    virtual ~generationLoop();

    virtual void populate();

    inline const word& varName() const
    {
        return type();
    }

    inline label resolution() const
        { return resolution_; }
    
    inline scalar start() const
        { return start_; }
    
    inline scalar end() const
        { return end_; }


    virtual double getBestConvergenceValue() const =0;
    virtual void setParameter(LDMconservedVars& p, double value) const =0;
    virtual double getParameter(LDMconservedVars& p) const =0;
    virtual autoPtr<generationLoop> clone() const =0;

    virtual autoPtr<LDM::FinalState> finalLowerState(const LDMconservedVars& p) const;
    virtual autoPtr<LDM::FinalState> finalUpperState(const LDMconservedVars& p) const;

    static autoPtr<LDM> performCalculation
    (
        const ChemicalSystem& sys,
        IOdictionary& generateIldmDict,
        LDMconservedVars& p,        
        chemistryTable& tab,
        bool& interpolated,
        autoPtr<LDM> lastsuccessful,
        const autoPtr<LDM::FinalState>&
    );

    void performLoop
    (
        IOdictionary& generateIldmDict,
        LDMconservedVars&, 
        chemistryTable&
    );

    virtual void appendTableDescriptions
    (
        ListDescription& pvdesc,
        ListDescription& contentdesc
    );

    virtual void setupTable
    (
        ListDescription& pvdesc,
        ListDescription& contentdesc,
        label* num,
        double* start,
        double* delta
    );


    static void generateTable
    (
        const ChemicalSystem& sys,
        IOdictionary& generateIldmDict,
        LDMconservedVars&, 
        chemistryTable&
    );

    virtual void writeAdditionalInformationToTableFile(Ostream& f);
};

}

#endif
