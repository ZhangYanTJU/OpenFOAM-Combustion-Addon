#include "generationLoop.H"
#include "mixingLine.H"

using namespace Foam;


defineTypeNameAndDebug(generationLoop, 0);
defineRunTimeSelectionTable(generationLoop, dictionary);



autoPtr<generationLoop> generationLoop::New
(
    const ChemicalSystem& sys,
    dictionary& dict
)
{
    word variableName(dict.lookup("variable"));    

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(variableName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "generationLoop::New(dictionary& dict)"
        )   << "Unknown loop variable " << variableName
            << endl
            << "Valid loop variables are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<generationLoop>(cstrIter()(sys, dict));

}

generationLoop::loopList::iterator generationLoop::myIndex()
{
    loopList::iterator it=allLoops_.begin();
  
    for (; it!=allLoops_.end(); it++)
        if ((*it)==this) return it;

    return allLoops_.end();
}

void generationLoop::registerMe()
{
    allLoops_.push_back(this);
}

void generationLoop::unregisterMe()
{
    loopList::iterator i=myIndex();
    if (i!=allLoops_.end()) allLoops_.erase(i);
}

generationLoop::loopList generationLoop::allLoops_;

generationLoop::generationLoop(const generationLoop& o)
    : sys_(o.sys_),
      start_(o.start_),
      end_(o.end_),
      resolution_(resolution_),
      upwardValueList_(o.upwardValueList_),
      downwardValueList_(o.downwardValueList_)
{
    registerMe();
}
    
generationLoop::generationLoop(const ChemicalSystem& sys, dictionary& dict)
    : sys_(sys),
      start_(readScalar(dict.lookup("start"))),
      end_(readScalar(dict.lookup("end"))),
      resolution_(readLabel(dict.lookup("resolution")))
{
    //Info<<"Loop over "<<variable()<<" defined."<<endl;

    registerMe();
}

void generationLoop::populate()
{

    // ======= populate value list ==========
    if (start_>end_)
    {
        double tmp=start_;
        start_=end_;
        end_=tmp;
    }

    double sv=getBestConvergenceValue();
    int svi=0;
  
    double delta=(end_-start_)/scalar(resolution_-1);
    double min=1e100;
    valueList tmpvlist_;
    for (int i=0; i<resolution_; i++)
    {
        double v=start_ + double(i)*delta;
        if (mag(v-sv)<min) { min=mag(v-sv); svi=i; }
        tmpvlist_.push_back(v);
    }

    for (int i=svi; i<resolution_; i++)
    {
        upwardValueList_.push_back(tmpvlist_[i]);
    }
    for (int i=svi-1; i>=0; i--)
    {
        downwardValueList_.push_back(tmpvlist_[i]);
    }

    /*
        for (int i=0;i<resolution_;i++)
        Info<<tmpvlist_[i]<<endl;
        */

}

generationLoop::~generationLoop()
{
  unregisterMe();
}


autoPtr<LDM> generationLoop::performCalculation
(
    const ChemicalSystem& sys,
    IOdictionary& generateIldmDict,
    LDMconservedVars& cv,
    chemistryTable& tab,
    bool& interpolated,
    autoPtr<LDM> lastsuccessful,
    const autoPtr<LDM::FinalState>& final
)
{    

    Info<<"Caluclating ILDM with"
        <<" z="<<cv.z()
        <<" p="<<cv.p()
        <<" h="<<cv.h()
        <<endl;

    dictionary dict(generateIldmDict);

    if (!interpolated)
    {
        try
        {
            autoPtr<LDM> ldm = LDM::New( sys, cv, dict );
            
            ldm().calculate();
            
            ldm().addToTable( tab, cv );
            
            return ldm;
        }
        catch(CanteraError)
        {
            showErrors(cout);
        }
        catch(std::string s)
        {
            cout<<s<<std::endl;
        }
        catch(...)
        {
        }
    }

    Info<<"Calculation failed. Trying interpolation..."<<endl;

    autoPtr<LDM> ldm = LDM::New( lastsuccessful, cv, final(), dict );
    
    ldm().addToTable( tab, cv );
    
    interpolated=true;
    return lastsuccessful;

}

autoPtr<LDM::FinalState> generationLoop::finalLowerState(const LDMconservedVars& cv) const
{
    return autoPtr<LDM::FinalState>();
}


autoPtr<LDM::FinalState> generationLoop::finalUpperState(const LDMconservedVars& cv) const
{
    return autoPtr<LDM::FinalState>();
}



void generationLoop::performLoop
(
    IOdictionary& generateIldmDict,
    LDMconservedVars& cv, 
    chemistryTable& tab
)
{
    
    loopList::iterator i=myIndex();
    i++;
    
    if (i!=allLoops_.end())
        (*i)->performLoop(generateIldmDict, cv, tab);
    else
    {

            #warning !!FIX: conversion of fuel into products at low enthalpy is always on
        mixingLine initialLDM( sys_, cv, generateIldmDict, true );
        initialLDM.calculate();

        autoPtr<LDM> lastsuccessful;

        lastsuccessful.reset(new mixingLine(initialLDM));
        bool interpolated=false;
        Info<<"Looping upward"<<endl;
        for (valueList::const_iterator curval=upwardValueList_.begin();
             curval!=upwardValueList_.end(); curval++)
        {
            setParameter(cv, *curval);
            lastsuccessful=performCalculation
                (
                    sys_,
                    generateIldmDict, cv, tab,
                    interpolated, lastsuccessful,
                    finalUpperState(cv)
                );
        }

        lastsuccessful.reset(new mixingLine(initialLDM));
        interpolated=false;
        Info<<"Looping downward"<<endl;
        for (valueList::const_iterator curval=downwardValueList_.begin();
             curval!=downwardValueList_.end(); curval++)
        {
            setParameter(cv, *curval);
            lastsuccessful=performCalculation
                (
                    sys_,
                    generateIldmDict, cv, tab,
                    interpolated, lastsuccessful,
                    finalLowerState(cv)
                );
        }
    }
}



void generationLoop::generateTable
(
    const ChemicalSystem& sys,
    IOdictionary& generateIldmDict,
    LDMconservedVars& cv, 
    chemistryTable& tab
)
{
    if (allLoops_.size()>0) 
        allLoops_[0]->performLoop(generateIldmDict, cv, tab);
    else
    {
        bool interpolated=false;
        performCalculation
            (
                sys,
                generateIldmDict,
                cv,
                tab,
                interpolated,
                autoPtr<LDM>(),
                autoPtr<LDM::FinalState>()
            );
    }
}


void generationLoop::appendTableDescriptions
(
    ListDescription& pvdesc,
    ListDescription& contentdesc
)
{
    pvdesc.appendEntry(varName());
}

void generationLoop::setupTable
(
    ListDescription& pvdesc,
    ListDescription& contentdesc,
    label* n,
    double* s,
    double* d
)
{
    label i=pvdesc[varName()];

    n[i]=resolution();
    s[i]=start();
    d[i]=(end()-start()) / scalar(resolution()-1);
}


void generationLoop::writeAdditionalInformationToTableFile(Ostream& f)
{
}
