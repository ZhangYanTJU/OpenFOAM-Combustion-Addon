#include "LDM.H"

#include "cantera/equilibrium.h"

#include "octave/dRowVector.h"

namespace Foam
{

double LDMconservedVars::h() const 
{ 
    if (!unburntTfixed_)
        return h_; 
    else
    {
        return sys_.setGasToMixFracWithProducts(z_, p_, 0, unburntT_);
    }
}

void LDMconservedVars::seth(double h) 
{ 
    if (unburntTfixed_)
        FatalErrorIn
            (
                "LDMconservedVars::seth()"
            )   << "Cannot change enthalpy because unburnt temperature is fixed to "<<unburntT_
                << endl
                << exit(FatalError);
    else
        h_=h; 
}

    
LDMconservedVars::LDMconservedVars
(
    const ChemicalSystem& sys,
    const HashTable<double, word>& sym
)
    : sys_(sys),
      z_(sys.stoichZ()),
      p_(1e5),
      h_(0.0),
      unburntTfixed_(false)
{
    HashTable<label,word> used_syms;

    if (sym.found("p")) 
    {
        p_=sym["p"];
        used_syms.insert("p", 1);
    }
    
    if (sym.found("z")) 
    {
        z_=sym["z"];
        used_syms.insert("z", 1);
    }
    else if (sym.found("phi"))
    {
        z_=1.0/(1.0+sys.stoichOF()/sym["phi"]);
        used_syms.insert("phi", 1);
    }
    else if (sym.found("lambda"))
    {
        z_=1.0/(1.0+sys.stoichOF()*sym["lambda"]);
        used_syms.insert("lambda", 1);
    }

    if (sym.found("h"))
    {
        h_=sym["h"];
        used_syms.insert("h", 1);
    }
    else if (sym.found("T"))
    {
        h_=sys.setGasToMixFracWithProducts(z_, p_, 0, sym["T"]);
        used_syms.insert("T", 1);
    }
    else if (sym.found("Tunburnt"))
    {
        unburntTfixed_=true;
        unburntT_=sym["Tunburnt"];
        h_=sys.setGasToMixFracWithProducts(z_, p_, 0, unburntT_);
        used_syms.insert("Tunburnt", 1);
    }

    HashTable<double, word> leftsym(sym);
    for (HashTable<double, word>::const_iterator it=sym.begin();
         it!=sym.end(); it++)
        if (used_syms.found(it.key())) leftsym.erase(it.key());

    if (leftsym.size()>0)
        FatalErrorIn
            (
                "LDMconservedVars::LDMconservedVars()"
            )   << "Error in specification of default conserved variables."
                << endl
                << "The following variables are redundant: "
                << leftsym.toc()
                << endl
                << exit(FatalError);
}




LDMconservedVars::LDMconservedVars(const ChemicalSystem& sys, double z, double p, double h)
    : sys_(sys),
      z_(z),
      p_(p),
      h_(h),
      unburntTfixed_(false),
      unburntT_(-1)
{
}

double LDMconservedVars::getByName(const word& name) const
{
    if (name=="z") return z();
    if (name=="p") return p();
    if (name=="h") return h();

    FatalErrorIn
        (
            "LDMconservedVars::getByName()"
        )   << "Unknown variable name " << name
            << endl
            << exit(FatalError);

    return -1;
}






defineTypeNameAndDebug(LDM, 0);
defineRunTimeSelectionTable(LDM, dictionary);
defineRunTimeSelectionTable(LDM, interpolate);

autoPtr<LDM> LDM::New
(
    const ChemicalSystem& s, 
    LDMconservedVars cv,
    dictionary& dict
)
{
    word type(dict.lookup("ILDMgeneratorType"));    

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "LDM::New()"
        )   << "Unknown LDM generator " << type
            << endl
            << "Valid generators are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<LDM>(cstrIter()(s, cv, dict));
}  

autoPtr<LDM> LDM::New
(
        const LDM& lowest,
        LDMconservedVars cv,
        const FinalState& finalState,
        dictionary& dict
)
{
    word type(dict.lookup("ILDMgeneratorType"));

    interpolateConstructorTable::iterator cstrIter =
        interpolateConstructorTablePtr_->find(type);

    if (cstrIter == interpolateConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "LDM::New()"
        )   << "Unknown LDM generator " << type
            << endl
            << "Valid generators are :" << endl
            << interpolateConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<LDM>(cstrIter()(lowest, cv, finalState));
}


void LDM::append(double pvval, const Cantera::compositionMap& compo)
{
    operator[](pvval)=toVector(gas(), compo);
}

void LDM::append(double pvval, const ColumnVector& compo)
{
    operator[](pvval)=compo;
}

string LDM::uniqueDescription() const
{
    std::ostringstream name;
    name<<"flame__h_"<<h()<<"kJkg__p_"<<p()<<"Pa__z_"<<z()<<".cache";
    return name.str();
}



LDM::LDM(const LDM& l)
    : std::map<double, ColumnVector>(l),
      sys_(l.sys_),
      pv_(l.pv_),
      pvidx_(l.pvidx_),
      cv_(l.cv_),
      nSlow_(l.nSlow_),
      P_(l.P_)
{
}


LDM::LDM
(
    const ChemicalSystem& s, 
    const std::string& pvname,
    LDMconservedVars cv,
    label resolution
)
    : sys_(s),
      pv_(pvname),
      pvidx_(s.gas().speciesIndex(pvname)),
      cv_(cv),
      nSlow_(gas().nElements()+1),
      P_(nSlow_, sys_.gas().nSpecies(), 0.0),
      resolution_(resolution)

{
    calcP();
}

LDM::LDM
(
    const ChemicalSystem& s, 
    LDMconservedVars cv,
    dictionary& dict
)
    : sys_(s),
      pv_(word(dict.lookup("progressVariableName"))),
      pvidx_(s.gas().speciesIndex(pv_)),
      cv_(cv),
      nSlow_(gas().nElements()+1),
      P_(nSlow_, sys_.gas().nSpecies(), 0.0),
      resolution_(readLabel(dict.lookup("resolution")))
{
    calcP();
}
 
void LDM::calcP()
{

    // element mass fractions
    const array_fp& Wa=gas().atomicWeights();
    for (int i=0; i<gas().nSpecies(); i++)
    {
        for (int j=0; j<gas().nElements(); j++)
            P_(j, i)=gas().nAtoms(i, j)*(Wa[j]/gas().molecularWeight(i));
    }
    // progress variable
    P_(nSlow_-1, gas().speciesIndex(pv_))=1.0;

}

double LDM::pveq() const
{
    if (!pveq_.valid())
    {
        Info<<"Calculating equilibrium value of progress variable"<<endl;
        
        sys_.setGasTo(z(), p(), h(), true);
        
        try 
        {
            equilibrate(gas(), "HP");
        }
        catch (CanteraError) 
        {
            showErrors(cout);
            {
                using namespace Foam;
                FatalErrorIn
                    (
                        "LDM::calcP()"
                    )
                    << "Could not solve for equilibrium composition." <<endl
                        << "Debug output from cantera displayed above"
                        << exit(FatalError);
            }
        }

        ColumnVector Yeq(gas().nSpecies(), 0.0);
        gas().getMassFractions(Yeq.fortran_vec());
        
        pveq_.reset(new double(P_.row(nSlow_-1) * Yeq));
    }

    return pveq_();
}

LDM::~LDM()
{
}

std::ofstream& operator<<(std::ofstream& f, const LDM& ildm)
{
    for (LDM::const_iterator it=ildm.begin();
         it!=ildm.end(); it++)
    {
        f<<it->first<<" ";
        for (int i=0; i<it->second.length(); i++)
         f<<it->second(i)<<" ";
        f<<std::endl;
    }
    return f;
}

std::ifstream& operator>>(std::ifstream& f, LDM& ildm)
{
    ildm.clear();
    while (!f.eof())
    {

        std::string linestr;
        getline(f, linestr);
        if (f.fail()) break;
        std::istringstream line(linestr);

        double pvval;
        line>>pvval;

        //ColumnVector Y(ildm.gas().nSpecies(), 0.0);
        std::vector<double> data;
        while (!line.eof())
        {
            double y;
            line>>y;
            if (line.fail()) break;
            data.push_back(y);
        }
        ColumnVector Y(data.size(), 0.0);
        for (int i=0; i< data.size(); i++)
            Y(i)=data[i];

        ildm[pvval]=Y;
    }
    return f;
}


void LDM::addToTable
(
    chemistryTable& tab,
    const LDMconservedVars& cv
)
{
    // get table index to current parameter set          
    int idx[tab.dimensionality()];
    for (int dim=0; dim<tab.dimensionality()-1; dim++)
    {
        idx[dim]=
            tab.indexAt
            (
                dim, 
                cv.getByName
                (
                    tab.nameOfPV(dim)
                )
            );
    }

    // add 1D LDM to table
    if (size() > 1)
    {
        Spline ldm_interpolator(*this);

        label pvdim=tab.dimensionality()-1;
        for (
            idx[pvdim]=0;
            idx[pvdim]<tab.nElements(pvdim);
            idx[pvdim]++
        )
        {
            scalar pv= tab.valueAt(pvdim, idx[pvdim]) * pveq();
            Info<<"adding to table: pv="<<pv<<endl;
            
            tab.access(idx)=
                extendToTableEntry
                (
                    tab,
                    ldm_interpolator.interpolate(pv)
                );
            
        }
    }
    else
    {
        label pvdim=tab.dimensionality()-1;
        for (
            idx[pvdim]=0;
            idx[pvdim]<tab.nElements(pvdim);
            idx[pvdim]++
        )
        {
            scalar pv= tab.valueAt(pvdim, idx[pvdim]) * pveq();
            Info<<"adding to table: pv="<<pv<<endl;
            
            // always insert the same state
            tab.access(idx)=
                extendToTableEntry
                (
                    tab,
                    this->begin()->second
                );
            
        }
    }

}


void LDM::appendTableDescriptions
(
    ListDescription& pvdesc,
    ListDescription& contentdesc
)
{
    pvdesc.appendEntry(pv_);

    contentdesc.appendEntry("Ydot"&pv_);
    
    contentdesc.appendEntry("T");
    contentdesc.appendEntry("rho");
    contentdesc.appendEntry("lambda");
    contentdesc.appendEntry("mu");
    contentdesc.appendEntry("Cp");

    contentdesc.appendEntry("Y"&pv_);
}

void LDM::setupTable
(
    ListDescription& pvdesc,
    ListDescription& contentdesc,
    label* num,
    double* start,
    double* delta
)
{
    label i=pvdesc[pv_];
    num[i]=resolution_;
    start[i]=0.0;
    delta[i]=1.0/scalar(resolution_-1);
}


chemicalSystemState LDM::extendToTableEntry
(
    const chemistryTable& tab,
    const ColumnVector& Y
)
{
    const ListDescription& desc=tab.contentDescription();
    chemicalSystemState s(desc);

    sys_.setGasTo(Y.data(), p(), h());

    s[desc["Y"&pv_]]=Y(pvidx_);
    s[desc["T"]]=gas().temperature();
    s[desc["rho"]]=gas().density();
    s[desc["Cp"]]=gas().cp_mass();;
    s[desc["mu"]]=sys_.tr().viscosity();
    s[desc["lambda"]]=sys_.tr().thermalConductivity();

    int nsp=gas().nSpecies();

    doublereal DcDt[nsp], c[nsp];
    gas().getNetProductionRates(DcDt);
    gas().getConcentrations(c);

    double rho=gas().density();

    double DrhoDt=0.0;
    for (int i=0; i<nsp; i++) 
        DrhoDt+=gas().molarMass(i)*DcDt[i];

    s[desc["Ydot"&pv_]]=gas().molarMass(pvidx_)*
        (rho*DcDt[pvidx_]-c[pvidx_]*DrhoDt)/(rho*rho);
    
    return s;
}


}
