
#include "FGM.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "OFstream.H"

#include <cantera/Cantera.h>
#include <cantera/IdealGasMix.h>
#include <cantera/onedim.h>
#include <cantera/equilibrium.h>
#include <cantera/transport.h>

namespace Foam
{

defineTypeNameAndDebug(FGM, 0);
addToRunTimeSelectionTable(LDM, FGM, dictionary);
addToRunTimeSelectionTable(LDM, FGM, interpolate);

FGM::FGM
(
    const ChemicalSystem& s,
    const std::string& pv,
    const std::string& solverName,
    LDMconservedVars cv,
    int resolution
)
    : LDM(s, pv, cv, resolution),
      initialResolution_(5),
      domainLength_(0.1),
      onlyCache_(false)
{
}


FGM::FGM
(
    const ChemicalSystem& s, 
    LDMconservedVars cv,
    dictionary& dict
)
    : LDM(s, cv, dict),
      initialResolution_(readLabel(dict.lookup("laminarFlameInitialDomainResolution"))),
      domainLength_(readScalar(dict.lookup("laminarFlameDomainLength"))),
      onlyCache_(false)
{
    if (dict.found("onlyFromCache"))
      onlyCache_=Switch(dict.lookup("onlyFromCache"));
}

    FGM::FGM
    (
        const LDM& lowest,
        LDMconservedVars cv,
        const FinalState& finalState
    ) : LDM(lowest, cv, finalState),
       onlyCache_(true)
   {
   }




void FGM::calculate()
{

    std::ifstream cache(uniqueDescription().c_str());

    if (!cache.fail())
    {
        //read results from cache
        cout << "Reading from cache file " << uniqueDescription() << endl;
        cache >> (*this);
        return;
    }
    else
    {        

	if (onlyCache_) throw std::string("only from cache: calculation forbidden by user");

        int i;
    
        sys().setGasTo(z(), p(), h(), true);

        double temp=gas().temperature();
        double pressure=p();
  
        doublereal rho_in=gas().density();
    
        int nsp = gas().nSpecies();
      
        double *yin=new double[nsp];
        gas().getMassFractions(yin);
    
        double *x=new double[nsp];
        gas().getMoleFractions(x);
    
        equilibrate(gas(),"HP");

        double *yout=new double[nsp];
        gas().getMassFractions(yout);
        doublereal rho_out = gas().density();
        doublereal Tad=gas().temperature();
   
        cout<<" adiabatic flame temperatur T_ad="<<Tad<<endl; 

        //=============  build each domain ========================


        //-------- step 1: create the flow -------------

        FreeFlame flow(&gas());

        int nz=initialResolution_; //50;//5;
        doublereal lz=domainLength_; //0.12;
        //doublereal shift=flameshift;

        double x0=0, x1=lz, xm=0.5*lz;

        // create an initial grid
        double uin=0.3;
      
        doublereal *z=new double[nz+1];
        doublereal dz=lz/(doublereal(nz-1));
        for(int iz=0;iz<nz;iz++)
        {
            z[iz]=(doublereal(iz))*dz;
        }
        
        //add one node onto end of domain to help with zero gradient at outlet
        z[nz]=lz*1.05;
        nz++;
        
        flow.setupGrid(nz, z);      
      
        // specify the objects to use to compute kinetic rates and 
        // transport properties
   
        flow.setTransport(sys().tr());
        flow.setKinetics(gas());
        flow.setPressure(pressure);
      
        //------- step 2: create the inlet  -----------------------
      
        Inlet1D inlet;
      
        inlet.setMoleFractions(DATA_PTR(x));
        doublereal mdot=uin*rho_in;
        inlet.setMdot(mdot);
        inlet.setTemperature(temp);


        //------- step 3: create the outlet  ---------------------

        Outlet1D outlet;

        //=================== create the container and insert the domains =====

        std::vector<Domain1D*> domains;
        domains.push_back(&inlet);
        domains.push_back(&flow);
        domains.push_back(&outlet);

        //    OneDim flamesim(domains);

        Sim1D flame(domains);

        int steps[]={1,2,5,10,20,40};
        flame.setTimeStep(1e-5, 6, steps);

        flame.setMaxTimeStep(100.0);

        //----------- Supply initial guess----------------------

        vector_fp locs;
        vector_fp value;
          
        locs.resize(3);
        value.resize(3);
            
        //ramp values from inlet to adiabatic flame conditions 
        //  over 70% of domain and then level off at equilibrium
        double z1=0.7;

        double uout;
        uout=inlet.mdot()/rho_out;
        uin=inlet.mdot()/rho_in;
        locs[0]=0.0; locs[1]=z1; locs[2]=1.0;
        value[0]=uin; value[1]=uout; value[2]=uout;
        flame.setInitialGuess("u",locs,value);
        
        value[0]=temp; value[1]=Tad; value[2]=Tad;
        flame.setInitialGuess("T",locs,value);
        
        for(i=0;i<nsp;i++)
        {
            value[0]=yin[i]; value[1]=yout[i]; value[2]=yout[i];
            flame.setInitialGuess(gas().speciesName(i), locs, value);
        }
      
        inlet.setMoleFractions(DATA_PTR(x));
        inlet.setMdot(mdot);
        inlet.setTemperature(temp);
      
        flame.showSolution();

        int flowdomain=1;
        double ratio=10.0;
        double slope=0.2;
        double curve=0.02;
        double prune=-0.00005;
        flame.setRefineCriteria(flowdomain,ratio,slope,curve,prune);

        int loglevel=1;
        bool refine_grid = true;
        flow.fixTemperature();
        flame.setFixedTemperature(0.5*Tad);

        /* Solve species*/             
        refine_grid=false;
        //    flame.setAdiabaticFlame();
        flame.solve(loglevel,refine_grid);
      
        ratio=3.0;
        slope=0.1;
        curve=0.2;
        flame.setRefineCriteria(flowdomain,ratio,slope,curve);
        refine_grid = true;
        flow.solveEnergyEqn();
        flame.solve(loglevel,refine_grid);

        int np=flow.nPoints();
            
        for (int i=0; i<np; i++)
        {

            // 2 additional entries for gradient and x
            ColumnVector Y(nsp+2, 0.0);

            for (int j=0;j<nsp;j++)
            {
                Y(j)=flame.value(flowdomain, 4+j, i);
                Y(nsp+1)=flow.grid(i);
            }

            (*this)[Yp(Y)]=Y;
        }

        // calculate gradient of progress variable
        {
            // forward difference at first point
            ColumnVector& Y0 = begin()->second;
            ColumnVector& Y1 = (++begin()) ->second;
            Y0(nsp)=(Y1(pvidx_)-Y0(pvidx_))/(Y1(nsp+1)-Y0(nsp+1));
        }
        {
            // backward difference at last point
            ColumnVector& Y0=( --   end() )   ->second;
            ColumnVector& Y1=( --(--end() ) ) ->second;
            Y1(nsp)=(Y1(pvidx_)-Y0(pvidx_))/(Y1(nsp+1)-Y0(nsp+1));
        }
        for (iterator p=(++begin()); p!=--(--end()); p++)
        {
            // central difference in between
            ColumnVector& Y0 = (--p) -> second;
            ColumnVector& Y1 = (++p) -> second;
            p->second(nsp)=(Y1(pvidx_)-Y0(pvidx_))/(Y1(nsp+1)-Y0(nsp+1));
        }
      
        cout << "Writing to cache file " << uniqueDescription() << endl;
        std::ofstream out(uniqueDescription().c_str());
        out << (*this);
      
    }
}

void FGM::appendTableDescriptions
(
    ListDescription& pvdesc,
    ListDescription& contentdesc
)
{
    LDM::appendTableDescriptions(pvdesc, contentdesc);

    contentdesc.appendEntry("grad("&pv_&")");
    contentdesc.appendEntry("x");
}


chemicalSystemState FGM::extendToTableEntry
(
    const chemistryTable& tab,
    const ColumnVector& Y
)
{
    chemicalSystemState s = LDM::extendToTableEntry(tab, Y);

    const ListDescription& desc=tab.contentDescription();

    s[desc["grad("&pv_&")"]] = Y(gas().nSpecies()); 
    s[desc["x"]] = Y(gas().nSpecies()+1);
    
    return s;
}

}
