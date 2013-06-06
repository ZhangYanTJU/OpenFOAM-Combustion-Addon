#include "mixingLine.H"
#include "cantera/equilibrium.h"

namespace Foam
{

defineTypeNameAndDebug(mixingLine, 0);

mixingLine::mixingLine
(
    const ChemicalSystem& s, 
    const std::string& pvname,
    LDMconservedVars cv,
    bool convert
)
    : LDM(s, pvname, cv, 2),
      convert_(convert)
{
}
    
mixingLine::mixingLine
(
    const ChemicalSystem& s, 
    LDMconservedVars cv,
    dictionary& dict,
    bool convert
)
    : LDM(s, cv, dict),
      convert_(convert)
{
}

mixingLine::mixingLine(const mixingLine&o)
    : LDM(o),
      convert_(o.convert_)
{
}
    
mixingLine::~mixingLine()
{
}

void mixingLine::calculate()
{
    ColumnVector Y(gas().nSpecies(), 0.0);

    sys().setGasTo(z(), p(), h(), convert_);
    gas().getMassFractions(Y.fortran_vec()); // unburnt composition
    (*this)[Yp(Y)]=Y;

    equilibrate(gas(), "HP");
    Info<<"adiabatic flame T="<<gas().temperature()<<endl;

    gas().getMassFractions(Y.fortran_vec()); // Y: equilibrium composition
    (*this)[Yp(Y)]=Y;
}

}
