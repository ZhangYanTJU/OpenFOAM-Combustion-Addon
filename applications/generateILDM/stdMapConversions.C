using namespace std;

#include "stdMapConversions.H"
#include <climits>

ColumnVector toVector
(
    const Cantera::IdealGasMix& gas, 
    const Cantera::compositionMap& Y, 
    std::string leaveOut
)
{
    int size=gas.nSpecies();

    int li=INT_MAX;
    if (leaveOut!="") 
    {
        li=gas.speciesIndex(leaveOut);
        size--;
    }

    ColumnVector vY(size, 0.0);
    for (Cantera::compositionMap::const_iterator it=Y.begin();
         it!=Y.end(); it++)
    {
        int i=gas.speciesIndex(it->first);
        if (i!=li)
            vY( i>li ? i-1 : i )=it->second;
    }
    return vY;
}


Cantera::compositionMap toMap
(
    const Cantera::IdealGasMix& gas, 
    const ColumnVector& vY, 
    std::string leaveOut
)
{
    int li=INT_MAX;
    if (leaveOut!="") li=gas.speciesIndex(leaveOut);

    Cantera::compositionMap Y;
    for (int i=0;i<vY.length();i++)
        if (vY(i)!=0.0 && li!=i) Y[gas.speciesName(i)]
            = vY( i>li ? i-1 : i );
    return Y;
}

void toArray
(
    const Cantera::IdealGasMix& gas, 
    const Cantera::compositionMap& Y, 
    double* y
)
{
    for (Cantera::compositionMap::const_iterator it=Y.begin();
         it!=Y.end(); it++)
    {
        y[gas.speciesIndex(it->first)]=it->second;
    }
}

Cantera::compositionMap toMap
(
    const Cantera::IdealGasMix& gas, 
    double* y
)
{
    Cantera::compositionMap Y;
    for (int i=0;i<gas.nSpecies();i++)
        Y[gas.speciesName(i)] = y[i];
    return Y;
}

void write
(
    std::ofstream& f,
    const Cantera::IdealGasMix& gas,
    const Cantera::compositionMap& mY
)
{
    double y[gas.nSpecies()];
    toArray(gas, mY, y);
    for (int i=0;i<gas.nSpecies();i++)
        f<<" "<<y[i];
}
