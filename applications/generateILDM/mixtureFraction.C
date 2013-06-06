
#include "mixtureFraction.H"

// i is index of atom on which mixture fractions definition is based
double Z(int i, Cantera::IdealGasMix& gas, const Cantera::compositionMap& Y)
{
  double z=0;
  for (Cantera::compositionMap::const_iterator it=Y.begin();
       it!=Y.end(); it++)
  {
    int j=gas.speciesIndex(it->first);
    z+=gas.nAtoms(j, i) * (gas.atomicWeight(i)/gas.molecularWeight(j)) * it->second;
    /*
    cout<<it->first<<": "
        <<gas.nAtoms(j, i)<<" "
        <<gas.atomicWeight(i)<<" "
        <<gas.molecularWeight(j)<<" "
        <<it->second<<endl;
    */
  }
  return z;
}

