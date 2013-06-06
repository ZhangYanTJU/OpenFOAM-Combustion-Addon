
#include "compositionMapArithmetics.H"

Cantera::compositionMap operator*(double c, const Cantera::compositionMap& Y)
{
  Cantera::compositionMap Yn(Y);
  for (Cantera::compositionMap::iterator it=Yn.begin();
       it!=Yn.end();it++)
    it->second*=c;
  return Yn;
}

Cantera::compositionMap operator+
    (
      const Cantera::compositionMap& Y1,
      const Cantera::compositionMap& Y2
    )
{
  Cantera::compositionMap Yn(Y1);
  for (Cantera::compositionMap::const_iterator it=Y2.begin();
       it!=Y2.end();it++)
  {
    if (Yn.find(it->first)==Yn.end())
      Yn[it->first]=it->second;
    else
      Yn[it->first]+=it->second;
  }
  return Yn;
}

