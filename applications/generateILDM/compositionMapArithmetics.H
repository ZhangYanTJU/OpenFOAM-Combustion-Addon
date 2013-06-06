#ifndef COMPOSITIONMAPARITHMETICS_H
#define COMPOSITIONMAPARITHMETICS_H

#include <map>
#include <string>

namespace Cantera
{
    typedef std::map<std::string, double> compositionMap;
}


Cantera::compositionMap operator*(double c, const Cantera::compositionMap& Y);

Cantera::compositionMap operator+
    (
      const Cantera::compositionMap& Y1,
      const Cantera::compositionMap& Y2
    );


#endif
