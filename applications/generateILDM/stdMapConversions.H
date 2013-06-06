#ifndef STDMAPCONVERSIONS_H
#define STDMAPCONVERSIONS_H

#include <map>
#include <cantera/Cantera.h>
#include <cantera/IdealGasMix.h>
#include "octave/config.h"
#include "octave/dColVector.h"

ColumnVector toVector
(
    const Cantera::IdealGasMix& gas, 
    const Cantera::compositionMap& Y, 
    std::string leaveOut=""
);

Cantera::compositionMap toMap
(
    const Cantera::IdealGasMix& gas, 
    const ColumnVector& vY, 
    std::string leaveOut=""
);

void toArray
(
    const Cantera::IdealGasMix& gas, 
    const Cantera::compositionMap& Y, 
    double* y
);

Cantera::compositionMap toMap
(
    const Cantera::IdealGasMix& gas, 
    double* y
);

void write
(
    std::ofstream& f,
    const Cantera::IdealGasMix& gas,
    const Cantera::compositionMap& mY
);

#endif
