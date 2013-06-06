#ifndef MIXTUREFRACTION_H
#define MIXTUREFRACTION_H

#include <cantera/Cantera.h>
#include <cantera/IdealGasMix.h>

double Z
(
    int i, 
    Cantera::IdealGasMix& gas, 
    const Cantera::compositionMap& Y
);

#endif
