#ifndef TIMEDECAYINGDAMPINGFUNCTION_H
#define TIMEDECAYINGDAMPINGFUNCTION_H

#include "wallQuenchingDampingFunction.H"

namespace Foam
{

class timeDecayingDampingFunction
: public wallQuenchingDampingFunction
{

    autoPtr<wallQuenchingDampingFunction> dampFunc_;

    scalar endTime_;
    scalar rampDuration_;

public:
    TypeName("timeDecaying");

    //- Construct from Istream
    timeDecayingDampingFunction(const fvMesh& mesh, const dictionary&);


    virtual tmp<volScalarField> operator()(const volScalarField& fm) const;
};

}

#endif
