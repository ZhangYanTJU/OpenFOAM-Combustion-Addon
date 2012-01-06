#ifndef COLINWRINKLING_H
#define COLINWRINKLING_H

#include "efficiencyFunction.H"

namespace Foam
{

    class colinWrinkling
    : public efficiencyFunction
    {

        volScalarField uPrime_;
        scalar alpha_;
        dimensionedScalar deltal0_, sl0_;

    public:
        TypeName("colinWrinkling");

        colinWrinkling(const volVectorField& U, const scalar& TF);
        ~colinWrinkling();

        virtual void correct();
    };

}

#endif
