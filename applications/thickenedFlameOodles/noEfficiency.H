#ifndef NOEFFICIENCY_H
#define NOEFFICIENCY_H

#include "efficiencyFunction.H"

namespace Foam
{

    class noEfficiency
    : public efficiencyFunction
    {

    public:
        TypeName("none");

        noEfficiency(const volVectorField& U, const scalar& TF);
        ~noEfficiency();

        virtual void correct();
    };

}

#endif
