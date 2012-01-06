#ifndef AIRMIX_H
#define AIRMIX_H

#include "sourceTerm.H"

namespace Foam
{

class airmix
: public sourceTerm
{

    scalar A_;
    scalar TA_;
    scalar MF_;
    scalar nuF_;
    scalar nuO_;
    scalar phi_;
    scalar stOF_;

public:
    TypeName("airmix");

    airmix(const hCombustionThermo& thermo);
    ~airmix();

    virtual void correct(const volScalarField& T);
};

}

#endif
