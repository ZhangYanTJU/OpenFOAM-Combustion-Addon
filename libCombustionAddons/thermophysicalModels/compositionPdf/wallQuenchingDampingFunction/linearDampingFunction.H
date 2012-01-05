#ifndef LINEARDAMPINGFUNCTION_H
#define LINEARDAMPINGFUNCTION_H

#include "wallQuenchingDampingFunction.H"

namespace Foam
{

class linearDampingFunction
: public wallQuenchingDampingFunction
{

    dimensionedScalar dist_;
    dimensionedScalar fac_;

public:
    TypeName("linear");

    //- Construct from Istream
    linearDampingFunction(const fvMesh& mesh, const dictionary&);


    virtual tmp<volScalarField> operator()(const volScalarField& fm) const;
};

}

#endif
