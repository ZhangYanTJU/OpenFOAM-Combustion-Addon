#include "linearDampingFunction.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(linearDampingFunction, 0);
addToRunTimeSelectionTable(wallQuenchingDampingFunction, linearDampingFunction, dictionary);


linearDampingFunction::linearDampingFunction(const fvMesh& mesh, const dictionary& d)
    : wallQuenchingDampingFunction(mesh, d),
      dist_(d.lookup("distance")), // 5e-3
      fac_(d.lookup("scale")) // [1, -3, -1, 0, 0], 1e2
{
}


tmp<volScalarField> linearDampingFunction::operator()(const volScalarField& fm) const
{
    return -fac_ * fm * neg(wallY() - dist_);
}


}
