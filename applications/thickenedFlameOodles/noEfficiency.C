#include "noEfficiency.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(noEfficiency, 0);
addToRunTimeSelectionTable(efficiencyFunction, noEfficiency, dictionary);


noEfficiency::noEfficiency
(
    const volVectorField& U,
    const scalar& TF
)
    : efficiencyFunction(U, TF)
{
}

noEfficiency::~noEfficiency()
{
}

void noEfficiency::correct()
{
    // leave default (1)
}

}
