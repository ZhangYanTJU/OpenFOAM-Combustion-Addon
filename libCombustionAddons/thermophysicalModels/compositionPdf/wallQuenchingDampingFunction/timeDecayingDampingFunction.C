#include "timeDecayingDampingFunction.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(timeDecayingDampingFunction, 0);
addToRunTimeSelectionTable(wallQuenchingDampingFunction, timeDecayingDampingFunction, dictionary);


timeDecayingDampingFunction::timeDecayingDampingFunction(const fvMesh& mesh, const dictionary& d)
    : wallQuenchingDampingFunction(mesh, d),
      dampFunc_(wallQuenchingDampingFunction::New(mesh, d, "dampingFunctionType")),
      endTime_(readScalar(d.lookup("endTime"))),
      rampDuration_(readScalar(d.lookup("rampDuration")))
{
}


tmp<volScalarField> timeDecayingDampingFunction::operator()(const volScalarField& fm) const
{
    scalar deltaT=endTime_ - cellDistFuncs::mesh().time().value();

    scalar fac=1.0;
    if ( (deltaT < rampDuration_) )
    {
        fac = deltaT / rampDuration_;
    }
    if (deltaT < 0) fac=0.0;

    return fac*dampFunc_()(fm);
}


}
