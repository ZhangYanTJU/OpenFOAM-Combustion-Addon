#include "ignitionSite.H"

ignitionSite::ignitionSite(const fvMesh& mesh, Istream& is)
{
 dictionary dict(is);
 pvid_=readLabel(dict.lookup("progressVariableID"));
 location_=vector(dict.lookup("location"));
 diameter_=readScalar(dict.lookup("diameter"));
 duration_=readScalar(dict.lookup("duration"));
 if (dict.found("rampDuration"))
  rampDuration_=readScalar(dict.lookup("rampDuration"));
 else
  rampDuration_=duration_;
 starttime_=readScalar(dict.lookup("start"));
 patchval_=readScalar(dict.lookup("patchValue"));

 forAll(mesh.C(), I)
  if ( mag(mesh.C()[I]-location_) < diameter_) igncells_.append(I);
 
 igncells_.shrink();
 //Pout<<"found cells:"<<igncells_<<endl;
}

bool ignitionSite::patchProgressVariable(const Time& runTime, volScalarField& pvf) const
{
    if ( runTime.value()-starttime_ < duration_ )
    {
        scalar fac=min(runTime.value()-starttime_, rampDuration_) / rampDuration_;
        forAll(igncells_, I)
        {
            //Pout<<"cell #"<<igncells_[I]<<": "<<pvf[igncells_[I]]<<endl;
            pvf[igncells_[I]]=
                max(fac*patchval_, pvf[igncells_[I]]);
        }
        return true;
    } else return false;
}
