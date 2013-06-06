#include "LDM.H"

namespace Foam
{
/*
    defineTypeNameAndDebug(interpolatedLDM, 0);
*/
    LDM::FinalState::FinalState(const ChemicalSystem& sys, LDMconservedVars fcv, bool unburnt)
        : cv(fcv)
    {
       if (unburnt)
          Y=sys.unburntComposition(fcv.z(), fcv.p(), fcv.h(), true);
       else
          Y=sys.equilibriumComposition(fcv.z(), fcv.p(), fcv.h(), true);
    }

    LDM::LDM
    (
        const LDM& lowest,
        LDMconservedVars cv,
        const FinalState& finalState
    )
    : sys_(lowest.sys()),
      pv_(lowest.pv()),
      pvidx_(sys_.gas().speciesIndex(pv_)),
      cv_(cv),
      nSlow_(gas().nElements()+1),
      P_(nSlow_, sys_.gas().nSpecies(), 0.0),
      resolution_(lowest.resolution())

{
    calcP();
        scalar diffz=lowest.z()-z();
        scalar diffp=lowest.p()-p();
        scalar diffh=lowest.h()-h();
        if (cv.unburntTfixed()) diffh=0.0;
        label num=0;
        if (mag(diffz)>SMALL) num++;
        if (mag(diffp)>SMALL) num++;
        if (mag(diffh)>SMALL) num++;
        if (num!=1)
        {
            FatalErrorIn
                (
                    "interpolatedLDM::interpolatedLDM()"
                )
                << "Interpolation is possible only in one direction."<<endl
                    << "diffz="<<diffz<<" diffp="<<diffp<<" diffh="<<diffh
                      << exit(FatalError);
        }
        scalar diff=max(diffz, max(diffp, diffh));
        scalar diff2z=lowest.z()-finalState.cv.z();
        scalar diff2p=lowest.p()-finalState.cv.p();
        scalar diff2h=lowest.h()-finalState.cv.h();
        scalar difftot=max(diff2z, max(diff2p, diff2h));

        for(LDM::const_iterator it=lowest.begin();
            it!=lowest.end(); it++)
        {
            scalar weight1=(difftot-diff)/difftot;
            scalar weight2=diff/difftot;

            ColumnVector finalLDMPoint(it->second.length(), 0.0);
            finalLDMPoint.insert(finalState.Y, 0);

            scalar c = weight1*(it->first) + weight2*Yp(finalLDMPoint);

            (*this)[c] = weight1*(it->second) + weight2*finalLDMPoint;
        }
    }
/*

    interpolatedLDM::~interpolatedLDM()
    {
    }

    void interpolatedLDM::calculate()
    {
    }
*/
}
