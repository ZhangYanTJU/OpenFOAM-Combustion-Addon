#include "integrationScheduleEntry.H"

namespace Foam
{


bool integrationScheduleEntry::operator!=(const integrationScheduleEntry& o) const
{
    FatalErrorIn("integrationScheduleEntry::operator!=()")
        <<"Not implemented."<<abort(FatalError);
    return true;
}


Ostream& operator<<(Ostream& os, const integrationScheduleEntry& e)
{
    os<<e.piCursor<<nl<<e.integratedState<<nl<<e.ix<<nl<<e.constants;
        /*
    FatalErrorIn("integrationScheduleEntry::operator<<()")
        <<"Not implemented."<<abort(FatalError);
            */
    return os;
}

Istream& operator>>(Istream& os, integrationScheduleEntry& e)
{
    os>>e.piCursor>>e.integratedState>>e.ix>>e.constants;
    /*
    FatalErrorIn("integrationScheduleEntry::operator>>()")
        <<"Not implemented."<<abort(FatalError);
        */
    return os;
}

}
