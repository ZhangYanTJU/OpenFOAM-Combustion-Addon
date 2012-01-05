#ifndef INTEGRATION_SCHEDULE_ENTRY_H
#define INTEGRATION_SCHEDULE_ENTRY_H

#include "chemicalSystemState.H"
#include "integrableFunction.H"

namespace Foam
{
    struct integrationScheduleEntry
    {

        List<int> piCursor;
        chemicalSystemState integratedState;
        integrationParameter ix;
        integrationParameter constants;

        bool operator!=(const integrationScheduleEntry& o) const;
    };

    Ostream& operator<<(Ostream& os, const integrationScheduleEntry& e);

    Istream& operator>>(Istream& os, integrationScheduleEntry& e);
}

#endif
