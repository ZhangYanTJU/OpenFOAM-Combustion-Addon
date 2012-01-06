#ifndef IGNITIONSITE_H
#define IGNITIONSITE_H

#include "fvCFD.H"

class ignitionSite
{

    vector location_;
    scalar diameter_;
    scalar duration_;
    scalar rampDuration_;
    scalar starttime_;

    label pvid_;
    scalar patchval_;

    DynamicList<label> igncells_;

public:

    class iNew
    {
        const fvMesh& mesh_;

    public:

        iNew(const fvMesh& mesh)
            :
            mesh_(mesh)
            {}

        autoPtr<ignitionSite> operator()(Istream& is) const
            {
                return autoPtr<ignitionSite>(new ignitionSite(mesh_, is));
            }
    };

    inline const vector& location() const { return location_; }

    ignitionSite(const fvMesh&, Istream&);

    inline label progressVariableID() const { return pvid_; };

    bool patchProgressVariable(const Time&, volScalarField&) const;

    autoPtr<ignitionSite> clone() const
        {
            return autoPtr<ignitionSite>(new ignitionSite(*this));
        }
 
};

#endif
