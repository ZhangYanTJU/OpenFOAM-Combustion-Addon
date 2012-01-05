#ifndef FLAMESURFACEDENSITY_H
#define FLAMESURFACEDENSITY_H

#include "scalar.H"
#include "beta.H"
#include "chemistryTable.H"
#include "Ostream.H"
#include "Istream.H"
#include "interpolationTable.H"

namespace Foam
{

class flameSurfaceDensity
: public beta
{

public:
    TypeName("flameSurfaceDensity");

    class Gdelta
    {

        scalar delta_;
        interpolationTable<scalar> Gdelta_;

        scalar a_, b_, d_;

    public:
        Gdelta(const Gdelta&);
        Gdelta(const dictionary&);

        inline scalar a() const { return a_; }
        inline scalar b() const { return b_; }
        inline scalar d() const { return d_; }

        scalar operator()(scalar) const;
    };

    class Parameters
      : public integrableFunction<scalar>
    {
        word parameterName_;
        scalar alpha_, beta_, sigma_;

        DynamicList< diracLocation1D > diracDeltas_;

        bool degenerated_;

        const Gdelta& gdelta_;
        
    public:

        Parameters(const Parameters& p);

        Parameters
        (
            const word& pv,
            scalar, scalar,
            const Gdelta&
        );

        inline bool degenerated() const
        {
          return degenerated_;
        }


        inline scalar interior
        (
            scalar x
        ) const
            {
                return sigma_/stabilise(gdelta_(x), SMALL);
            }

        scalar valueAt(const integrationParameter&) const;

        integrableFunction<scalar>* clone() const;
        
        friend class flameSurfaceDensity;

    };

    autoPtr<Gdelta> Gdelta_;

public:
    flameSurfaceDensity(const flameSurfaceDensity&);
    flameSurfaceDensity(const word&);
    flameSurfaceDensity(const dictionary&);

    virtual ~flameSurfaceDensity(); 

    //void setParameters(const Parameters&);
    void setParameters(scalar, scalar);

    autoPtr<pdf> clone() const;


};

}


#endif
