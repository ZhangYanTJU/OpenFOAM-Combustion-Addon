#ifndef BETA_H
#define BETA_H

#include "scalar.H"
#include "pdf.H"
#include "chemistryTable.H"
#include "Ostream.H"


namespace Foam
{

class beta
  : public pdf
{

public:
    TypeName("beta");

    class Parameters
      : public integrableFunction<scalar>
    {

        word parameterName_;
        scalar alpha_, beta_, lG_;

        DynamicList< diracLocation1D > diracDeltas_;

        bool degenerated_;

        static const scalar epsilon/*=1e-3*/;
        static const scalar boundary/*=1e-3*/;

    public:
      	Parameters(const Parameters&);
        Parameters(const word&, scalar, scalar);
        ~Parameters();

        inline scalar alpha() const
        {
            return alpha_;
        }
        
        inline scalar beta() const
        {
            return beta_;
        }
        
        inline bool degenerated() const
        {
          return degenerated_;
        }

        inline scalar interior
        (
            scalar x
        ) const
        {
          if (!degenerated())
            return
              pow( x, alpha_-1. ) *
              pow( 1.-x, beta_-1. ) *
              lG_;
          else
            return 0.0;
        }

        scalar valueAt(const integrationParameter&) const;
        
        //label nDim() const;
        
        integrableFunction<scalar>* clone() const;
        
        friend class beta;
    };


public:

    beta(const beta&);
    beta(const word&);
    beta(const dictionary&);

    virtual ~beta(); 

    void setParameters(const Parameters&);
    void setParameters(scalar, scalar);

    autoPtr<pdf> clone() const;
};

}
#endif
