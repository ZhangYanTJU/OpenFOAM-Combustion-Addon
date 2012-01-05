#ifndef CLIPPEDGAUSSIAN_H
#define CLIPPEDGAUSSIAN_H

#include "scalar.H"
#include "pdf.H"
#include "chemistryTable.H"
#include "Ostream.H"

namespace Foam
{

    struct clippedGaussianParameterMapping
    {
      scalar mug_;
      scalar sigmag_;

      void operator+=(const clippedGaussianParameterMapping&);
      void operator/=(scalar);

    };

    clippedGaussianParameterMapping operator*
    (
     scalar,
     const clippedGaussianParameterMapping&
    );

    class clippedGaussianParameterMappingTable
      : public MultidimensionalLookupTable<clippedGaussianParameterMapping>
    {
    public:
      clippedGaussianParameterMappingTable();
      ~clippedGaussianParameterMappingTable();
    };

    class clippedGaussian
    : public pdf
    {
    public:
        TypeName("clippedGaussian");

        class Parameters
        : public integrableFunction<scalar>
        {
            word parameterName_;
            scalar mug_, sigmag_;
            scalar Al_, Ar_;
	    bool degenerated_;

            DynamicList< diracLocation1D > diracDeltas_;

        public:
  	    static clippedGaussianParameterMappingTable mapping_;

            Parameters(const Parameters&);
            Parameters(const word&, scalar, scalar);
            ~Parameters();

            inline scalar mu() const { return mug_; }
            inline scalar sigma() const { return sigmag_; }
            inline scalar Al() const { return Al_; }
            inline scalar Ar() const { return Ar_; }

            inline scalar interior
            (
                scalar x
            ) const
                {
                    if (degenerated_)
                     return 0.0;
                    else
                     return 
                        ( 0.3989422804014327/sigma() )*
                        exp(-0.5*pow( (x-mu())/sigma(),2));
                }
            
            scalar valueAt(const integrationParameter&) const;
        
            //label nDim() const;
        
            integrableFunction<scalar>* clone() const;

  	    friend class clippedGaussian;
        };


    public:

        clippedGaussian(const clippedGaussian&);
        clippedGaussian(const word&);
        clippedGaussian(const dictionary&);
        
        virtual ~clippedGaussian(); 
        
        void setParameters(const Parameters&);
        void setParameters(scalar, scalar);

        autoPtr<pdf> clone() const;


        static inline scalar Al(scalar mu, scalar sigma)
            {
                scalar z=-mu/sigma;
                return 1.0 + erf(z/1.414213562);
            }

        static inline scalar Ar(scalar mu, scalar sigma)
            {
                scalar z=(1.0-mu)/sigma;
                return 1.0 - erf(z/1.414213562);
            }

        /*
        static scalar integrate
        (
            const Parameters&,
            const integrableFunction&
        );
            */

    };

}

#endif
