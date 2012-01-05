#ifndef RECTANGLE_H
#define RECTANGLE_H

#include "scalar.H"
//#include "compositionPdf.H"
#include "chemistryTable.H"

namespace Foam
{

    class rectangle
    {
    public:

        class Parameters
        {
            scalar a_, b_;
	  bool singleDirac_;

        public:
            Parameters(scalar, scalar);

            inline scalar a() const { return a_; }
            inline scalar b() const { return b_; }

	  friend class rectangle;
        };

    private:
        static const Parameters* Ipdf_;
        static const integrableFunction* I_;
        static double integrand(double x);

    public:

        static scalar integrate
        (
            const Parameters&,
            const integrableFunction&
        );

    };

}

#endif

