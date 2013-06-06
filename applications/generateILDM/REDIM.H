#ifndef REDIM_H
#define REDIM_H

#include "ILDM.H"

class Spline;

namespace Foam
{

class REDIM
: public ILDM
{
    int nsp_;

    double constant_;
    double h_;
    int iterations_;

    std::auto_ptr<Spline> spline
    (
        const std::vector<ColumnVector>& Y
    ) const;

    auto_ptr<std::vector<ColumnVector> > splineInterpolate
    (
        const std::vector<ColumnVector>& Y
    ) const;

    void splineDerivatives
    (
        const std::vector<ColumnVector>& Y,
        std::vector<Matrix>& DPsiDtheta,
        std::vector<Matrix>& DPsiDtheta2
    ) const;


public:
    TypeName("REDIM");

    REDIM
    (
        const ChemicalSystem& s, 
        const std::string& pv,
        const std::string& solverName,
        double z, double p, double h,
        int resolution,
        double c
    );

    REDIM
    (
        const ChemicalSystem& s, 
        double z,
        double p,
        double h,
        dictionary& dict
    );


    ~REDIM();

    void calculate();

};

}

#endif
