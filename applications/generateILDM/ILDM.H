#ifndef ILDM_H
#define ILDM_H

#include "LDM.H"

namespace Foam
{

class ILDM
: public LDM
{

    ColumnVector tau_;

    //Cantera::compositionMap predict
    ColumnVector predict
    (
        //const Cantera::compositionMap& Ystart,
        const ColumnVector& Ystart,
        double step
    );


    ColumnVector Yunburnt, Ybound;
    Matrix Porth;
    Matrix bdryT;

    word solverName_;

public:
    TypeName("ILDM");

    ILDM
    (
        const ChemicalSystem& s, 
        const std::string& pv,
        const std::string& solverName,
        double z, double p, double h,
        int resolution=20
    );

    ILDM
    (
        const ChemicalSystem& s, 
        double z,
        double p,
        double h,
        dictionary& dict
    );

    void calculate();

    void outputEigenvalueSpectrum(std::ofstream&);

    inline const ColumnVector& tau() const { return tau_; }

    ColumnVector integrand
    (
        const ColumnVector& Y,
        double
    ) const;

    ColumnVector objective(const ColumnVector& Y) const;

    Matrix jacobian(const ColumnVector& Y) const;

    bool onILDM
    (
        const ColumnVector&Y
    );
    /*
    virtual void addToTable
    (
        chemistryTable& tab,
        const LDMconservedVars& cv
    );

    virtual void appendTableDescriptions
    (
        ListDescription& pvdesc,
        ListDescription& contentdesc
    );

    virtual void setupTable
    (
        ListDescription& pvdesc,
        ListDescription& contentdesc,
        label* num,
        double* start,
        double* delta
    );

    virtual chemicalSystemState extendToTableEntry
    (
        const chemistryTable& tab,
        const ColumnVector& Y
    );
        */

};

}

#endif
