#ifndef CHEMICALSYSTEM_H
#define CHEMICALSYSTEM_H

#include <cantera/Cantera.h>
#include <cantera/IdealGasMix.h>
#include <cantera/transport.h>

#include "octave/config.h"
#include "octave/dMatrix.h"
#include "octave/dColVector.h"
#include "octave/dbleSCHUR.h"
//#include "octave/NLEqn.h"

#include "compositionMapArithmetics.H"

#include "chemicalSystemState.H"

#include <memory>


Matrix nullspace(Matrix A);

class ChemicalSystem
{
    mutable Cantera::IdealGasMix gas_;
    mutable std::auto_ptr<Cantera::Transport> tr_;

    Cantera::compositionMap oxi_;
    Cantera::compositionMap fuel_;
    Cantera::compositionMap prod_;

    std::string mixFracAtom_;

    double stoichOF_;

    double cBr_, cProd_;

public:
    ChemicalSystem
    (
        const std::string& fname,
        const std::string& id,
        const Cantera::compositionMap& fuel,
        const Cantera::compositionMap& oxi,
        const Cantera::compositionMap& prod,
        const std::string& ctype,
        const std::string& mfatom,
        double stoichOF
    );

    inline double stoichOF() const { return stoichOF_; }
    inline double stoichZ() const { return 1.0/(1.0+stoichOF_); }

    inline Cantera::compositionMap stoichY() const
        {
            return 
                (1.0/(1.0+stoichOF_))*(fuel_)
                +
                (stoichOF_/(1.0+stoichOF_))*(oxi_);
        }

    inline Cantera::compositionMap Y(double z) const
        {
            return
                z*(fuel_)
                +
                (1.0-z)*(oxi_);
      }

    inline Cantera::IdealGasMix& gas() const { return gas_; }
    inline Cantera::Transport& tr() const { return *tr_; }

    double hmin(double z, double p) const;
    double hzero(double z, double p) const;

/*
    Foam::chemicalSystemState unburntState
    (
        const Foam::HashTable<Foam::label,Foam::word>& desc,
        double z, double p, double h
    ) const;

    Foam::chemicalSystemState equilibriumState
    (
        const Foam::HashTable<Foam::label,Foam::word>& desc,
        double z, double p, double h
    ) const;
*/
    ColumnVector equilibriumComposition( double z, double p, double h, bool convert=false ) const;
    ColumnVector unburntComposition( double z, double p, double h, bool convert=false ) const;

    double setGasToMixFracWithProducts(double z, double p, double Yprod, double temp) const;
    void setGasToTemp(double z, double p, double T) const;
    void setGasTo(double z, double p, double h, bool convert=false) const;
    void setGasTo(const double* Y, double p, double h) const;


    ColumnVector DYDt
    (
        const double* Y,
        double h,
        double p,
        double* rho_copy=NULL
    ) const;

    Matrix F
    (
        const double* Y,
        double h,
        double p,
        ColumnVector* DYDt_copy=NULL,
        double* rho_copy=NULL
    ) const;

    void subSpace
    (
        const Matrix& J,
        const ColumnVector& DYDt,
        int nSlow,
        Matrix* Zdach_fast=NULL,
        Matrix* Zslow=NULL,
        Matrix* Zdach_slow=NULL,
        Matrix* Zfast=NULL,
        Matrix* D=NULL
    ) const;

    ColumnVector calcSchur
    (
        const Matrix& J,
        const ColumnVector& DYDt,
        int nSlow,
        Matrix* Zdach_copy=NULL
    ) const;

};

#endif
