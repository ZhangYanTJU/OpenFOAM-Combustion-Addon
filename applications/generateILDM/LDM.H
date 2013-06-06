#ifndef LDM_H
#define LDM_H

#include "ChemicalSystem.H"
#include "stdMapConversions.H"
#include "typeInfo.H"
#include "Spline.H"
#include "chemistryTable.H"

namespace Foam
{

class LDMconservedVars
{

    const ChemicalSystem& sys_;
    double z_, p_, h_;

    bool unburntTfixed_;
    double unburntT_;

public:
    LDMconservedVars
    (
        const ChemicalSystem& sys,
        const HashTable<double, word>& sym
    );
    LDMconservedVars
    (
        const ChemicalSystem& sys,
        double z, double p, double h
    );

    inline bool unburntTfixed() const { return unburntTfixed_; }

    inline double z() const { return z_; }
    inline double p() const { return p_; }    
    double h() const;

    inline void setz(double z) { z_=z; }
    inline void setp(double p) { p_=p; }    
    void seth(double h);

    double getByName(const word& name) const;
};

class LDM;

std::ofstream& operator<<(std::ofstream&, const LDM&);
std::ifstream& operator>>(std::ifstream&, LDM&);



class LDM
: public std::map<double, ColumnVector>
{

    mutable autoPtr<double> pveq_; // computed on request

protected:
    const ChemicalSystem& sys_;

    std::string pv_;
    int pvidx_;

    //double z_, p_, h_;
    LDMconservedVars cv_;
    int nSlow_;
    Matrix P_;
    label resolution_;

    void calcP();

    void append(double pvval, const Cantera::compositionMap& compo);
    void append(double pvval, const ColumnVector& compo);

    string uniqueDescription() const;
    
public:
    
    //- Runtime type information
    TypeName("LDM");
    
    struct FinalState
    {
        LDMconservedVars cv;
        ColumnVector Y;

        FinalState(const ChemicalSystem& sys, LDMconservedVars fcv, bool unburnt=false);
    };

    // Declare run-time constructor selection table
    
    declareRunTimeSelectionTable
    (
        autoPtr,
        LDM,
        dictionary,
        (
            const ChemicalSystem& s, 
            LDMconservedVars cv,
            dictionary& dict
        ),
        (s, cv, dict)
    );

    declareRunTimeSelectionTable
    (
        autoPtr,
        LDM,
        interpolate,
        (
         const LDM& lowest,
         LDMconservedVars cv,
         const FinalState& finalState
        ),
        (lowest, cv, finalState)
    );


    //- Return a reference to the selected surfaceWriter
    static autoPtr<LDM> New
    (
        const ChemicalSystem& s, 
        LDMconservedVars cv,
        dictionary& dict
    );

    static autoPtr<LDM> New
    (
        const LDM& lowest,
        LDMconservedVars cv,
        const FinalState& finalState,
        dictionary& dict
    );

    LDM
    (
        const LDM& lowest,
        LDMconservedVars cv,
        const FinalState& finalState
    );

    LDM(const LDM&);

    LDM
    (
        const ChemicalSystem& s, 
        const std::string& pvname,
        LDMconservedVars cv,
        label resolution=20
    );

    LDM
    (
        const ChemicalSystem& s, 
        LDMconservedVars cv,
        dictionary& dict
    );

    virtual ~LDM();

    double pveq() const;

    virtual void calculate() =0;

    inline Cantera::IdealGasMix& gas() const { return sys_.gas(); }
    inline const ChemicalSystem& sys() const { return sys_; }

    inline double Yp(const ColumnVector& Y) const
        {
            return Y(pvidx_);
        }

    inline const std::string& pv() const { return pv_; }
    inline double z() const { return cv_.z(); };
    inline double p() const { return cv_.p(); };
    inline double h() const { return cv_.h(); };
    inline const int& nSlow() const { return nSlow_; };
    inline const label& resolution() const { return resolution_; };

    inline const Matrix& P() const { return P_; }

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

    friend
    std::ofstream& operator<<(std::ofstream&, const LDM&);

    friend
    std::ifstream& operator>>(std::ifstream&, LDM&);

};

}

#endif
