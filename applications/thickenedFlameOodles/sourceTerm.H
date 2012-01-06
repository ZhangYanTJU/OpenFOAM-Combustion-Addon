#ifndef SOURCETERM_H
#define SOURCETERM_H

#include "fvCFD.H"
#include "hCombustionThermo.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class sourceTerm
: public IOdictionary
{

protected:

    dictionary coeffsDict_;
    const hCombustionThermo& thermo_;

    volScalarField omega_;

public:
    TypeName("sourceTerm");

    declareRunTimeSelectionTable
    (
        autoPtr,
        sourceTerm,
        dictionary,
        (
            //            const volScalarField& b
            const hCombustionThermo& thermo
        ),
        (thermo)
    );

    //- Return a reference to the selected LES model
    static autoPtr<sourceTerm> New
    (
        //const volScalarField& b
        const hCombustionThermo& thermo
    );



    sourceTerm
    (
        const word& type,
        //const volScalarField& b
        const hCombustionThermo& thermo
    );
    ~sourceTerm();

    inline const volScalarField& omega() const
        { return omega_; }

    // calculates omega=dY/dt from T and normalised fuel mass fraction b=(Yf-Yfb)/(Yfu-Yfb)
    virtual void correct(const volScalarField& T) =0;
    
};

}

#endif
