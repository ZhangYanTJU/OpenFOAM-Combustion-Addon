#ifndef EFFICIENCYFUNCTION_H
#define EFFICIENCYFUNCTION_H

#include "fvCFD.H"
#include "LESdelta.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class efficiencyFunction
: public IOdictionary
{

protected:
    const volVectorField& U_;
    const scalar& TF_;

    volScalarField efficiency_;

    autoPtr<LESdelta> delta_;

public:
    TypeName("efficiencyFunction");

    declareRunTimeSelectionTable
    (
        autoPtr,
        efficiencyFunction,
        dictionary,
        (
            const volVectorField& U,
            const scalar& TF
        ),
        (U, TF)
    );

    //- Return a reference to the selected LES model
    static autoPtr<efficiencyFunction> New
    (
        const volVectorField& U,
        const scalar& TF
    );



    efficiencyFunction(const volVectorField& U, const scalar& TF);
    ~efficiencyFunction();

    inline const volScalarField& eff() const 
        { return efficiency_; }

    //- Access function to filter width
    inline const volScalarField& delta() const
        {
            return delta_();
        }    
 

    virtual void correct() =0;
};

}

#endif
