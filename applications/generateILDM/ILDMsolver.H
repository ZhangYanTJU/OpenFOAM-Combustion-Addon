#ifndef ILDMSOLVER_H
#define ILDMSOLVER_H

#include <string>
#include "word.H"
#include "autoPtr.H"
#include "HashTable.H"
#include "runTimeSelectionTables.H"

#include "ChemicalSystem.H"

namespace Foam
{

class ILDM;


class ILDMsolverError
{
    std::string errorMessage_;
public:
    ILDMsolverError(const std::string&);

    inline const std::string& errorMessage() const
        { return errorMessage_; }
};



class ILDMsolver
{

protected:
    const ILDM& ildm_;

public:
    TypeName("ILDMsolver");


    // Declare run-time constructor selection table
    
    declareRunTimeSelectionTable
    (
        autoPtr,
        ILDMsolver,
        dictionary,
        (
            ILDM& ildm
        ),
        (ildm)
    );

    static Foam::autoPtr<ILDMsolver> New
    (
        ILDM& ildm,
        const word& solverName
    );

    ILDMsolver(ILDM&);
    virtual ~ILDMsolver();

    virtual ColumnVector solve
    (
        double pv,
        const ColumnVector& inital
    ) =0;
};

}

#endif
