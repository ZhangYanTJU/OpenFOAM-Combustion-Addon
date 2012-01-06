#include "airmix.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

    defineTypeNameAndDebug(airmix, 0);
    addToRunTimeSelectionTable(sourceTerm, airmix, dictionary);
 

    airmix::airmix(/*const volScalarField& b*/ const hCombustionThermo& thermo)
        : sourceTerm(typeName, thermo),
          A_(readScalar(coeffsDict_.lookup("A"))),
          TA_(readScalar(coeffsDict_.lookup("TA"))),
          MF_(readScalar(coeffsDict_.lookup("MF"))),
          nuF_(readScalar(coeffsDict_.lookup("nuF"))),
          nuO_(readScalar(coeffsDict_.lookup("nuO"))),
          phi_(0.0),
          stOF_(0.0)
    {
        dimensionedScalar stof(thermo.lookup("stoichiometricAirFuelMassRatio"));
        stOF_=stof.value();
        if (!thermo_.composition().contains("ft"))
        {
            phi_=readScalar(coeffsDict_.lookup("phi"));
        }
    }

    airmix::~airmix()
    {
    }

    void airmix::correct(const volScalarField& T)
    {

        const scalar MO2=32;

        const volScalarField& b_ = thermo_.composition().Y("b");
        const volScalarField& rho = //thermo_.rho(); //thermo.rho has uncorrected BC's! Do not use
           T.db().lookupObject<volScalarField>("rho"); //lookup returns rho field from top level solver

        if (thermo_.composition().contains("ft"))
        {
            const volScalarField& ft=thermo_.composition().Y("ft");

            forAll(omega_, I)
            {
                scalar maxYF= ft[I];
                scalar YF= b_[I]*ft[I]
                    +(1.0 - b_[I])*max(thermo_.composition().fres(ft[I], stOF_), 0.0);
                scalar YO2= 0.233005 * (1.0 - ft[I] - (ft[I] - YF)*stOF_);

                omega_[I]=maxYF>SMALL ? 1e3* // from cgs
                    A_ * nuF_ * MF_
                    *pow( 1e-3*rho[I]*YF / MF_, nuF_ )  // rho is kg/m^3, change to cgs
                    *pow( 1e-3*rho[I]*YO2 / MO2, nuO_ )
                    *exp(-TA_/T[I])
                    /maxYF : 0.0;
            }

            forAll(omega_.boundaryField(), bI)
                forAll(omega_.boundaryField()[bI], fI)
            {
                scalar maxYF= ft.boundaryField()[bI][fI];
                scalar YF= b_.boundaryField()[bI][fI]*ft.boundaryField()[bI][fI]
                    +(1.0 - b_.boundaryField()[bI][fI])*
                    max(thermo_.composition().fres(ft.boundaryField()[bI][fI], stOF_), 0.0);
                scalar YO2= 0.233005 * (1.0 - ft.boundaryField()[bI][fI] 
                           - (ft.boundaryField()[bI][fI] - YF)*stOF_);

                omega_.boundaryField()[bI][fI]=maxYF > SMALL ? 1e3*
                    A_ * nuF_ * MF_
                    *pow( 1e-3*rho.boundaryField()[bI][fI]*YF / MF_, nuF_ )
                    *pow( 1e-3*rho.boundaryField()[bI][fI]*YO2 / MO2, nuO_ )
                    *exp(-TA_/T.boundaryField()[bI][fI])
                    /maxYF : 0.0;
            }
        }
        else
        {

            scalar maxYF=1.0/((stOF_/phi_)+1.0);
            scalar YLex=1.0 - maxYF - stOF_*maxYF;

            forAll(omega_, I)
            {
                scalar YF  =                   maxYF  * b_[I];
                scalar YO2 = 0.233005 * (1.0 - maxYF) * b_[I]
                           + 0.233005 * YLex * (1.0 - b_[I]);
                omega_[I]=1e3* // from cgs
                    A_ * nuF_ * MF_
                    *pow( 1e-3*rho[I]*YF / MF_, nuF_ )  // rho is kg/m^3, change to cgs
                    *pow( 1e-3*rho[I]*YO2 / MO2, nuO_ )
                    *exp(-TA_/T[I])
                    /maxYF;
            }

            forAll(omega_.boundaryField(), bI)
                forAll(omega_.boundaryField()[bI], fI)
            {
                scalar YF  =                   maxYF  * b_.boundaryField()[bI][fI];
                scalar YO2 = 0.233005 * (1.0 - maxYF) * b_.boundaryField()[bI][fI]
                           + 0.233005 * YLex * (1.0 - b_.boundaryField()[bI][fI]);
                omega_.boundaryField()[bI][fI]=1e3*
                    A_ * nuF_ * MF_
                    *pow( 1e-3*rho.boundaryField()[bI][fI]*YF / MF_, nuF_ )
                    *pow( 1e-3*rho.boundaryField()[bI][fI]*YO2 / MO2, nuO_ )
                    *exp(-TA_/T.boundaryField()[bI][fI])
                    /maxYF;
            }
        }
    }

}
