
#include "ILDM.H"
#include "REDIM.H"

#include "fvCFD.H"
#include "IOdictionary.H"
#include "OFstream.H"
#include "chemistryTable.H"
#include "generationLoop.H"
#include "mixtureFractionGenerationLoop.H"
#include <fstream>


using namespace Foam;


std::map<std::string, double> toStdMap(Foam::HashTable<double>& tab)
{
  std::map<std::string, double> copy;
  for (Foam::HashTable<double>::const_iterator it=tab.begin();
       it!=tab.end();it++)
    copy[it.key()]=it();
  return copy;
}





int main(int argc, char* argv[])
{
    try
    {
        // setRootCase.H causes errors in cantera!?
        Foam::Time runTime
            (
                Foam::Time::controlDictName,
                argv[1],
                argv[2]
            );

        IOdictionary generateIldmDict
            (
                IOobject
                (
                    "generateIldmDict",
                    runTime.system(),
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

        HashTable<scalar> f(generateIldmDict.lookup("fuel"));
        HashTable<scalar> o(generateIldmDict.lookup("oxidant"));
        HashTable<scalar> p(generateIldmDict.lookup("majorProducts"));
  
        ChemicalSystem chemsys
            (
                Foam::string(generateIldmDict.lookup("inputFile")),
                Foam::string(generateIldmDict.lookup("gasID")),
                toStdMap(f),
                toStdMap(o),
                toStdMap(p),
                word(generateIldmDict.lookup("constituentsGivenAs")),
                word(generateIldmDict.lookup("mixtureFractionBasedOn")),
                readScalar(generateIldmDict.lookup("stoichiometricOxidantFuelMassRatio"))
            );




        Info << " Creating loops"<<endl;

        // create loops over conserved variables
        PtrList<generationLoop> loops
            (
                generateIldmDict.lookup("conservedVariableLoops"),
                generationLoop::iNew(chemsys)
            );

        forAll(loops, J)
            loops[J].populate();



        Info<<"Creating chemistry table"<<endl;

        // **************** create chemistry table ********************

        autoPtr<LDM> ldm = LDM::New
            (
                chemsys, LDMconservedVars(chemsys, 0, 0, 0), generateIldmDict
            );

        ListDescription pvdesc, contentdesc;
        forAll(loops, J)
            loops[J].appendTableDescriptions(pvdesc, contentdesc);

        // progress variable after conserved variables!
        ldm().appendTableDescriptions(pvdesc, contentdesc);

        label dims=pvdesc.size();

        label num[dims];
        scalar start[dims];
        scalar delta[dims];

        forAll(loops, J)
            loops[J].setupTable(pvdesc, contentdesc, num, start, delta);
        ldm().setupTable(pvdesc, contentdesc, num, start, delta);

        Info<<"Allocating "<<dims<<"-dimensional chemistry table"<<endl;

        Info<<"PV's"<<pvdesc<<endl;
        Info<<"content"<<contentdesc<<endl;

        chemistryTable table
            (
                dims,
                num, start, delta,
                contentdesc, pvdesc
            );



        Info<<"Start generation of chemistry table"<<endl;

        LDMconservedVars defaultcv
            (
                chemsys, 
                generateIldmDict.lookup("defaultConservedVariables")
            );


        generationLoop::generateTable(chemsys, generateIldmDict, defaultcv, table);

  

        Info<<"Writing LDMtable"<<endl;

        IOdictionary tableDict
            (
                IOobject
                (
                    "LDMtable",
                    runTime.constant(),
                    runTime,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                )
            );

        OFstream tableDictFile(tableDict.objectPath());
        tableDict.writeHeader(tableDictFile);

        forAll(loops, J)
            loops[J].writeAdditionalInformationToTableFile(tableDictFile);

        // look for loop over z
        const mixtureFractionGenerationLoop* zloop = NULL;
        forAll(loops, I)
            if (loops[I].varName()=="z")
            { zloop=(const mixtureFractionGenerationLoop*) loops(I); break; }
        if (zloop)
        {

            scalar zmax=zloop->end();
            table.rescale(table.indexOfPV("z"), 1.0/zmax);

            tableDictFile
                <<"constantNormValue_z"
                    <<token::SPACE
                    <<zmax
                    <<token::END_STATEMENT
                    <<endl;

            List<scalar> c(table.nProgressVariables(), 1.0);

            graph normValues=table.slice
                ( c, "z", "Y"&ldm().pv() );

            tableDictFile
                <<"normValues_"
                    <<word(ldm().pv())
                    <<token::SPACE
                    <<normValues.title()
                    <<token::SPACE
                    <<normValues.xName()
                    <<token::SPACE
                    <<normValues.yName()
                    <<token::SPACE
                    <<"("
                    <<normValues
                    <<")"
                    <<token::END_STATEMENT
                    <<endl;
        }
        else
        {
            tableDictFile
                <<"constantNormValue_" << word(ldm().pv())
                    <<token::SPACE
                    <<ldm().Yp
                (
                    chemsys.equilibriumComposition
                    (
                        defaultcv.z(),
                        defaultcv.p(),
                        defaultcv.h(),
                        true
                    )
                )
                    <<token::END_STATEMENT
                    <<endl;
        }
      
        tableDictFile
            <<"tableData"
                <<token::SPACE
                <<table
                <<token::END_STATEMENT
                <<endl;


    }
    catch (CanteraError)
    {
        showErrors(std::cout);
    }
    return 0;
}

