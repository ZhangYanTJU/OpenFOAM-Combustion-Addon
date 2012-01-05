#include "pdf.H"

namespace Foam
{


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(pdf, 0);
defineRunTimeSelectionTable(pdf, dictionary);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<pdf> pdf::New
(
    const dictionary& dict
)
{
    word pdfTypeName;

    dict.lookup("PDFtype") >> pdfTypeName;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(pdfTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
            (
                "pdf::select()"
            )   << "Unknown pdf type " << pdfTypeName
                << endl << endl
                << "Valid pdf types are :" << endl
                << dictionaryConstructorTablePtr_->toc()
                << exit(FatalError);
    }

    Info << "Selecting pdf " << pdfTypeName << endl;
    return autoPtr<pdf>(cstrIter()(dict));

}

pdf::pdf(const pdf& o)
  : sumOfIntegrableFunctions<scalar>(o),
    name_(o.name_),
    parameterName_(o.parameterName_),
    resFM_(o.resFM_),
    resSM_(o.resSM_)                                       
{
}

pdf::pdf(const dictionary& dict)
  : parameterName_(dict.lookup("variable")),
    resFM_(readLabel(dict.lookup("preIntegrationResolutionFirstMoment"))),
    resSM_(readLabel(dict.lookup("preIntegrationResolutionSecondMoment")))
{
  if (dict.found("name"))
    name_=word(dict.lookup("name"));
  else
  {
      name_="";
  }
}

pdf::pdf(const word& pn, label res1, label res2)
: name_(pn),
  parameterName_(pn),
  resFM_(res1),
  resSM_(res2)
{
}

pdf::pdf(const word& n, const word& pn, label res1, label res2)
  : name_(n),
    parameterName_(pn),
    resFM_(res1),
    resSM_(res2)
{
}

pdf::~pdf()
{
}

}
