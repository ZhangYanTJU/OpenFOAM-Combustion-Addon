#ifndef PDF_H
#define PDF_H

#include "integrableFunction.H"

namespace Foam
{

class pdf
: public sumOfIntegrableFunctions<scalar>
{

  word name_;
  word parameterName_;

  label resFM_;
  label resSM_;

public:
  TypeName("pdf");

  // Declare run-time constructor selection table

  declareRunTimeSelectionTable
  (
      autoPtr,
      pdf,
      dictionary,
      (
          const dictionary& dict
      ),
      (dict)
  );

  // Selectors

  //- Return a reference to the selected LES model
  static autoPtr<pdf> New
  (
      const dictionary& dict
  );

  class iNew
  {

  public:

      iNew()

      {}

      autoPtr<pdf> operator()(Istream& is) const
      {
          return pdf::New(dictionary(is));
      }
  };

  pdf(const pdf&);
  pdf(const dictionary&);
  pdf(const word&, label res1=20, label res2=20);
  pdf(const word&, const word&, label res1=20, label res2=20);
  virtual ~pdf();
  
  inline const word& name() const
  { return name_; }
  
  inline const word& parameterName() const
  { return parameterName_; }

  inline const label& tabResolutionFM() const
  { return resFM_; }
  inline const label& tabResolutionSM() const
  { return resSM_; }

  virtual void setParameters(scalar, scalar) =0;

  virtual autoPtr<pdf> clone() const =0;

};

}

#endif
