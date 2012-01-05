#include "chemicalSystemState.H"

namespace Foam
{

    chemicalSystemState::chemicalSystemState()
        : scalarField(0),
          hasValidContent_(false)
    {
    }

    chemicalSystemState::chemicalSystemState
    (
        const HashTable<label,word>& description
    )
        : scalarField(description.size(), 0.0),
          hasValidContent_(false)
    {
    }


    chemicalSystemState::chemicalSystemState
    (
        const HashTable<label,word>& description,
        Istream& f
    )
        : scalarField(f),
          hasValidContent_(readBool(f))
    {
    }

    chemicalSystemState::~chemicalSystemState()
    {
    }

    void chemicalSystemState::operator=(const chemicalSystemState& e)
    {
        scalarField::operator=(e);
        hasValidContent_=e.hasValidContent_;
    }

  void chemicalSystemState::operator+=(const chemicalSystemState& e)
  {
    scalarField::operator+=(e);
    hasValidContent_=hasValidContent_&&e.hasValidContent_;
  }

  chemicalSystemState operator*
  (
   const scalar& sc,
   const chemicalSystemState& st
   )
  {
    chemicalSystemState s(st);
    s.scalarField::operator*=(sc);
    return s;
  }
  
  chemicalSystemState operator/
  (
   const chemicalSystemState& st,
   const scalar& sc
  )
  {
    chemicalSystemState s(st);
    s.scalarField::operator/=(sc);
    return s;
  }
  
  chemicalSystemState operator+
  (
   const chemicalSystemState& st1,
   const chemicalSystemState& st2
  )
  {
    chemicalSystemState s(st1);
    s.scalarField::operator+=(st2);
    s.hasValidContent_=st1.hasValidContent_&&st2.hasValidContent_;
    return s;
  }

}
