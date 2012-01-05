#ifndef CHEMICALSYSTEMSTATE_H
#define CHEMICALSYSTEMSTATE_H

#include "HashTable.H"
#include "Istream.H"
#include "scalarField.H"

namespace Foam
{

  class chemicalSystemState;

  chemicalSystemState operator*
  (
   const scalar&,
   const chemicalSystemState&
  );
  
  chemicalSystemState operator/
  (
   const chemicalSystemState&,
   const scalar&
  );
  
  chemicalSystemState operator+
  (
   const chemicalSystemState&,
   const chemicalSystemState&
  );
  
  class chemicalSystemState
    : public scalarField
  {
    /**
       marks if table entry is already calculated or not
    */
    bool hasValidContent_;
    
  public:
    chemicalSystemState();
    chemicalSystemState(const HashTable<label,word>&);
    chemicalSystemState(const HashTable<label,word>&, Istream&);
    ~chemicalSystemState();
    
    void operator=(const chemicalSystemState&);

    inline bool hasValidContent() const { return hasValidContent_; }
    
    inline void setValid() { hasValidContent_=true; }

    void operator+=(const chemicalSystemState&);
    
    friend chemicalSystemState operator*
    (
     const scalar&,
     const chemicalSystemState&
     );
  
    friend chemicalSystemState operator/
    (
     const chemicalSystemState&,
     const scalar&
     );
  
    friend chemicalSystemState operator+
    (
     const chemicalSystemState&,
     const chemicalSystemState&
    );
  
  };
  
}

#endif
