#ifndef CHEMISTRYTABLE_H
#define CHEMISTRYTABLE_H

#include "scalar.H"
#include "chemicalSystemState.H"
#include "MultidimensionalLookupTable.H"
#include "integrableFunction.H"
#include "graph.H"

namespace Foam
{

class chemistryTable;

Ostream& operator<< ( Ostream&, const chemistryTable& );


class chemistryTableEntryIntegrand
  : public integrableFunction<scalar>
{
  const chemistryTable& table_;
  label contentVariable_;

public:

  chemistryTableEntryIntegrand
  (
    const chemistryTableEntryIntegrand&
  );
  
  chemistryTableEntryIntegrand
  (
      const chemistryTable&,
      const word&
  );

  scalar valueAt ( const integrationParameter& ) const;

  integrableFunction<scalar>* clone() const;
};



class ListDescription
: public HashTable<label,word>
{
public:
    ListDescription();
    ListDescription(const ListDescription&);
    ListDescription(const HashTable<label,word>&);

    void appendEntry(const word&);
};



class chemistryTable
: public MultidimensionalLookupTable<chemicalSystemState>
{

  /**
  hashtable containing the names of the progress variables
  and their indices
  */
  ListDescription pvdesc_;

  /**
  hashtable containing the names of all fields in a
  single table entry. conventions:
  Y_*: mass fractions (if present)
  Ydot_*: reaction rates
  M: mixture molecular weight
  rho: mixture density
  */
  ListDescription contentdesc_;


public:

  chemistryTable
  (
      const label,
      const label*,
      const scalar*,
      const scalar*,
      const ListDescription&,
      const ListDescription&
  );

  chemistryTable ( Istream& );

  ~chemistryTable();

  chemistryTableEntryIntegrand entry(const word&) const;


  const ListDescription& contentDescription() const
  { return contentdesc_; }

  const ListDescription& progressVariableDescription() const
  { return pvdesc_; }

  inline label indexOfPV ( const word& n ) const
  {
    return pvdesc_[n];
  }

  inline label indexOf ( const word& n ) const
  {
    return contentdesc_[n];
  }

    inline label nProgressVariables() const
        {
            return pvdesc_.size();
        }

    inline label nContentElements() const
        {
            return contentdesc_.size();
        }

  word nameOf(label) const;

  word nameOfPV(label) const;

    chemicalSystemState reverseLookup
    (
        label pvidx, 
        const List<scalar>& cursor, 
        label contentidx, 
        scalar contentval,
        scalar& pvval,
        string& error
    ) const;

    scalar findPV
    (
        label pvidx, 
        const List<scalar>& cursor, 
        label contentidx, 
        scalar contentval
    ) const;

    graph slice(const List<scalar>& cursor, const word& along, const word& entry) const;

  friend Ostream& operator<< ( Ostream&, const chemistryTable& );

};


}

#endif
