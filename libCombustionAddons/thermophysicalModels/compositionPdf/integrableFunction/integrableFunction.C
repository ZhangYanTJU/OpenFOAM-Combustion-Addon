#include "integrableFunction.H"
#include "interpolateXY.H"

namespace Foam
{
  
  integrationParameter::integrationParameter(const integrationParameter& p)
  : HashTable<scalar>(p)
  {
  }

    void integrationParameter::operator=(const integrationParameter& o)
    {
        HashTable<scalar>::operator=(o);
    }

  
  integrationParameter::integrationParameter(const parameterList& l, const scalar& def)
  : HashTable<scalar>()
  {
    for (parameterList::const_iterator it=l.begin();
         it!=l.end(); it++)
      insert(it.key(), def);
  }

  
  integrationParameter::integrationParameter
  (
    const integrationParameter& p1,
    const integrationParameter& p2
  )
  : HashTable<scalar>(p1)
  {
    for (const_iterator it=p2.begin();it!=p2.end();it++)
    {
      if (found(it.key()))
      {
        FatalErrorIn("integrationParameter(p1, p2)") <<endl
            << "Error in merging parameter sets: duplicate variable "<<it.key()<<endl
            << abort(FatalError);
      }
      else
      {
        insert(it.key(), it());
      }
    }
  }
  
  integrationParameter::integrationParameter(const word& name, const scalar& value)
  {
    insert(name, value);
  }

  integrationParameter::integrationParameter()
  {
  }
  
  void integrationParameter::setValuesFromArray(const scalar values[])
  {
    label i=0;
    for (iterator it=begin();it!=end();it++)
      it()=values[i++];
  }
  
  void integrationParameter::writeValuesToArray(scalar values[])
  {
    label i=0;
    for (iterator it=begin();it!=end();it++)
      values[i++]=it();
  }


  
  
  
  bool operator==(const integrationParameter& p1, const integrationParameter& p2)
  {
    bool equal=true;
    for (integrationParameter::const_iterator it=p1.begin();
         it!=p1.end(); it++)
      if (p2.found(it.key())) equal=equal && ( p2[it.key()] == it() );
    return equal;
  }
  
  
  parameterList::parameterList()
  {
  }
  
  void parameterList::operator=(const integrationParameter& x)
  {
    clear();
    label i=0;
    for (integrationParameter::const_iterator it=x.begin();
         it!=x.end(); it++)
      insert(it.key(), i++);
  }
  
  
  parameterList merge(const parameterList& l1, const parameterList& l2)
  {
    parameterList ret(l1);
    for (parameterList::const_iterator it=l2.begin();it!=l2.end();it++)
    {
      if (!l1.found(it.key()))
      {
        ret.insert(it.key(), 0);
      }
    }
    return ret;
  }
  
  
  
  
  template<>
  label toArray<scalarField>::size(const scalarField& f)
  {
    return f.size();
  }
  
  template<>
  void toArray<scalarField>::convert(const scalarField& value, scalar array[])
  {
    for (label i=0;i<value.size();i++)
      array[i]=value[i];
  }

  template<>
  void toArray<scalarField>::convertBack(const scalar array[], scalarField& value)
  {
    for (label i=0;i<value.size();i++)
      value[i]=array[i];
  }
  
  
  
  
  
  void diracDelta::updateLocation(const integrationParameter& x)
  {
    location_=x;
    parameters_=location_;
  }
  
  diracDelta::diracDelta(const diracDelta& d)
  : integrableFunction<scalar>(d),
    location_(d.location_),
    diracFactor_(d.diracFactor_)
  {
  }

    void diracDelta::operator=(const diracDelta& o)
    {
        integrableFunction<scalar>::operator=(o);
        location_=o.location_;
        diracFactor_=o.diracFactor_;
    }

  diracDelta::diracDelta(const word& pn, const diracLocation1D& loc)
  : location_(pn, loc.x_),
    diracFactor_(loc.factor_)
  {
  }

  
  diracDelta::diracDelta(const integrationParameter& x, const scalar& v)
  : location_(x),
    diracFactor_(v)
  {
    //Info<<"new Dirac at "<<location_<<endl<<"with factor="<<diracFactor_<<endl;
    updateLocation(x);
  }
  
  void diracDelta::operator*=(const diracDelta& od)
  {
    integrationParameter np(location_, od.location_);
    updateLocation(np);
    diracFactor_*=od.diracFactor_;
  }

    void diracDelta::modifyParameterList(parameterList& l) const
    {
        //Info<<"unmodified parameter list:"<<l<<endl;
        for (integrationParameter::const_iterator it=location_.begin();
             it!=location_.end(); it++)
            if (l.found(it.key())) l.erase(it.key());

        //Info<<"modified parameter list by diracDelta:"<<l<<endl;
    }
    
    void diracDelta::modifyParameter(integrationParameter& x) const
    {
        //Info<<"unmodified parameter:"<<x<<endl;
        for (integrationParameter::const_iterator it=location_.begin();
             it!=location_.end(); it++)
        {
            if (!x.found(it.key())) x.insert(it.key(), it());
            else x[it.key()]=0.5*(x[it.key()]+it());
        }
        //Info<<"modified parameter by diracDelta:"<<x<<endl;
    }
        
  
  scalar diracDelta::valueAt(const integrationParameter& x) const
  {
      /*if (x==location_) return diracFactor_;
      else return 0.0;
          */
      return diracFactor_;
  }
  
  label diracDelta::nDim() const
  {
    return 0;
  }
  
  integrableFunction<scalar>* diracDelta::clone() const
  {
    return new diracDelta(*this);
  }
  
  
  singleParameter::singleParameter(const singleParameter& p)
  : integrableFunction<scalar>(p),
    name_(p.name_)
  {
  }

    void singleParameter::operator=(const singleParameter& o)
    {
        integrableFunction<scalar>::operator=(o);
        name_=o.name_;
    }

  singleParameter::singleParameter(const word& name)
  : name_(name)
  {
    parameters_.insert(name_, 0);
  }


  scalar singleParameter::valueAt(const integrationParameter& x) const
  {
    return x[name_];
  }
  
  integrableFunction<scalar>* singleParameter::clone() const
  {
    return new singleParameter(name_);
  }



  
  integrableGraph::integrableGraph
  (
    const integrableGraph& o
  )
      : integrableFunction<scalar>(o),
        name_(o.name_),
        x_(o.x_),
        y_(o.y_)
  {
  }

    void integrableGraph::operator=(const integrableGraph& o)
    {
        integrableFunction<scalar>::operator=(o);
        name_=o.name_;
        x_=o.x_;
        y_=o.y_;
    }
  
  integrableGraph::integrableGraph
  (
    const word& n,
    const scalarField& x,
    const scalarField& y
  )
      : integrableFunction<scalar>(),
        name_(n),
        x_(x),
        y_(y)
  {
    parameters_.insert(name_, 0);
  }

  scalar integrableGraph::valueAt ( const integrationParameter& x ) const
  {
    return interpolateXY(x[name_], x_, y_);
  }

  integrableFunction<scalar>* integrableGraph::clone() const
  {
    return new integrableGraph(*this);
  }


  inverseIntegrableFunction::inverseIntegrableFunction
  (
    const integrableFunction<scalar>& f
  )
   : f1_(f)
  {
    this->parameters_=f.parameters();
  }

  scalar inverseIntegrableFunction::valueAt(const integrationParameter& x) const
  {
    return 1.0/f1_.valueAt(x);
  }
    
  integrableFunction<scalar>* inverseIntegrableFunction::clone() const
  {
    return new inverseIntegrableFunction(f1_);
  }

  inverseIntegrableFunction inv(const integrableFunction<scalar>& f)
  {
    return inverseIntegrableFunction(f);
  }
}
