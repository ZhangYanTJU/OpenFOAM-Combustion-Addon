#ifndef INTEGRABLEFUNCTION_H
#define INTEGRABLEFUNCTION_H

#include "scalar.H"
#include "scalarField.H"
#include "Map.H"
#include "refCount.H"
#include "tmp.H"

#include <vector>

namespace Foam
{
  
class parameterList;
  
class integrationParameter
 : public HashTable<scalar>
{
public:
  
  /**
  copy constructor
  */
  integrationParameter(const integrationParameter&);
  
  /**
  initialize with names from parameterList and default value
  */
  integrationParameter(const parameterList&, const scalar&);
  
  /**
  merge two parameter sets
  */
  integrationParameter
      (
       const integrationParameter&,
       const integrationParameter&
      );
  
  /**
  create single parameter
  */
  integrationParameter(const word&, const scalar&);
  
  /**
  empty parameter set
  */
  integrationParameter();
  
  void operator=(const integrationParameter&);
  
  /**
  set values from array
  attention: ordering of entries must match
  */
  void setValuesFromArray(const scalar[]);
  
  /**
  write values to array
  */
  void writeValuesToArray(scalar[]);
  
};

class parameterList
  : public HashTable<label,word> 
{
  public:
    parameterList();
    
    void operator=(const integrationParameter&);
};

bool operator==(const integrationParameter&, const integrationParameter&);

parameterList merge(const parameterList&, const parameterList&);

/**
generic class for converting arbitrary types to arrays
(needed for component-by-component integration)
*/
template<class T>
class toArray
{
public:
  static label size(const T&)
  {
    return pTraits<T>::nComponents;
  }
  
  static void convert(const T& value, scalar array[])
  {
    for (label i=0;i < pTraits<T>::nComponents;i++)
      array[i]=component(value, i);
  }
  
  static void convertBack(const scalar array[], T& value)
  {
    for (label i=0;i < pTraits<T>::nComponents;i++)
     setComponent(value, i)=array[i];
  }
};

template<>
label toArray<scalarField>::size(const scalarField&);

template<>
void toArray<scalarField>::convert(const scalarField& value, scalar array[]);

template<>
void toArray<scalarField>::convertBack(const scalar array[], scalarField& value);



/**
 A integrand function of arbitrary dimensionality
*/
template<class T>
class integrableFunction
  : public refCount
{
  
  protected:
    parameterList parameters_;

  public:
    
    integrableFunction()
      : refCount()
    {}
    
    integrableFunction(const integrableFunction& f)
      : refCount(),
        parameters_(f.parameters_)
    {}
    
    virtual ~integrableFunction()
    {}

    void operator=(const integrableFunction& f)
        {
            //refCount::operator=(f);
            parameters_=f.parameters_;
        }
    
    inline const parameterList& parameters() const
    { return parameters_; }

    /**
        modify parameterList of a product of this function and another
        (useful in case this represents a dirac delta)
        */
    virtual void modifyParameterList(parameterList&) const
        {
        }

    /**
        modify parameter of a product of this function and another
        (useful in case this represents a dirac delta)
        */
    virtual void modifyParameter(integrationParameter&) const
        {
        }

    virtual T valueAt(const integrationParameter&) const =0;

    virtual label nDim() const
    {
      return parameters_.size();
    }
    
    virtual integrableFunction<T>* clone() const =0;

};

struct diracLocation1D
{
  scalar x_;
  scalar factor_;

  diracLocation1D()
  : x_(0.0), factor_(1.0)
  {}
  
  diracLocation1D(scalar x, scalar f)
  : x_(x), factor_(f)
  {}
};



class diracDelta
  : public integrableFunction<scalar>
{
  integrationParameter location_;
  scalar diracFactor_;
  
  void updateLocation(const integrationParameter&);

public:
  diracDelta(const diracDelta&);
  diracDelta(const word&, const diracLocation1D&);
  diracDelta(const integrationParameter&, const scalar&);

    void operator=(const diracDelta&);
  
  inline const scalar& factor() const { return diracFactor_; }
  
  inline const integrationParameter& x() const { return location_; }
  inline integrationParameter& x() { return location_; }
  
    virtual void modifyParameterList(parameterList&) const;
    virtual void modifyParameter(integrationParameter&) const;

  scalar valueAt(const integrationParameter&) const;
  
  void operator*=(const diracDelta&);
  
  label nDim() const;
  
  integrableFunction<scalar>* clone() const;
};



class singleParameter
 : public integrableFunction<scalar>
{
  word name_;
  
public:
  singleParameter(const singleParameter&);
  singleParameter(const word&);

    void operator=(const singleParameter&);
    
  scalar valueAt(const integrationParameter&) const;

  integrableFunction<scalar>* clone() const;

};


template<class T>
class Constant
 : public integrableFunction<T>
{
    T value_;
public:
    Constant(const Constant& o)
        : integrableFunction<T>(o),
          value_(o.value_)
        {
        }

    Constant(const T& v)
        : value_(v)
        {
        }

    void operator=(const Constant& o)
        {
            integrableFunction<T>::operator=(o);
            value_=o.value_;
        }
    
    T valueAt(const integrationParameter& x) const
        {
            return value_;
        }

    integrableFunction<T>* clone() const
        {
            return new Constant<T>(*this);
        }

};




class integrableGraph
  : public integrableFunction<scalar>
{
  word name_;
  scalarField x_;
  scalarField y_;

  
  public:

    integrableGraph
        (
         const integrableGraph&
        );
  
    integrableGraph
        (
         const word&,
         const scalarField&,
         const scalarField&
        );

    void operator=(const integrableGraph&);

    scalar valueAt ( const integrationParameter& ) const;

    integrableFunction<scalar>* clone() const;
};




template<class T, class T2>
class multipliedIntegrableFunction
 : public integrableFunction<typename outerProduct <T, T2>::type>
{
protected:
    const integrableFunction<T>* f1_;
    const integrableFunction<T2>* f2_;
    
    void createParameterList()
    {

      this->parameters_ = merge(f1_->parameters(), f2_->parameters());

      modifyParameterList(this->parameters_);

      /*
      if (typeid(f1_)==typeid(diracDelta))
      {
        const diracDelta &d=dynamic_cast<const diracDelta&>(f1_);
        for (integrationParameter::const_iterator it=d.x().begin();
             it!=d.x().end(); it++)
          if (this->parameters_.found(it.key())) this->parameters_.erase(it.key());
      }

      if (typeid(f2_)==typeid(diracDelta))
      {
        const diracDelta &d=dynamic_cast<const diracDelta&>(f2_);
        for (integrationParameter::const_iterator it=d.x().begin();
             it!=d.x().end(); it++)
          if (this->parameters_.found(it.key())) this->parameters_.erase(it.key());
      }
          */
    }
    
    /*
    integrationParameter filterParameters(const integrationParameter& x, const diracDelta& dirac) const
    {
      integrationParameter xnew;

      // Reduce dimensionality of integrand function
      for (integrationParameter::const_iterator it=dirac.x().begin();
           it!=dirac.x().end(); it++)
        xnew.insert(it.key(), it());
      
      for (integrationParameter::const_iterator it=x.begin();
           it!=x.end(); it++)
        if (!dirac.x().found(it.key())) xnew.insert(it.key(), it());

      return xnew;
    }
        */
    
public:

    multipliedIntegrableFunction()
        {
        }

    multipliedIntegrableFunction
    (
      const multipliedIntegrableFunction& m
    )
      : integrableFunction<typename outerProduct <T, T2>::type>(m),
        f1_(m.f1_),
        f2_(m.f2_)
    {
    }
    
    multipliedIntegrableFunction
    (
        const integrableFunction<T>& f1,
        const integrableFunction<T2>& f2
    )
    : f1_(&f1),
      f2_(&f2)
    {
        setFunctions(f1, f2);
    }

    void setFunctions
    (
        const integrableFunction<T>& f1,
        const integrableFunction<T2>& f2
    )
        {
            f1_=&f1;
            f2_=&f2;
            createParameterList();
        }

    void operator=(const multipliedIntegrableFunction& o)
        {
            integrableFunction<typename outerProduct <T, T2>::type>::operator=(o);
            f1_=o.f1_;
            f2_=o.f2_;
        }

    virtual void modifyParameterList(parameterList& l) const
        {
            f1_->modifyParameterList(l);
            f2_->modifyParameterList(l);
        }

    virtual void modifyParameter(integrationParameter& x) const
        {
            f1_->modifyParameter(x);
            f2_->modifyParameter(x);
        }

    typename outerProduct<T, T2>::type valueAt(const integrationParameter& x) const
    {

        integrationParameter X(x);
        //Info<<"X bef:"<<X<<endl;
        modifyParameter(X);
        //Info<<"X after:"<<X<<endl;

        return f1_->valueAt(X) * f2_->valueAt(X);

        /*
      if ((typeid(f1_) == typeid(diracDelta)) && (typeid(f2_) != typeid(diracDelta)))
      {
        //Info<<"f1>d"<<endl;
        const diracDelta &d=dynamic_cast<const diracDelta&>(f1_);
        return d.factor() * f2_.valueAt(filterParameters(x, d));
      }
      else if ((typeid(f1_) != typeid(diracDelta)) && (typeid(f2_) == typeid(diracDelta)))
      {
        //Info<<"f2>d"<<endl;
        const diracDelta &d=dynamic_cast<const diracDelta&>(f2_);
        return d.factor() * f1_.valueAt(filterParameters(x, d));
      }
      else if ((typeid(f1_) == typeid(diracDelta)) && (typeid(f2_) == typeid(diracDelta)))
      {
        //Info<<"f1,f2>d"<<endl;
        diracDelta d(dynamic_cast<const diracDelta&>(f1_));
        d*=dynamic_cast<const diracDelta&>(f2_);
        return d.valueAt(x);
      }
      else
      {
        //Info<<"no d"<<endl;
        //Info<<this->parameters_<<endl;
        return f1_.valueAt(x) * f2_.valueAt(x);
      }
            */
    }
    
    integrableFunction<typename outerProduct <T, T2>::type>* clone() const
    {
      return new multipliedIntegrableFunction<T, T2>(*this);
    }

};


class inverseIntegrableFunction
  : public integrableFunction<scalar>
{
  const integrableFunction<scalar>& f1_;

  public:
    inverseIntegrableFunction
    (
      const integrableFunction<scalar>& f
    );

    scalar valueAt(const integrationParameter& x) const;

    integrableFunction<scalar>* clone() const;

};


template<class T>
class dividedByConstantIntegrableFunction
: public multipliedIntegrableFunction<scalar, T>
{
    Constant<scalar> c_;
    inverseIntegrableFunction cinv_;

public:
    dividedByConstantIntegrableFunction
    (
        const dividedByConstantIntegrableFunction& o
    )
        : multipliedIntegrableFunction<scalar, T>(),
          c_(o.c_),
          cinv_(c_)
    {
        setFunctions(cinv_, *o.f2_);
    }

    dividedByConstantIntegrableFunction
    (
        const integrableFunction<T>& f,
        const scalar& c
    )
    : c_(c),
      cinv_(c_)
    {
        setFunctions(cinv_, f);
    }

    void operator=(const dividedByConstantIntegrableFunction& o)
        {
            c_=o.c_;
            cinv_=inverseIntegrableFunction(c_);
            setFunctions(cinv_, *o.f2_);
        }

    integrableFunction<T>* clone() const
    {
      return new dividedByConstantIntegrableFunction<T>(*this);
    }

};

template<class T>
class multipliedByConstantIntegrableFunction
: public multipliedIntegrableFunction<scalar, T>
{

    Constant<scalar> c_;

public:
    multipliedByConstantIntegrableFunction(const multipliedByConstantIntegrableFunction& o)
        : multipliedIntegrableFunction<scalar, T>(),
          c_(o.c_)
        {
            setFunctions(c_, *o.f2_);
        }

    multipliedByConstantIntegrableFunction
    (
        const scalar& c,
        const integrableFunction<T>& f
    )
        : c_(c)
    {
        setFunctions(c_, f);
    }
    
    void operator=(const multipliedByConstantIntegrableFunction& o)
        {
            c_=o.c_;
            setFunctions(c_, *o.f2_);
        }

    integrableFunction<T>* clone() const
    {
        return new multipliedByConstantIntegrableFunction<T>(*this);
    }

};



template<class T>
class sumOfIntegrableFunctions
 : public std::vector< tmp< integrableFunction<T> > >
{

  public:
    sumOfIntegrableFunctions()
    {
    }
    
    sumOfIntegrableFunctions(const sumOfIntegrableFunctions& o)
      : std::vector< tmp< integrableFunction<T> > >(o)
    {
    }
    
    sumOfIntegrableFunctions(const integrableFunction<T>& f)
    {
      operator+=(f);
    }

    sumOfIntegrableFunctions
    (
      const integrableFunction<T>& f1,
      const integrableFunction<T>& f2
    )
    {
      operator+=(f1);
      operator+=(f2);
    }

    void operator=(const sumOfIntegrableFunctions& o)
        {
            std::vector< tmp< integrableFunction<T> > >::operator=(o);
        }

    void operator+=(const integrableFunction<T>& f)
    {
      push_back(tmp<integrableFunction<T> >( f.clone() ));
    }

};



/* ------------------- sum of integrable functions ---------------------- */

template<class T>
const sumOfIntegrableFunctions<T> operator+
(
    const integrableFunction<T>& f1,
    const integrableFunction<T>& f2
)
{
  return sumOfIntegrableFunctions<T>(f1, f2);
}

template<class T>
const sumOfIntegrableFunctions<T> operator+
(
  const integrableFunction<T>& f1,
  const sumOfIntegrableFunctions<T>& f2
)
{
  sumOfIntegrableFunctions<T> s(f2);
  s+=f1;
  return s;
}


template<class T1, class T2>
const sumOfIntegrableFunctions<typename outerProduct <T1, T2>::type> operator*
(
    const integrableFunction<T1>& f1,
    const sumOfIntegrableFunctions<T2>& f2
)
{
  sumOfIntegrableFunctions<typename outerProduct <T1, T2>::type> ret;
  for //(size_type i=0; i<f2.size(); i++)
  (typename sumOfIntegrableFunctions<T2>::const_iterator i=f2.begin();
   i!=f2.end(); i++)
    ret += f1 * /*f2[i]*/(*i)();
  return ret;
}

template<class T, class T2>
const sumOfIntegrableFunctions<typename outerProduct <T, T2>::type> operator*
(
    const sumOfIntegrableFunctions<T>& f1,
    const sumOfIntegrableFunctions<T2>& f2
)
{
  sumOfIntegrableFunctions<typename outerProduct <T, T2>::type> ret;
  for (typename sumOfIntegrableFunctions<T>::const_iterator i=f1.begin(); i!=f1.end(); i++)
    for (typename sumOfIntegrableFunctions<T>::const_iterator j=f2.begin(); j!=f2.end(); j++)
      ret += (*i)() * (*j)();
  return ret;
}

/* ------------ operations between integrable functions ----------------- */

template<class T, class T2>
const multipliedIntegrableFunction<T, T2> operator*
(
    const integrableFunction<T>& f1,
    const integrableFunction<T2>& f2
)
{
  return multipliedIntegrableFunction<T, T2>(f1, f2);
}

template<class T>
const dividedByConstantIntegrableFunction<T> operator/
(
    const integrableFunction<T>& f,
    const scalar& c
)
{
  return dividedByConstantIntegrableFunction<T>(f, c);
}
    
template<class T>
const multipliedByConstantIntegrableFunction<T> operator*
(
    const scalar& c,
    const integrableFunction<T>& f
)
{
    return multipliedByConstantIntegrableFunction<T>(c, f);
}
   

template<class T2>
const sumOfIntegrableFunctions<T2> operator*
(
    scalar c,
    const sumOfIntegrableFunctions<T2>& f2
)
{
  sumOfIntegrableFunctions<T2> ret;
  for (typename sumOfIntegrableFunctions<T2>::const_iterator j=f2.begin(); j!=f2.end(); j++)
      ret += c * (*j)();
  return ret;
}

inverseIntegrableFunction inv(const integrableFunction<scalar>&);

}

#endif
