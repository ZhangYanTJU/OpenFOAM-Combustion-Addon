#include "chemistryTable.H"

namespace Foam
{


chemistryTableEntryIntegrand::chemistryTableEntryIntegrand
(
  const chemistryTableEntryIntegrand& o
)
  : integrableFunction<scalar>(o),
    table_(o.table_),
    contentVariable_(o.contentVariable_)
{
}

chemistryTableEntryIntegrand::chemistryTableEntryIntegrand
(
  const chemistryTable& tab,
  const word& entryName
)
  : table_(tab),
    contentVariable_(table_.indexOf(entryName))
{
  for (HashTable<label,word>::const_iterator it=table_.progressVariableDescription().begin();
       it!=table_.progressVariableDescription().end(); it++)
    this->parameters_.insert(it.key(), 0);
}

scalar chemistryTableEntryIntegrand::valueAt(const integrationParameter& x) const
{
  scalarField cF(table_.dimensionality(), 0.0);

  /*
  for (integrationParameter::const_iterator it=x.begin();
       it!=x.end(); it++)
    cF[table_.indexOfPV(it.key())]=it();
  */
  for (HashTable<label,word>::const_iterator i=table_.progressVariableDescription().begin();
       i!=table_.progressVariableDescription().end(); i++)
  {
    cF[i()]=x[i.key()];
    //Info<<i.key()<<"("<<i()<<")="<<x[i.key()]<<endl;
  }

  return
    table_.MultidimensionalLookupTable<chemicalSystemState>
    ::lookup(cF)[contentVariable_];
}

integrableFunction<scalar>* chemistryTableEntryIntegrand::clone() const
{
  return new chemistryTableEntryIntegrand(*this);
}



ListDescription::ListDescription()
{
}

ListDescription::ListDescription(const ListDescription& o)
    : HashTable<label,word>(o)
{
}

ListDescription::ListDescription(const HashTable<label,word>& o)
    : HashTable<label,word>(o)
{
}

void ListDescription::appendEntry(const word& name)
{
    label lastidx=-1;
    for (const_iterator it=begin(); it!=end(); it++)
        if (it()>lastidx) lastidx=it();
    
    insert(name, ++lastidx);
}




chemistryTable::chemistryTable
(
    const label dim,
    const label* n,
    const scalar* s,
    const scalar* d,
    const ListDescription& content,
    const ListDescription& pvdesc
)
    : MultidimensionalLookupTable<chemicalSystemState>(dim, n, s, d),
      pvdesc_(pvdesc),
      contentdesc_(content)//,
    //cursor_(dim, 0.0),
    //integrationVariable_(0)
{
}

chemistryTable::chemistryTable
(
Istream& f
)
    : MultidimensionalLookupTable<chemicalSystemState>(f),
      pvdesc_(f),
      contentdesc_(f)//,
//cursor_(pvdesc_.size(), 0.0),
//integrationVariable_(0)
{
}

chemistryTable::~chemistryTable()
{
}


chemistryTableEntryIntegrand chemistryTable::entry(const word& entryName) const
{
  return chemistryTableEntryIntegrand(*this, entryName);
}

word chemistryTable::nameOf(label i) const
{
  for (ListDescription::const_iterator it=contentdesc_.begin();
       it!=contentdesc_.end(); it++)
    if (it()==i) return it.key();

  return "";
}

word chemistryTable::nameOfPV(label i) const
{
  for (ListDescription::const_iterator it=pvdesc_.begin();
       it!=pvdesc_.end(); it++)
    if (it()==i) return it.key();

  return "";
}

scalar chemistryTable::findPV
(
    label pvidx, 
    const List<scalar>& cursor, 
    label contentidx, 
    scalar contentval
) const
{
    List<scalar> tmpcur(cursor);

    // find via bisection
    scalar lo, hi, mid;
    scalar diff_lo, diff_hi, diff_mid;
    lo=startValue(pvidx);
    hi=lo+deltaValue(pvidx)*(nElements(pvidx)-1);
    
    //Info<<"CUR:"<<cursor<<endl;
    label count=0;
    do 
    {
        count++;

        mid=0.5*(hi+lo);        
        //Info<<lo<<" "<<mid<<" "<<hi<<endl;

        tmpcur[pvidx]=lo;
        diff_lo=lookup(tmpcur)[contentidx] - contentval;
        if (mag(diff_lo)<1e-3) 
        {
            //Info<<"RETURN"<<endl;
            return lo;
        }
        
        tmpcur[pvidx]=hi;
        diff_hi=lookup(tmpcur)[contentidx] - contentval;
        if (mag(diff_hi)<1e-3) 
        {
            //Info<<"RETURN"<<endl;
             return hi;
        }
        tmpcur[pvidx]=mid;
        diff_mid=lookup(tmpcur)[contentidx] - contentval;
        if (mag(diff_mid)<1e-3) 
        {
            //Info<<"RETURN"<<endl;
            return mid;
        }
        
        //Info<<diff_lo<<" "<<diff_mid<<" "<<diff_hi<<endl;

        
        if ( sign(diff_hi) != sign(diff_mid) ) lo=mid; 
        else
        if ( sign(diff_mid) != sign(diff_lo) ) hi=mid;
        else
        {
            WarningIn("chemistryTable::findPV()")
                << "no solution in any interval: target="<<contentval<<" cursor="<<cursor<< endl
                    <<"target="<<contentval<<" cursor: "<<cursor<<endl
                    <<lo<<" "<<mid<<" "<<hi<<endl
                    <<diff_lo<<" "<<diff_mid<<" "<<diff_hi<<endl
                ;//<< exit(FatalError);
            if (mag(diff_lo)<mag(diff_hi)) 
            { Info<<"returning lo val"<<endl; return lo; }
            else 
            { Info<<"returning hi val"<<endl; return hi; }
        }
        
        if (count>1000)
        {
            FatalErrorIn("chemistryTable::reverseLookup()")
                << "bisection did not converge within 1000 iterations."<< endl
                    << abort(FatalError);
        }
    } 
    while (mag(diff_mid)>1e-3);
    //Info<<"RETURN"<<endl;
    return mid;
}


//#define REVERSELOOKUP_NEWTON
#undef REVERSELOOKUP_NEWTON

#ifdef REVERSELOOKUP_NEWTON


/*
  Table Search via Newton-Method
    */
chemicalSystemState chemistryTable::reverseLookup
(
    label pvidx, 
    const List<scalar>& cursor, 
    label contentidx, 
    scalar contentval,
    scalar& pvval,
    string& error
) const
{
    List<scalar> tmpcur(cursor);
    
    // find via newton method

    scalar minpv=startValue(pvidx);
    scalar maxpv=minpv+deltaValue(pvidx)*(nElements(pvidx)-1);

    scalar delta=deltaValue(pvidx);

    scalar estim=tmpcur[pvidx];
    if (estim<minpv || estim>maxpv) estim=0.5*(minpv+maxpv);

    scalar pvnew=estim;
    scalar pvtol=max(1e-3,mag(1e-3*estim));

    int iter=0;
    do
    {
        estim=pvnew;

        scalar lo=estim-delta;
        if (lo<minpv) lo=minpv;
        scalar hi=estim+delta;
        if (hi>maxpv) hi=maxpv;

        tmpcur[pvidx]=lo;
        chemicalSystemState vlo=lookup(tmpcur);
        tmpcur[pvidx]=hi;
        chemicalSystemState vhi=lookup(tmpcur);
        scalar deriv=(vhi[contentidx] - vlo[contentidx]) / (hi-lo);

        if (mag(deriv)<1e-3) 
        {
            //Info<<"Derivative zero!: "<<deriv<<endl;
            // Very probably we hit a row in the table where the quantity is constant
            // Doesn't matter which value to return !?
            pvval=estim;
            break;
        }

        tmpcur[pvidx]=estim;
        chemicalSystemState vest=lookup(tmpcur);

        pvnew=estim - 0.98*((vest[contentidx] - contentval)/deriv);

        if (iter>10) Info<<"Iter #"<<iter<<": "<<estim<<" "<<pvnew<<" "<<deriv<<" "<<lo<<" "<<hi<<endl;

        if (pvnew>maxpv) 
        {
            pvval=maxpv;
            break;
        }
        if (pvnew<minpv) 
        {
            pvval=minpv;
            break;
        }

        if (iter++ > 100)
        {
            FatalErrorIn
                (
                    "chemistryTable::reverseLookup()"
                )   << "Maximum number of iterations exceeded"
                    << abort(FatalError);
        }

    }
    while (mag(pvnew - estim) > pvtol);
    pvval=pvnew;

    tmpcur[pvidx]=pvval;
    return lookup(tmpcur);
}


#else

/*
  Table Search via bisection method
    */
chemicalSystemState chemistryTable::reverseLookup
(
    label pvidx, 
    const List<scalar>& cursor, 
    label contentidx, 
    scalar contentval,
    scalar& pvval,
    string& error
) const
{
    List<scalar> tmpcur(cursor);

    error="";

    // find via bisection
    scalar lo, hi, mid;
    scalar diff_lo, diff_hi, diff_mid;
    chemicalSystemState vlo, vmid, vhi;
    lo=startValue(pvidx);
    hi=lo+deltaValue(pvidx)*(nElements(pvidx)-1);
    
    tmpcur[pvidx]=lo;
    vlo=lookup(tmpcur);
    diff_lo=vlo[contentidx] - contentval;
    tmpcur[pvidx]=hi;
    vhi=lookup(tmpcur);
    diff_hi=vhi[contentidx] - contentval;

    if (sign(diff_lo) == sign(diff_hi))
    {
        // no solution in interval, return boundary value
	//Info<<"Out of bounds for cv("<<contentidx<<")="<<contentval<<": "<<diff_lo<<" "<<diff_hi<<endl;

        if (mag(diff_lo)<mag(diff_hi))
        {
            pvval=lo;
            return vlo;
        }
        else
        {
            pvval=hi;
            return vhi;
        }
    }


    label count=0;
    do 
    {
        count++;
        
        tmpcur[pvidx]=lo;
        vlo=lookup(tmpcur);
        diff_lo=vlo[contentidx] - contentval;
        if (mag(diff_lo)<1e-3) 
        {
            pvval=lo;
            return vlo;
        }
        
        tmpcur[pvidx]=hi;
        vhi=lookup(tmpcur);
        diff_hi=vhi[contentidx] - contentval;
        if (mag(diff_hi)<1e-3)
        {
            pvval=hi;
            return vhi;
        }
                    
        mid=0.5*(hi+lo);
        tmpcur[pvidx]=mid;
        vmid=lookup(tmpcur);
        diff_mid=vmid[contentidx] - contentval;
        if (mag(diff_mid)<1e-3)
        {
            pvval=mid;
            return vmid;
        }
        
        /*
        Info<<lo<<" "<<mid<<" "<<hi<<endl;
        Info<<diff_lo<<" "<<diff_mid<<" "<<diff_hi<<endl;
            */ 
        
        if ( sign(diff_lo) != sign(diff_mid) ) hi=mid; 
        else
        if ( sign(diff_mid) != sign(diff_hi) ) lo=mid;
        else
        {
            FatalErrorIn("chemistryTable::reverseLookup()")
                //OStringStream errmsg;
                //errmsg
                << "no solution in any interval while searching for "<<contentval
                    << endl
                    <<lo<<" "<<mid<<" "<<hi<<endl
                    <<vlo[contentidx]<<" "<<vmid[contentidx]<<" "<<vhi[contentidx]<<endl
                    <<diff_lo<<" "<<diff_mid<<" "<<diff_hi<<endl
                << exit(FatalError);
            /*
            if (mag(diff_lo)<mag(diff_hi)) 
            { 
                errmsg<<"returning lo val "<<lo<<endl; 
                error=errmsg.str(); 
                pvval=lo; 
                return vlo; 
            }
            else 
            { 
                errmsg<<"returning hi val "<<hi<<endl; 
                error=errmsg.str(); 
                pvval=hi; 
                return vhi; 
            }
                */
        }
        
        if (count>1000)
        {
            FatalErrorIn("chemistryTable::reverseLookup()")
                << "bisection did not converge within 1000 iterations."<< endl
                    << abort(FatalError);
        }
    } 
    while (mag(diff_mid)>1e-3);
    pvval=mid;
    return vmid;
}

#endif

graph chemistryTable::slice
(
    const List<scalar>& cursor, 
    const word& along, 
    //const List<word>& entries
    const word& entry
) const
{

    label ialong=indexOfPV(along);
    List<scalar> c=cursor;

    scalarField x(nElements(ialong), 0.0);
    scalarField val(nElements(ialong), 0.0);

    for (label i=1;i<nElements(ialong);i++) 
    {
        x[i]=valueAt(ialong, i);
        c[ialong]=x[i];
        chemicalSystemState state=lookup(c);
        val[i]=state[indexOf(entry)];
    }

    return graph("slice", along, entry, x, val);

}



Ostream& operator<<(Ostream& f, const chemistryTable& tab)
{
  tab.MultidimensionalLookupTable<chemicalSystemState>::write(f);
  f<<token::SPACE;
  f<<tab.pvdesc_;
  f<<token::SPACE;
  f<<tab.contentdesc_;
  return f;
}

}
