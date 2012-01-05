#include "clippedGaussian.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

    defineTypeNameAndDebug(clippedGaussian, 0);
    addToRunTimeSelectionTable(pdf, clippedGaussian, dictionary);

  #include "parameterMapData.H"

  clippedGaussianParameterMappingTable
    ::clippedGaussianParameterMappingTable()
      : MultidimensionalLookupTable<clippedGaussianParameterMapping>
         (2, numbers, starts, deltas)
  {
    delete elements;
    elements=(clippedGaussianParameterMapping*) &clipgaussmap;
  }

  clippedGaussianParameterMappingTable
    ::~clippedGaussianParameterMappingTable()
  {
    elements=new clippedGaussianParameterMapping;
  }

  void clippedGaussianParameterMapping::operator+=
  (
   const clippedGaussianParameterMapping& o
  )
  {
    mug_+=o.mug_;
    sigmag_+=o.sigmag_;
  }
      
  void clippedGaussianParameterMapping::operator/=
  (
   scalar s
  )
  {
    mug_/=s;
    sigmag_/=s;
  }

  clippedGaussianParameterMapping operator*
  (
   scalar s,
   const clippedGaussianParameterMapping& o
  )
  {
    clippedGaussianParameterMapping sn(o);
    sn.mug_*=s;
    sn.sigmag_*=s;
    return sn;
  }


  scalar clipexp(scalar x)
  {
    if (x<100) return exp(x);
    else {
      return 1e10;
    }
  }

  scalar clippow(scalar x, scalar y)
  {
    if (x==0.0) return 1e10;
    else return pow(x, y);
  }

  clippedGaussianParameterMappingTable 
     clippedGaussian::Parameters::mapping_;

    clippedGaussian::Parameters::Parameters
    (
     const Parameters& o
    )
    : integrableFunction<scalar>(o),
      parameterName_(o.parameterName_),
      mug_(o.mug_),
      sigmag_(o.sigmag_),
      Al_(o.Al_), 
      Ar_(o.Ar_),
      degenerated_(o.degenerated_),
      diracDeltas_(o.diracDeltas_)
    {}


    clippedGaussian::Parameters::Parameters
    (
        const word& v,
        scalar firstMoment, 
        scalar secondMoment
    )
        : parameterName_(v)
    {
        parameters_.insert(parameterName_, 0);
        
        if (
            ( secondMoment - firstMoment > SMALL )
            ||
            ( /*secondMoment*/ SMALL < sqr(firstMoment)-secondMoment )
        )
        {
            scalar newsecondMoment=min(max(secondMoment, sqr(firstMoment)), firstMoment);
            WarningIn("clippedGaussian::Parameters()") << endl
                << "Error: prescribed second moment (" << secondMoment
                << ") is  greater than firstMoment (" << firstMoment
                << ") or smaller than sqr(first moment) (" << sqr(firstMoment) << ")"
                << endl
                <<"Bounding second moment to " <<newsecondMoment <<endl;
            //<< abort(FatalError);
            secondMoment=newsecondMoment;
        }


        if ( mag(secondMoment-firstMoment*firstMoment) < 1e-5) 
        {
            mug_=firstMoment;
            diracDeltas_.append(diracLocation1D(mug_, 1.0));
            degenerated_=true;
            //Info<<"mu="<<mug_<<endl;
	} 
        else
        {
            degenerated_=false;
            scalarField p(2,0);
            p[1]=firstMoment;
            p[0]=(secondMoment-firstMoment*firstMoment)/
                (firstMoment*(1.0-firstMoment));
            
            const clippedGaussianParameterMapping& m=
                mapping_.lookup(p);
            
            mug_=m.mug_;
            sigmag_=m.sigmag_;
            
            Al_=clippedGaussian::Al(mug_, sigmag_);
            Ar_=clippedGaussian::Ar(mug_, sigmag_);
            
            diracDeltas_.append(diracLocation1D(0.0, 0.5*Al_));
            diracDeltas_.append(diracLocation1D(1.0, 0.5*Ar_));

            //Info<<"mu_g="<<mug_<<" sigma_g="<<sigmag_<<endl;
            //Info<<"Al="<<Al_<<" Ar="<<Ar_<<endl;
        }
    }

    clippedGaussian::Parameters::~Parameters()
    {
    }

    scalar clippedGaussian::Parameters::valueAt
    (
        const integrationParameter& x
    ) const
    {
        return interior(x[parameterName_]);
    }


    integrableFunction<scalar>* clippedGaussian::Parameters::clone() const
    {
        return new Parameters(*this);
    }

    clippedGaussian::clippedGaussian(const clippedGaussian& o)
        : pdf(o)     
    {
    }
    
    clippedGaussian::clippedGaussian(const word& pn)
        : pdf(pn)
    {
    }

    clippedGaussian::clippedGaussian(const dictionary& dict)
        : pdf(dict)
    {
    }

    void clippedGaussian::setParameters(const Parameters& p)
    {
        clear();
        operator+=(p);
        for (label i=0; i<p.diracDeltas_.size(); i++)
            operator+=(diracDelta(parameterName(), p.diracDeltas_[i]));
    }

    void clippedGaussian::setParameters(scalar m1, scalar m2)
    {
        //Info<<m1<<" "<<m2<<endl;
        setParameters(Parameters(parameterName(), m1, m2));
    }

    clippedGaussian::~clippedGaussian()
    {}
    
    autoPtr<pdf> clippedGaussian::clone() const
    {
        return autoPtr<pdf>(new clippedGaussian(*this));
    }
    

    /*
    double clippedGaussian::integrand(double x)
    {
        return 
            (
            I_->valueAt(x)
            *
            clippedGaussian::interior(*Ipdf_, x)
            )
            + 
            0.5*Ipdf_->Al()
             *I_->valueAt(0.)
            + 
            0.5*Ipdf_->Ar()
             *I_->valueAt(1.);

    }

  void clippedGaussian::savePlot(Ostream& f, const Parameters& p)
  {
    for (scalar x=1.0/25.0;x<1.0;x+=1.0/25.0)
      f<<x<<" "<<interior(p, x)<<endl;
  }

    scalar clippedGaussian::integrate
    (
        const Parameters& p,
        const integrableFunction& integ
    )
    {
      if (p.singleDirac_)
	{
	  return integ.valueAt(p.mug_);
	}
      else
	{
	  Ipdf_=&p;
	  I_=&integ;
	  DefQuad integrator(&integrand, 0, 1);
	  int ier, neval;
	  double abserr=1e-9;
	  double I=integrator.do_integrate(ier, neval, abserr);
	  //Info<<ier<<" "<<neval<<" "<<abserr<<endl;
	  return I;
	}
    }

        */

}
