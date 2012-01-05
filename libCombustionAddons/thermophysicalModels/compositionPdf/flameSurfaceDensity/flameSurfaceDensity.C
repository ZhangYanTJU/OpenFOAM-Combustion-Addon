#include "flameSurfaceDensity.H"

#include "messageStream.H"
#include "interpolateXY.H"

#include "OFstream.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(flameSurfaceDensity, 0);
addToRunTimeSelectionTable(pdf, flameSurfaceDensity, dictionary);


    flameSurfaceDensity::Gdelta::Gdelta(const Gdelta& gd)
        : delta_(gd.delta_),
          Gdelta_(gd.Gdelta_),
          a_(gd.a_),
          b_(gd.b_),
          d_(gd.d_)
    {
    }

    flameSurfaceDensity::Gdelta::Gdelta(const dictionary& dict)
        : delta_(readScalar(dict.lookup("delta")))
    {

        typedef List<Tuple2<scalar,scalar> > L;
        interpolationTable<scalar> X_C(dict);
        label Np=X_C.size();

        // normalize
        scalar maxC=0.0;
        forAll(X_C, I)
            if (X_C[I].second()> maxC) maxC=X_C[I].second();
        forAll(X_C, I)
            X_C.L::operator[](I).second()/=maxC;

        // compute gradient of progress variable
        interpolationTable<scalar> X_GradC(X_C);
        // first and last point
        X_GradC.L::operator[](0).second()=(X_C[1].second()-X_C[0].second())/(X_C[1].first()-X_C[0].first());
        X_GradC.L::operator[](Np-1).second()=(X_C[Np-1].second()-X_C[Np-2].second())/(X_C[Np-1].first()-X_C[Np-2].first());
        for (label I=1; I<Np-1; I++)
            X_GradC.L::operator[](I).second()=(X_C[I+1].second()-X_C[I-1].second())/(X_C[I+1].first()-X_C[I-1].first());
if (debug)
{
OFstream f("X_GradC");
f<<X_GradC<<endl;
}

        // compute filtered gradient of progress variable
        interpolationTable<scalar> X_GradC_filtered(X_C);
        forAll(X_GradC_filtered, I)
        {
            scalar x=X_GradC_filtered[I].first();
            scalar Ix=0.0;
            X_GradC_filtered.L::operator[](I).second()=0.0;
            for (label J=0; J<X_GradC.size()-1; J++)
            {
                scalar x0=X_GradC[J].first();
                scalar x1=X_GradC[J+1].first();
                bool inside=false;
                if (( mag(x0-x)<0.5*delta_ ) && ( mag(x1-x)<0.5*delta_ ))
                {
                    // completely inside filter width
                    inside=true;
                } else 
                if (( mag(x0-x)>0.5*delta_ ) && ( mag(x1-x)<0.5*delta_ ))
                {
                   // intervall cuts left boundary
                   x0=x-0.5*delta_;
                   inside=true;
                } else 
                if (( mag(x0-x)<0.5*delta_ ) && ( mag(x1-x)>0.5*delta_ ))
                {
                   // intervall cuts right boundary
                   x1=x+0.5*delta_;
                   inside=true;
                } else 
                if (( x-x0>0.5*delta_ ) && ( x1-x>0.5*delta_ ))
                {
                   // intervall encloses boundary
                   x0=x-0.5*delta_;
                   x1=x+0.5*delta_;
                   inside=true;
                }
                if (inside)
                {
                    scalar d=x1-x0;
                    Ix+=d;
                    scalar y0=X_GradC(x0);
                    scalar y1=X_GradC(x1);
                    X_GradC_filtered.L::operator[](I).second()+=0.5*d*(y1+y0);
                }
            }
            if (Ix<SMALL)
             X_GradC_filtered.L::operator[](I).second()=X_GradC[I].second();
            else
             X_GradC_filtered.L::operator[](I).second()/=Ix;
        }
if (debug)
{
OFstream f("X_GradC_filtered");
f<<X_GradC_filtered<<endl;
}


        Gdelta_.setSize(X_GradC_filtered.size());
        forAll(Gdelta_, I)
        {
            Gdelta_.L::operator[](I).first()=X_C[I].second();
            Gdelta_.L::operator[](I).second()=X_GradC_filtered[I].second();
        }

if (debug)
{
OFstream f("Gdelta");
f<<Gdelta_<<endl;
}

        // integrate to get a, b and d
        a_=0.0;
        b_=0.0;
        d_=0.0;

        for (label i=1; i<Gdelta_.size(); i++)
        {
            scalar dpv=Gdelta_[i].first()-Gdelta_[i-1].first();
            a_+=0.5*( 
                (1.0/Gdelta_[i].second())
                +
                (1.0/Gdelta_[i-1].second()) 
                    )*dpv;
            b_+=0.5*( 
                (Gdelta_[i].first()/Gdelta_[i].second())
                +
                (Gdelta_[i-1].first()/Gdelta_[i-1].second()) 
                   )*dpv;
            d_+=0.5*( 
                (sqr(Gdelta_[i].first())/Gdelta_[i].second())
                +
                (sqr(Gdelta_[i-1].first())/Gdelta_[i-1].second()) 
                   )*dpv;
        }

        Info<<"a="<<a_<<", b="<<b_<<", d="<<d_<<endl;
    };
 
    scalar flameSurfaceDensity::Gdelta::operator()(scalar pvv) const
    {
        return Gdelta_(pvv);
    };

    flameSurfaceDensity::Parameters::Parameters(const Parameters& p)
        : integrableFunction<scalar>(p),
          parameterName_(p.parameterName_),
          alpha_(p.alpha_),
          beta_(p.beta_),
          sigma_(p.sigma_),
          diracDeltas_(p.diracDeltas_),
          degenerated_(p.degenerated_),
          gdelta_(p.gdelta_)
    {}


    flameSurfaceDensity::Parameters::Parameters
    (
        const word& pv,
        scalar firstMoment, scalar secondMoment,
        const Gdelta& gd
    )
        : parameterName_(pv),
          gdelta_(gd)
    {
        parameters_.insert(parameterName_, 0);

        if (
            ( secondMoment - firstMoment > SMALL )
            ||
            ( /*secondMoment*/ SMALL < sqr(firstMoment)-secondMoment )
        )
        {
            scalar newsecondMoment=min(max(secondMoment, sqr(firstMoment)), firstMoment);
            WarningIn("flameSurfaceDensity::Parameters()") << endl
                << "Error: prescribed second moment (" << secondMoment
                << ") is  greater than firstMoment (" << firstMoment
                << ") or smaller than sqr(first moment) (" << sqr(firstMoment) << ")"
                << endl
                <<"Bounding second moment to " <<newsecondMoment <<endl;
            //<< abort(FatalError);
            secondMoment=newsecondMoment;
        }

        scalar a=gdelta_.a();
        scalar b=gdelta_.b();
        scalar d=gdelta_.d();

        if 
            (
                (secondMoment < ((a-d)*firstMoment+d-b)/(a-b))
                ||
                (secondMoment < d*firstMoment/b)
            )
            throw string("Variance too low for FSD PDF");

        alpha_=(b-b*secondMoment+a*(secondMoment-firstMoment)+(firstMoment-1)*d)/(b-d);
        beta_=(b*secondMoment-firstMoment*d)/(b-d);
        sigma_=(firstMoment-secondMoment)/(b-d);

        Info<<"alpha="<<alpha_<<", beta="<<beta_<<", sigma="<<sigma_<<endl;

        diracDeltas_.append(diracLocation1D(0.0, alpha_));
        diracDeltas_.append(diracLocation1D(1.0, beta_));

    }



    scalar flameSurfaceDensity::Parameters::valueAt(const integrationParameter& x) const
    {
        return interior(x[parameterName_]);
    }

    integrableFunction<scalar>* flameSurfaceDensity::Parameters::clone() const
    {
        return new Parameters(*this);
    }
    /*
    flameSurfaceDensity::Parameters::Parameters
    (
        scalar firstMoment, 
        scalar secondMoment,
        scalar delta,
        const chemistryTable& table,
        const word& pv
    )
        : Gdelta_
          (
              "Gdelta", pv, "<grad("&pv&")>",
              scalarField(table.nElements(0), 0.0),
              scalarField(table.nElements(0), 0.0)
          ),
        singularities_(0, 0.0)
    {
        if (table.dimensionality()!=1)
            FatalErrorIn("flameSurfaceDensity::Parameters::Parameters()")
		<< "only one-dimensional FGM tables are supported."
                    << abort(FatalError);			      



        // filter gradient profile

        scalarField& c=Gdelta_.x();
        scalarField& G_filt=Gdelta_.y();

        label Ix=table.indexOf("x");
        label Igrad=table.indexOf("grad("&pv&")");
        for (label i=0;i<table.nElements(0);i++)
        {
            scalar Gfilt=0.0;

            // to the right
            scalar x0=table.const_access(&i)[Ix];
            scalar xjm1=x0;
            for (label j=i+1;j<table.nElements(0);j++)
            {
                scalar xj=table.const_access(&j)[Ix];
                label jm1=j-1;
                if (xj<x0+0.5*delta)
                {
                    Gfilt+=0.5*
                        (
                            table.const_access(&j)[Igrad]
                            +
                            table.const_access(&jm1)[Igrad]
                        )*(xj-xjm1);
                            
                } else
                {
                    scalar x=x0+0.5*delta;
                    scalar G=
                        (
                        table.const_access(&j)[Igrad]*(x-xjm1)
                        +
                        table.const_access(&jm1)[Igrad]*(xj-x)
                        )/(xj-xjm1);

                    Gfilt+=0.5*
                        (G + table.const_access(&jm1)[Igrad])*(x-xjm1);

                    break;
                }
                xjm1=xj;
            }            

            // to the left
            x0=table.const_access(&i)[Ix];
            scalar xjp1=x0;
            for (label j=i-1;j>=0;j--)
            {
                label jp1=j+1;
                scalar xj=table.const_access(&j)[Ix];
                if (xj>x0-0.5*delta)
                {
                    Gfilt+=0.5*
                        (
                            table.const_access(&jp1)[Igrad]
                            +
                            table.const_access(&j)[Igrad]
                        )*(xjp1-xj);
                            
                } else
                {
                    scalar x=x0-0.5*delta;
                    scalar G=
                        (
                        table.const_access(&jp1)[Igrad]*(x-xj)
                        +
                        table.const_access(&j)[Igrad]*(xjp1-x)
                        )/(xjp1-xj);

                    Gfilt+=0.5*
                        (table.const_access(&jp1)[Igrad] + G)*(xjp1-x);

                    break;
                }
                xjp1=xj;
            }

            c[i]=table.startValue(0)+table.deltaValue(0)*scalar(i);
            G_filt[i]=Gfilt/delta;
        }


	scalar a=0.0, b=0.0, d=0.0;
	for(label i=1;i<c.size();i++)
	{
	  a+=0.5*( (1.0/G_filt[i-1]) + (1.0/G_filt[i]) )*(c[i] - c[i-1]);
	  b+=0.5*( (c[i-1]/G_filt[i-1]) + (c[i]/G_filt[i]) )*(c[i] - c[i-1]);
	  d+=0.5*( (sqr(c[i-1])/G_filt[i-1]) + (sqr(c[i])/G_filt[i]) )
	      *(c[i] - c[i-1]);
	}

	scalar S=max
	  (
	   (secondMoment-sqr(firstMoment))/
	   (firstMoment*(1.0-firstMoment)),
	   max
	   (
	    1.0-((b-d)/(a-b))/firstMoment,
	    1.0-((b-d)/b)/(1.0-firstMoment)
           )
	  );

	alpha_=(1.0-firstMoment) - ((a-b)/(b-d))*(1.0-S)
	  *firstMoment*(1.0-firstMoment);

	beta_=firstMoment - (b/(b-d))*(1.0-S)
	  *firstMoment*(1.0-firstMoment);

	sigma_=(1.0/(b-d))*(1.0-S)*firstMoment*(1.0-firstMoment);

	Info<<firstMoment<<" "<<S<<endl;
	Info<<a<<" "<<b<<" "<<d<<" "<<alpha_<<" "<<beta_<<" "<<sigma_<<endl;
    }

    scalar flameSurfaceDensity::Parameters::Gdelta(scalar c) const
    {
        return interpolateXY(c, Gdelta_.x(), Gdelta_.y());
    }


    void flameSurfaceDensity::savePlot(Ostream& f, const Parameters& p)
    {
        for (scalar x=1.0/25.0;x<1.0;x+=1.0/25.0)
            f<<x<<" "<<interior(p, x)<<endl;
    }

    const flameSurfaceDensity::Parameters* 
       flameSurfaceDensity::Ipdf_=NULL;
    const integrableFunction* flameSurfaceDensity::I_=NULL;

    double flameSurfaceDensity::integrand(double x)
    {
        return 
	  (
	   I_->valueAt(x)
	    *
	   flameSurfaceDensity::interior(*Ipdf_, x)
	  );
    }

    scalar flameSurfaceDensity::integrate
    (
        const Parameters& p,
        const integrableFunction& func
    )
    {        
        /*
      if (p.singleDirac_)
	{
	  return table.lookup(p.alpha_, i);
	}
      else if (p.boundaryDirac_)
	{
	  return 
	    p.alpha_*table.lookup(0, i)
	    +
	    p.beta_*table.lookup(1, i);
	}
      else
	{	  
	  Ipdf_=&p;
	  Itable_=&table;
	  Iidx_=i;

	  scalar epsilon=1e-6;
	  DefQuad integrator
	    (
	     &integrand, 
	     epsilon, 1.0-epsilon
	    );

	  int ier, neval;
	  double abserr;
	  double I=
	    integrator.do_integrate(ier, neval, abserr)
	    +table.lookup(0, i)*pow(epsilon, p.alpha_)*p.lG_/p.alpha_
	    +table.lookup(1, i)*pow(epsilon, p.beta_)*p.lG_/p.beta_;

	  if (ier>5)
	    {
	      FatalErrorIn("beta::integrateTable()")
		<< "numerical integration of chemistry table over"
		<< " beta pdf failed."<<endl
		<< " IER="<<ier<<endl
		<< " NEVAL="<<neval<<endl
		<< " ABSERR="<<abserr<<endl
		<< " I="<<I<<endl
		<< abort(FatalError);			      
	    }
	  Info<<ier<<" "<<neval<<" "<<abserr<<" "<<I<<endl;
	  return I;
            */
    //    }

    flameSurfaceDensity::flameSurfaceDensity(const flameSurfaceDensity& fsd)
        : beta(fsd),
          Gdelta_(fsd.Gdelta_)
    {
    }

    flameSurfaceDensity::flameSurfaceDensity(const word& pn)
        : beta(pn)
    {
    }

    flameSurfaceDensity::flameSurfaceDensity(const dictionary& dict)
        : beta(dict)
    {
        Gdelta_.reset(new Gdelta(dict));
    }

    flameSurfaceDensity::~flameSurfaceDensity()
    {
    }

    /*
    void flameSurfaceDensity::setParameters(const Parameters&)
    {
        clear();
    }
        */

    void flameSurfaceDensity::setParameters(scalar m1, scalar m2)
    {
        clear();
        try
        {
            Parameters fsd_p(parameterName(), m1, m2, Gdelta_);
            operator+=(fsd_p);
            for (label i=0; i<fsd_p.diracDeltas_.size(); i++)
                operator+=(diracDelta(parameterName(), fsd_p.diracDeltas_[i]));
        }
        catch (...)
        {
            beta::setParameters(m1, m2);
        }

    }

    autoPtr<pdf> flameSurfaceDensity::clone() const
    {
        return autoPtr<pdf>(new flameSurfaceDensity(*this));
    }
}

