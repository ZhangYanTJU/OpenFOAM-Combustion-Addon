#include "rectangle.H"
#include "octave/config.h"
#include "octave/Quad.h"
#include "messageStream.H"

namespace Foam
{

    rectangle::Parameters::Parameters
    (
        scalar firstMoment, 
        scalar secondMoment
    )
    {
        a_=firstMoment;
        b_=min
	  (
	   Foam::sqrt(12.0*max(0,secondMoment-a_*a_)),
	   min
	   (
	    2.*a_,
	    2.*(1.-a_)
           )
	  );
	if (b_ < 1e-3) singleDirac_=true;
	else singleDirac_=false;
        Info<<"a="<<a_<<" b="<<b_<<endl;
    }


    const rectangle::Parameters* rectangle::Ipdf_=NULL;
    const integrableFunction* rectangle::I_=NULL;

    double rectangle::integrand(double x)
    {
        return 
            (
                (x>(Ipdf_->a() - 0.5*Ipdf_->b()))
                &&
                (x<(Ipdf_->a() + 0.5*Ipdf_->b())) ?
                I_->valueAt(x)
                /Ipdf_->b()
                :0.0
            );
    }

    scalar rectangle::integrate
    (
        const Parameters& p,
        const integrableFunction& func
    )
    {
      if (p.singleDirac_)
	{
	  return func.valueAt(p.a_);
	}
      else
	{	  
	  Ipdf_=&p;
	  I_=&func;
	  DefQuad integrator
	    (
	     &integrand, 
	     p.a() - 0.5*p.b(),
	     p.a() + 0.5*p.b()
	    );
	  int ier, neval;
	  double abserr;
	  double I=integrator.do_integrate(ier, neval, abserr);
	  if (ier>2)
	    {
	      FatalErrorIn("rectangle::integrateTable()")
		<< "numerical integration of chemistry table over"
		<< " rectangle pdf failed."<<endl
		<< " IER="<<ier<<endl
		<< " NEVAL="<<neval<<endl
		<< " ABSERR="<<abserr<<endl
		<< abort(FatalError);			      
	    }
	  Info<<ier<<" "<<neval<<" "<<abserr<<endl;
	  return I;
	}
    }

}
