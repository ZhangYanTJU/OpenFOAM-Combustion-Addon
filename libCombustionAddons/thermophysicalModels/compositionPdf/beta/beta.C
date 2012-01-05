#include "beta.H"
#include "addToRunTimeSelectionTable.H"
#include "messageStream.H"

namespace Foam
{

defineTypeNameAndDebug(beta, 0);
addToRunTimeSelectionTable(pdf, beta, dictionary);

const scalar beta::Parameters::epsilon=1e-3;
const scalar beta::Parameters::boundary=1e-3;

beta::Parameters::Parameters(const Parameters& p)
  : integrableFunction<scalar>(p),
    parameterName_(p.parameterName_),
    alpha_(p.alpha_),
    beta_(p.beta_),
    lG_(p.lG_),
    diracDeltas_(p.diracDeltas_),
    degenerated_(p.degenerated_)
{}

beta::Parameters::Parameters
(
    const word& v,
    scalar firstMoment,
    scalar secondMoment
)
  :parameterName_(v)
{
  parameters_.insert(parameterName_, 0);

    if (
      ( secondMoment - firstMoment > SMALL )
        ||
      ( /*secondMoment*/ SMALL < sqr(firstMoment)-secondMoment )
      )
    {
      scalar newsecondMoment=min(max(secondMoment, sqr(firstMoment)), firstMoment);
      WarningIn("beta::Parameters()") << endl
        << "Error: prescribed second moment (" << secondMoment
        << ") is  greater than firstMoment (" << firstMoment
        << ") or smaller than sqr(first moment) (" << sqr(firstMoment) << ")"
          << endl
          <<"Bounding second moment to " <<newsecondMoment <<endl;
          //<< abort(FatalError);
      secondMoment=newsecondMoment;
    }

  
    if (mag(secondMoment - firstMoment*firstMoment) < epsilon)
    {
        alpha_=firstMoment;
        diracDeltas_.append(diracLocation1D(alpha_, 1.0));
        degenerated_=true;
    }
    else if
    (
        mag(1.0 - (secondMoment - firstMoment*firstMoment)
            /(firstMoment*(1.-firstMoment))) < epsilon
    )
    {
        alpha_=1.0-firstMoment;
        beta_=firstMoment;
        diracDeltas_.append(diracLocation1D(0.0, alpha_));
        diracDeltas_.append(diracLocation1D(1.0, beta_));
        degenerated_=true;
    }
    else
    {
        alpha_= (firstMoment*(secondMoment-firstMoment))
                / (firstMoment*firstMoment-secondMoment);
        beta_=((firstMoment-1.0)*(firstMoment-secondMoment))
              / (firstMoment*firstMoment-secondMoment);

        scalar lGab, lGa, lGb;
        lGab = Foam::lgamma(alpha_+beta_);
        lGa  = Foam::lgamma(alpha_);
        lGb  = Foam::lgamma(beta_);
        lG_ = Foam::exp(lGab - lGa - lGb);

        diracDeltas_.append(diracLocation1D(0.0, pow(boundary, alpha_)*lG_/alpha_));
        diracDeltas_.append(diracLocation1D(1.0, pow(boundary, beta_)*lG_/beta_));

        degenerated_=false;
    }

    //Info<<"alpha="<<alpha_<<" beta="<<beta_<<" SLogGamma="<<lG_<<endl;
}

beta::Parameters::~Parameters()
{
}

scalar beta::Parameters::valueAt(const integrationParameter& x) const
{
    // transform integral to leave out boundaries
    double xstar=boundary + (1.0 - 2.0*boundary) * x[parameterName_];
    return interior(xstar) * 1.0/(1.0-2.0*boundary);
}


integrableFunction<scalar>* beta::Parameters::clone() const
{
  return new Parameters(*this);
}

beta::beta(const beta& o)
  : pdf(o)     
{
}

beta::beta(const word& pn)
: pdf(pn)
{
}

beta::beta(const dictionary& dict)
: pdf(dict)
{
}

void beta::setParameters(const Parameters& p)
{
  clear();
  operator+=(p);
  for (label i=0; i<p.diracDeltas_.size(); i++)
    operator+=(diracDelta(parameterName(), p.diracDeltas_[i]));
}

void beta::setParameters(scalar m1, scalar m2)
{
    //Info<<m1<<" "<<m2<<endl;
  setParameters(Parameters(parameterName(), m1, m2));
}

beta::~beta()
{}

autoPtr<pdf> beta::clone() const
{
  return autoPtr<pdf>(new beta(*this));
}

/*
void beta::savePlot(Ostream& f, const Parameters& p)
{
    for (scalar x=1.0/25.0;x<1.0;x+=1.0/25.0)
        f<<x<<" "<<interior(p, x)<<endl;
}
*/

}
