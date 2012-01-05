/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*----------------------------------------------------------------------------*/

#include "integrator.H"


#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>


#define EPSREL 1e-3
#define EPSABS 1e-12
#define VERBOSE 2
#define LAST 4
#define MINEVAL 0
#define MAXEVAL 50000
#define KEY 0
 
#define NSTART 1000
#define NINCREASE 500
  
#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN dim
#define NEXTRA 0


#include <stdio.h>
#include <math.h>

#include "cuba.h"


//#define DEBUG


namespace Foam
{

    template<class T>
    integrationParameter integrator<T>::X__;

    template<class T>
        integrationParameter integrator<T>::constants__;

    template<class T>
    const integrableFunction<T>* integrator<T>::F__=NULL;
    
    template<class T>
    double integrator<T>::gsl_integrand(double x, void * params)
    {
      X__.setValuesFromArray(&x);

      T res = F__->valueAt(integrationParameter(X__, constants__));

      scalar sres;
      toArray<T>::convert(res, &sres);

      return sres;
    }
    
    template<class T>
    double integrator<T>::octaveIntegrand(double x)
    {
      X__.setValuesFromArray(&x);

      T res = F__->valueAt(integrationParameter(X__, constants__));

      scalar sres;
      toArray<T>::convert(res, &sres);

      return sres;
    }
    
    template<class T>
    int integrator<T>::Integrand
    (
        const int *ndim, const double xx[],
        const int *ncomp, double ff[], void *userdata
    )
    {

      X__.setValuesFromArray(xx);

      T res = F__->valueAt(integrationParameter(X__, constants__));

      toArray<T>::convert(res, &ff[0]);

      return 0;
    }

    template<class T>
    T integrator<T>::integrate
    (
        const integrableFunction<T>& F,
        const integrationParameter& parameters,
        const integrationParameter& constants
    )
    {
      //X__=integrationParameter(parameters, constants),
      X__=parameters;
      constants__=constants;
      F__=&F;
      
      T finalresult=F.valueAt(integrationParameter(X__, constants__)); // Do a first evaluation to test size of return value
      
      label ncomp=toArray<T>::size(finalresult);
      
      label dim=parameters.size();//F.nDim();
      if (dim>1)
      {
        // Multidimensional: Use Cuba-Lib (does not work in one dimension)
        
#ifdef DEBUG
        Info<<"Multidimensional integrand"<<endl;
#endif
        
        // Step 1: integrate over given function
        int /*comp,*/ nregions, neval, fail;
        double integral[dim], error[dim], prob[dim];
  
        Cuhre
        //Divonne
            (
            dim, ncomp, 
            Integrand, NULL,
            EPSREL, EPSABS, 
#ifdef DEBUG
            VERBOSE | LAST, 
#else
            LAST, 
#endif
            MINEVAL, MAXEVAL, 
            
            KEY, //Cuhre
            
            //KEY1, KEY2, KEY3, MAXPASS, BORDER, MAXCHISQ, MINDEVIATION, //Divonne
              //NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
              
            &nregions, &neval, &fail, 
            integral, error, prob
            );
#ifdef DEBUG
        printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n",
        nregions, neval, fail);
        for( int comp = 0; comp < ncomp; ++comp )
            printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
            integral[comp], error[comp], prob[comp]);
#endif
        
        toArray<T>::convertBack(integral, finalresult);

        
      } else if (dim==1)
      {

        // One-dimensional: use Quadpack (GNU scientific library)
        
#ifdef DEBUG
        Info<<"Onedimensional integrand"<<endl;
#endif

       gsl_integration_workspace * w 
         = gsl_integration_workspace_alloc (10000);
       gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
 
       T result;
       double error;
       //double expected = 1.0;
       //double alpha = 1.0;
     
       gsl_function Func;
       Func.function = &gsl_integrand;
       //F.params = &alpha;
       int status=1;
       double relerr=1e-3, abserr=1e-6;
       while(status && (relerr<1))
       { 
        status=gsl_integration_qags (&Func, 0, 1, abserr, relerr, 10000,
                             w, &result, &error); 
        relerr*=1.1;
	abserr*=1.1;
        if (status)
        {
          Pout<<"result          = "<<result<<endl;
          Pout<<"estimated error = "<<error<<endl;
          Pout<<"Warning in integrator::integrate(): "
           <<"GSL integration failed, increasing tolerance to "<<relerr<<endl;
        }
       }

#ifdef DEBUG
       printf ("result          = % .18f\n", result);
       printf ("estimated error = % .18f\n", error);
#endif

       gsl_set_error_handler(old_handler);

       gsl_integration_workspace_free (w);

       finalresult=result;

      } else
      {
        // Zero-dimensional: may be after multiplication with dirac delta functions
#ifdef DEBUG
        Info<<"Dirac delta integrand"<<endl;
#endif
        
        //integrationParameter empty;
        finalresult=F.valueAt(constants);

#ifdef DEBUG
        printf ("result          = % .18f\n", finalresult);
#endif
        
      }
      
      return finalresult;

    }

    
    template<class T>
    T integrator<T>::integrate
    (
        const sumOfIntegrableFunctions<T>& F,
        const integrationParameter& parameters,
        const integrationParameter& constants
    )
    {
      T result=pTraits<T>::zero;
      /*
      for (label i=0; i<F.size(); i++)
      {
        integrationParameter param(F[i]().parameters(), 0.0);
        result+=integrate(F[i](), param, constants);
      }
      */
      for (typename sumOfIntegrableFunctions<T>::const_iterator i=F.begin(); i!=F.end(); i++)
      {
        integrationParameter param( (*i)().parameters(), 0.0);
        
        for (integrationParameter::const_iterator j=constants.begin();
             j!=constants.end(); j++)
          if (param.found(j.key())) param.erase(j.key());
        
        //Info<<param<<endl;
        result+=integrate( (*i)(), param, constants);
      }
      return result;
    }
    
}


// ************************************************************************* //
