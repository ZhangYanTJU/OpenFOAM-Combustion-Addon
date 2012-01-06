
    scalar nm1Coeff = averagingTime / ( averagingTime + runTime.deltaT().value() );

    averagingTime +=runTime.deltaT().value();

    scalar nCoeff = runTime.deltaT().value() / averagingTime;

    R += sqr(Umean);
    Umean = nm1Coeff*Umean + nCoeff*U;
    R = nm1Coeff*R + nCoeff*sqr(U) - sqr(Umean);
    Bmean = nm1Coeff*Bmean + nCoeff*turbulence->B();
    epsilonMean = nm1Coeff*epsilonMean + nCoeff*turbulence->epsilon();

    pPrime2Mean += sqr(pMean);
    pMean = nm1Coeff*pMean + nCoeff*p;
    pPrime2Mean = nm1Coeff*pPrime2Mean + nCoeff*sqr(p) - sqr(pMean);

    Tmean = nm1Coeff*Tmean + nCoeff*T;
