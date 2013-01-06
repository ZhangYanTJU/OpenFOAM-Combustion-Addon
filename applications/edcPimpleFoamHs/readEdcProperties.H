{
  IOdictionary edcProperties
    (
     IOobject
     (
      "edcProperties",
      runTime.constant(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
     )
    );

  //beta=readScalar(edcProperties.lookup("beta"));
  Cmag=readScalar(edcProperties.lookup("Cmag"));
  s=readScalar(edcProperties.lookup("stoichAirToFuelMassRatio"));

  fuel=word(edcProperties.lookup("fuel"));
  ox=word(edcProperties.lookup("ox"));
  prod=HashTable<scalar, word>(edcProperties.lookup("prod"));}
