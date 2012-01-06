
IOdictionary thickenedFlameProperties
(
    IOobject
    (
        "thickenedFlameProperties",
        runTime.constant(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

scalar TF = readScalar
    (
        thickenedFlameProperties.lookup("thickeningFactor")
    );
