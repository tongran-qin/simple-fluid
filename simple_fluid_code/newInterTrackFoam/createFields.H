
Info<< "\nReading field p\n" << endl << flush;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

/// RG-TQ Read Absolute pressure field
Info<< "\nReading field pAbs\n" << endl << flush;
volScalarField pAbs
(
    IOobject
    (
        "pAbs",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field U\n" << endl << flush;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Uliquid\n" << endl << flush;
volVectorField Uliquid
(
    IOobject
    (
        "Uliquid",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

/// RG-TQ Read Temperature field
Info<< "Reading field T\n" << endl << flush;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


/// RG-TQ Read concentration field of water vapor in vapor

Info<< "Reading field conVapor\n" << endl << flush;
volScalarField conVapor
(
    IOobject
    (
        "conVapor",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

/// RG-TQ Read concentration fields

Info<< "Reading field concentration\n" << endl << flush;
volScalarField con
(
    IOobject
    (
        "con",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#   include "createPhi.H"
phi.oldTime();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        //IOobject::MUST_READ,
        //IOobject::AUTO_WRITE
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimMass/dimVolume, 0)
);

volScalarField rho1
(
    IOobject
    (
        "rho1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField rhoVapor
(
    IOobject
    (
        "rhoVapor",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField rho1Vapor
(
    IOobject
    (
        "rho1Vapor",
        runTime.timeName(),
        mesh,
        //IOobject::MUST_READ,
        //IOobject::AUTO_WRITE
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimMass/dimVolume, 0)
);

volScalarField nc1
(
    IOobject
    (
        "nc1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
//    dimensionedScalar("0", dimless/dimVolume, 0)
);

volScalarField nt
(
    IOobject
    (
        "nt",
        runTime.timeName(),
        mesh,
        //IOobject::MUST_READ,
        //IOobject::AUTO_WRITE
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimless/dimVolume, 0)
);


///RG-TQ Add Temperature T
freeSurface interface(mesh, rho,rho1, nc1, U, p, T,con,phi);

volScalarField* fluidIndicatorPtr = NULL;

if(args.options().found("parallel"))
{
    fluidIndicatorPtr = new volScalarField
    (
        IOobject
        (
            "fluidIndicator",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
Info << "Reading fluidIndicator finished parallel" << endl;
}
else
{
    fluidIndicatorPtr = new volScalarField
    (
        IOobject
        (
            "fluidIndicator",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        interface.fluidIndicator()
    );
Info << "Reading fluidIndicator finished 1" << endl;
}

volScalarField& fluidIndicator = *fluidIndicatorPtr;
Info << "Reading fluidIndicator finished" << endl;
rho = fluidIndicator
   *(
        interface.rhoFluidA() 
      - interface.rhoFluidB()
    ) 
  + interface.rhoFluidB();
Info << "Reading rho finished" << endl;
/*
rho1 = fluidIndicator
   *(
        interface.rhoFluidA() 
      - interface.rhoFluidB1()
    ) 
  + interface.rhoFluidB1();

Uliquid = fluidIndicator * U;
rhoVapor = (1.0-fluidIndicator)*rho;
rho1Vapor = (1.0-fluidIndicator)*rho1;
*/
volScalarField kFluid = fluidIndicator
   *(
        interface.kFluidA() 
      - interface.c0().value()*interface.kFluidB1() - (1.0- interface.c0().value())*interface.kFluidB()
    ) 
  + interface.c0().value()*interface.kFluidB1() + (1.0- interface.c0().value())*interface.kFluidB();
Info << "Reading kA finished" << endl;
/*
dimensionedScalar TRefA = fvc::domainIntegrate(fluidIndicator*T)/fvc::domainIntegrate(fluidIndicator);
dimensionedScalar TRefB = fvc::domainIntegrate((1.0-fluidIndicator)*T)/fvc::domainIntegrate(1.0 - fluidIndicator);
Info << "Reading TrefA finished" << endl;
volVectorField buoForceA = interface.rhoFluidA()*interface.betaFluidA()*(T-TRefA)*interface.g();
volVectorField buoForceB = interface.rhoFluidB()*(T-TRefB)*interface.g()/T;
*/
Info<< "Reading field buoyancy\n" << endl << flush;
volVectorField buoForce
(
    IOobject
    (
        "buoForce",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//buoForce = fluidIndicator*buoForceA + (1.0 - fluidIndicator)*buoForceB;

//conVapor = (1.0-fluidIndicator) * con;

volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fluidIndicator
   *(
        interface.muFluidA() 
      -  interface.c0().value()*interface.muFluidB1() - (1.0- interface.c0().value())*interface.muFluidB()
    ) 
  + interface.c0().value()*interface.muFluidB1() + (1.0- interface.c0().value())*interface.muFluidB()
);

/*
///RG-TQ Read gas constant and psi
dimensionedScalar Rwater = interface.gasConstantWater();
dimensionedScalar Rair = interface.gasConstantAir();
//volScalarField psi = 1.0/(R*T);

///RG-TQ Initialize total mass, liquid mass and vapor mass
//dimensionedScalar M0 = fvc::domainIntegrate(rho);
dimensionedScalar Mv0 = fvc::domainIntegrate((1.0-fluidIndicator)*interface.rhoFluidB1());
dimensionedScalar Ml0 = fvc::domainIntegrate(fluidIndicator*interface.rhoFluidA());
dimensionedScalar Ma0 = fvc::domainIntegrate((1.0-fluidIndicator)*interface.rhoFluidB());

dimensionedScalar M0 = Mv0+Ml0+Ma0;

Info<< "Total mass: " << M0.value() << endl;
Info<< "Liquid mass: " << Ml0.value() << endl;
Info<< "Vapor mass: " << Mv0.value() << endl;
Info<< "Air mass: " << Ma0.value() << endl;
//Info<< min(con) << " < Initial concentration < " << max(con) << endl;
*/

/*
volScalarField betaFluidB = 1.0/T;

volScalarField betaVapor
(
    IOobject
    (
        "betaVapor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fluidIndicator
   *(
        interface.betaFluidA() 
      - betaFluidB
    ) 
  + betaFluidB
);
*/

///RG-TQ Read thermal diffusivity
volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fluidIndicator
   *(
        interface.kFluidA()/interface.rhoFluidA()/interface.CpFluidA() 
      - (interface.c0().value()*interface.kFluidB1() + (1.0- interface.c0().value())*interface.kFluidB())/interface.rhoFluidB()/(interface.c0().value()*interface.CpFluidB1()  + (1.0- interface.c0().value())*interface.CpFluidB())
    ) 
  + (interface.c0().value()*interface.kFluidB1() + (1.0- interface.c0().value())*interface.kFluidB())/interface.rhoFluidB()/(interface.c0().value()*interface.CpFluidB1()  + (1.0- interface.c0().value())*interface.CpFluidB())
//  + interface.kFluidB()/interface.rhoFluidB()/interface.CpFluidB()
);

///RG-TQ Read thermal diffusivity
volScalarField Df
(
    IOobject
    (
        "Df",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fluidIndicator
   *(
        interface.DfFluidA() 
      - interface.DfFluidB()
    ) 
   + interface.DfFluidB()
);

volScalarField Dfn
(
    IOobject
    (
        "Dfn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Df*nt
);


///RG-TQ Read specific heat capacity
volScalarField Cp
(
    IOobject
    (
        "Cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fluidIndicator
   *(
        interface.CpFluidA() 
      - interface.c0().value()*interface.CpFluidB1()  - (1.0- interface.c0().value())*interface.CpFluidB()
    ) 
  + interface.c0().value()*interface.CpFluidB1()  + (1.0- interface.c0().value())*interface.CpFluidB()
);
//scalarSquareMatrix A(3, 1.0);

///RG-TQ Define convergence criteria for pressure(dp) and temperature(dT)
dimensionedScalar dp("dp",dimless,0.0);
dimensionedScalar dT("dT",dimless,0.0);
dimensionedScalar dconVapor("dconVapor",dimless,0.0);
dimensionedScalar drho1("drho1",dimless,0.0);
dimensionedScalar dnc1("dnc1",dimless,0.0);
dimensionedScalar rhoSMALL("rhoSMALL",dimDensity,VSMALL);

#   include "createSf.H"

///RG-TQ Initialization of absolute pressure
//#   include "pAbs.H"

/// Read interpolators if present
if (interface.twoFluids())
{
    interface.interpolatorAB();
}

//Info<< "creatField1\n" << endl << flush;
#   include "setRefCell.H"
//Info<< "end of creatField\n" << endl << flush;
