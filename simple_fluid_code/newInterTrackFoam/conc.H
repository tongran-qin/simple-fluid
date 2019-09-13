///RG-TQ con.H

//dimensionedScalar Dfv = interface.DfFluidB();

/*
Df = fluidIndicator
   *(
        interface.DfFluidA() 
      - interface.DfFluidB()
    ) 
  + interface.DfFluidB();
*/

fvScalarMatrix conEqn
(
  fvm::ddt(con)
  + fvm::div(phi, con)
  - fvm::laplacian(Df, con)
);

conEqn.solve();

Info<<min(con)<<" < con < " << max(con) << endl;

con = fluidIndicator*(fvc::domainIntegrate((1.0-fluidIndicator)*(con))/fvc::domainIntegrate(1.0-fluidIndicator)) + (1.0-fluidIndicator)*(con);

nc1 = fluidIndicator*(fvc::domainIntegrate((1.0-fluidIndicator)*(nt*con))/fvc::domainIntegrate(1.0-fluidIndicator)) + (1.0-fluidIndicator)*(nt*con);

//partial pressure of component 1
volScalarField rho1 = nc1*interface.molarMass1();
volScalarField rhoD = (nt-nc1)*interface.molarMassD();
rho = (1.0-fluidIndicator)*(rho1 + rhoD)
	+fluidIndicator*interface.rhoFluidA()*(1.0 - interface.betaFluidA()*(T - interface.TRef()));

//volScalarField rho1 = con*pAbs/interface.gasConstant1()/T;
//volScalarField p1 = rho1*interface.gasConstant1()*T;
//rho = (1.0-fluidIndicator)*((pAbs*(1.0-con))/interface.gasConstantD()/T+rho1)
//	+fluidIndicator*interface.rhoFluidA()*(1.0 - interface.betaFluidA()*(T - interface.TRef()));
