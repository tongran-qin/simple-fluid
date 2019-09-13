///RG-TQ TEqn.H
/*
kappa = fluidIndicator
   *(
        interface.kFluidA()/interface.rhoFluidA()/interface.CpFluidA() 
      - interface.kFluidB()/interface.rhoFluidB()/interface.CpFluidB()
    ) 
  + interface.kFluidB()/interface.rhoFluidB()/interface.CpFluidB();
*/

fvScalarMatrix TEqn
(
  fvm::ddt(T)
  + fvm::div(phi, T)
  - fvm::laplacian(kappa, T)
);

TEqn.solve();
/*
	TRefA = fvc::domainIntegrate(fluidIndicator*T)/fvc::domainIntegrate(fluidIndicator);
	TRefB = fvc::domainIntegrate((1.0-fluidIndicator)*T)/fvc::domainIntegrate(1.0 - fluidIndicator);

	buoForceA = interface.rhoFluidA()*interface.betaFluidA()*(T-TRefA)*interface.g();
	buoForceB = interface.rhoFluidB()*(T-TRefB)*interface.g()/T;

	buoForce = fluidIndicator*buoForceA + (1.0 - fluidIndicator)*buoForceB;
*/
