///RG-TQ con.H

//Diffusion coefficient
/*
Df = fluidIndicator
   *(
        interface.DfFluidA() 
      - interface.DfFluidB()
    ) 
  + interface.DfFluidB();
*/
//advection-diffusion equation for rho1

//nt = pAbs/interface.gasConstant1()/interface.molarMass1()/T;
//nc1 = con*nt;
volScalarField Dfc = Df*con;
//volScalarField jdiff =  fvc::laplacian(Dfn, con);

fvScalarMatrix nc1Eqn
(
  fvm::ddt(nc1)
  + fvm::div(phi, nc1)
//  - fvc::laplacian(Dfn, con)
//  - fvm::laplacian(Df*(nc1/nt), con)
  - fvm::laplacian(Df, nc1)
  + fvc::laplacian(Dfc, nt)

);

nc1Eqn.solve();

//temperature dependence is considered for component 1, in the vapor phase
//volScalarField rho1_T = rho1 *(1.0 - (T - interface.TRef())/T );

//partial pressure of component 1
volScalarField rho1 = nc1*interface.molarMass1();
volScalarField p1 = rho1*interface.gasConstant1()*T;

//total density in vapor is calculated and temperature dependence in liquid density is considered
//update the density field
rho = (1.0-fluidIndicator)*((pAbs-p1)/interface.gasConstantD()/T+rho1)
	+fluidIndicator*interface.rhoFluidA()*(1.0 - interface.betaFluidA()*(T - interface.TRef()));

//volScalarField rho_T = (1.0-fluidIndicator)*((pAbs-rho1*interface.gasConstant1()*T)/interface.gasConstantD()/T+rho1)
//	+fluidIndicator*interface.rhoFluidA()*(1.0 - interface.betaFluidA()*(T - interface.TRef()));

//concentration field of component 1 is updated
con = fluidIndicator*(fvc::domainIntegrate((1.0-fluidIndicator)*(nc1/nt))/fvc::domainIntegrate(1.0-fluidIndicator)) + (1.0-fluidIndicator)*(nc1/nt);

//con = fluidIndicator*(fvc::domainIntegrate((1.0-fluidIndicator)*(p1/pAbs))/fvc::domainIntegrate(1.0-fluidIndicator)) + (1.0-fluidIndicator)*(p1/pAbs);

nc1 = fluidIndicator*(fvc::domainIntegrate((1.0-fluidIndicator)*(nc1))/fvc::domainIntegrate(1.0-fluidIndicator)) + (1.0-fluidIndicator)*(nc1);

        Info<<min(con)<<" < con < " << max(con) << endl;

//rho1 in the vapor, for visualization purposes only
//rho1Vapor = fluidIndicator*(fvc::domainIntegrate((1.0-fluidIndicator)*rho1)/fvc::domainIntegrate(1.0-fluidIndicator)) + (1.0-fluidIndicator)*rho1;

//******************************************************************
/*
    // reference density, average value in the two phases

    volScalarField rhoRef = fluidIndicator
   			*(fvc::domainIntegrate(fluidIndicator*rho_T)/fvc::domainIntegrate(fluidIndicator))
			+ (1.0-fluidIndicator)
			*(fvc::domainIntegrate((1.0-fluidIndicator)*rho_T)/fvc::domainIntegrate(1.0 - fluidIndicator)); 

    // buoyancy force
    buoForce = ( rho_T - rhoRef ) * interface.g();

//set rho as spatial avergaged value in two phases
//    rho = rhoRef;
*/
//******************************************************************
//volScalarField pD = (1.0-fluidIndicator)*(pAbs-p1)+fluidIndicator*pAbs;
