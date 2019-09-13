///RG-TQ pAbs.H

{
  dimensionedScalar L0("L0", dimLength, 1.0);
  dimensionedVector R0 = L0*interface.geometricMeanPosition();
  pAbs = p + rho*(interface.g() & mesh.C())
       - interface.rhoFluidA()*(interface.g() & R0)*fluidIndicator;
  
    dimensionedScalar Vv = fvc::domainIntegrate((1.0-fluidIndicator));

    Info<< "Vapor Volume: " << Vv.value() << endl;

//  dimensionedScalar Mair = fvc::domainIntegrate((1.0-fluidIndicator)*p*(1.0-con)/(interface.gasConstantD()*T));
/*
  dimensionedScalar pCorr = (interface.Mair0()-Mair)/
//  dimensionedScalar pCorr = (interface.Mair0())/
  (
  fvc::domainIntegrate((1.0-con)*(1.0-fluidIndicator)/(interface.gasConstantD()*T))
  );
*/

  dimensionedScalar Mair = fvc::domainIntegrate((1.0-fluidIndicator)*(pAbs-nt*con*interface.molarMass1()*interface.gasConstant1()*T)/(interface.gasConstantD()*T));

  Info << "Mass of air calculated from dynamic pressure: "<< Mair.value() << endl;

  dimensionedScalar pCorr = (interface.Mair0()-Mair)/
//  dimensionedScalar pCorr = (interface.Mair0())/
  (
  fvc::domainIntegrate((1.0-fluidIndicator)/(interface.gasConstantD()*T))
  );


//  pCorr = 4220.6*(pCorr/pCorr.value());//5037.2//8101.58//4555.94

  pAbs += pCorr;

  interface.updatePressureDensity(pCorr);

  Info<< "Pressure offset " << pCorr.value() << endl;


//  Info<< "Liquid mass " << Ml.value() << endl;
//  Info<< "rhoAir " << rhoAir.value() << endl;
//  Info<< "rhoVapor " << rhoVapor.value() << endl;
//  Info<< "Vapor+Air mass " << M0.value()-Ml.value() << endl;
//  Info<< "Pressure offset " << pCorr.value() << endl;
//  Info<< "Vapor density " << interface.rhoFluidB().value() << endl;

}
