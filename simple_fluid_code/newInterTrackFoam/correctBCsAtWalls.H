    forAll (p.boundaryField(), patchI)
    {
        if
        (
            (
                mesh.boundary()[patchI].type()
             == wallFvPatch::typeName
            )
         &&
            (
                p.boundaryField()[patchI].type()
             == fixedGradientFvPatchScalarField::typeName
            )
        )
        {

        Info << "Correct phi at wall: "
	     << mesh.boundary()[patchI].name() << endl;

        phi.boundaryField()[patchI] +=
            p.boundaryField()[patchI].snGrad()
           *mesh.magSf().boundaryField()[patchI]
	   /UA.boundaryField()[patchI];

        }
    }

//    TRefA = fvc::domainIntegrate(fluidIndicator*T)/fvc::domainIntegrate(fluidIndicator);
//    TRefB = fvc::domainIntegrate((1.0-fluidIndicator)*T)/fvc::domainIntegrate(1.0 - fluidIndicator);
//    buoForceA = interface.rhoFluidA()*interface.betaFluidA()*(T-TRefA)*interface.g();
//    buoForceB = interface.rhoFluidB()*(T-TRefB)*interface.g()/T;
//    buoForce = fluidIndicator*buoForceA + (1.0 - fluidIndicator)*buoForceB;
