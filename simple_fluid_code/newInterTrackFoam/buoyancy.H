
    // reference density, average value in the two phases

    volScalarField rhoRef = fluidIndicator
   			*(fvc::domainIntegrate(fluidIndicator*rho)/fvc::domainIntegrate(fluidIndicator))
			+ (1.0-fluidIndicator)
			*(fvc::domainIntegrate((1.0-fluidIndicator)*rho)/fvc::domainIntegrate(1.0 - fluidIndicator)); 

    // buoyancy force
    buoForce = ( rho - rhoRef )* interface.g() ;
/*
volTensorField gradU = fvc::grad(U);
scalarSquareMatrix A(3, 1.0);
fvVectorMatrix MA (fvm::laplacian(mu, U));
fvVectorMatrix MM = MA&MA;
*/

//set rho as spatial avergaged value in two phases
    rho = rhoRef;

/*
    buoForceA = interface.rhoFluidA()*interface.betaFluidA()*(T-interface.TRef())*interface.g();
    buoForceB = interface.rhoFluidB()*(T-interface.TRef())*interface.g()/T;
    buoForce = fluidIndicator*buoForceA + (1.0 - fluidIndicator)*buoForceB;
*/
//******************************************************************

/*
    // do wall boundary condition by hand
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

            fixedGradientFvPatchScalarField& pPatch =
                refCast<fixedGradientFvPatchScalarField>
                (
                    p.boundaryField()[patchI]
                );

            pPatch.gradient() =-
                (
                        mesh.boundary()[patchI].Sf()
                      & buoForce.boundaryField()[patchI]

                )
               /mesh.magSf().boundaryField()[patchI];

        Info << "Correct phi at wall: "
	     << mesh.boundary()[patchI].name() << endl;

        phi.boundaryField()[patchI] +=
            p.boundaryField()[patchI].snGrad()
           *mesh.magSf().boundaryField()[patchI]
	   /UA.boundaryField()[patchI];

        }
    }
*/
//    TRefA = fvc::domainIntegrate(fluidIndicator*T)/fvc::domainIntegrate(fluidIndicator);
//    TRefB = fvc::domainIntegrate((1.0-fluidIndicator)*T)/fvc::domainIntegrate(1.0 - fluidIndicator);
//    buoForceA = interface.rhoFluidA()*interface.betaFluidA()*(T-TRefA)*interface.g();
//    buoForceB = interface.rhoFluidB()*(T-TRefB)*interface.g()/T;
//    buoForce = fluidIndicator*buoForceA + (1.0 - fluidIndicator)*buoForceB;
