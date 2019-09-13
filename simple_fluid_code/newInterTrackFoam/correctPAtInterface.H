// if (!ddtPhiCorr)
{

    // do wall boundary condition by hand

    Info << "Correcting pressure gradient at wall" << endl;
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

//		buoyancy force changed signs
//            pPatch.gradient() =-
            pPatch.gradient() =
                ((
                        mesh.boundary()[patchI].Sf()
                      & buoForce.boundaryField()[patchI]

                )
               /mesh.magSf().boundaryField()[patchI]);
        }
    }


    Info << "Correcting pressure gradient at the interface" << endl;

        if
        (
            p.boundaryField()[interface.bPatchID()].type()
         == fixedGradientFvPatchField<scalar>::typeName
        )
        {
            fixedGradientFvPatchField<scalar>& pB =
                refCast<fixedGradientFvPatchField<scalar> >
                (
                    p.boundaryField()[interface.bPatchID()]
                );

    	pB.gradient() += 
	mesh.boundary()[interface.bPatchID()].nf() 
	& buoForce.boundaryField()[interface.bPatchID()];

        }

/*
    Info << "Correcting interface flux" << endl;

    phi.boundaryField()[interface.aPatchID()] +=
        (
            (
                U.boundaryField()[interface.aPatchID()].snGrad()
               /mesh.boundary()[interface.aPatchID()].deltaCoeffs()
            )
          & mesh.Sf().boundaryField()[interface.aPatchID()]
        );
*/

}

