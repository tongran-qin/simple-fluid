// if (!ddtPhiCorr)
{
/*
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
*/

    Info << "Correcting interface flux" << endl;

    phi.boundaryField()[interface.aPatchID()] +=
        (
            (
                U.boundaryField()[interface.aPatchID()].snGrad()
               /mesh.boundary()[interface.aPatchID()].deltaCoeffs()
            )
          & mesh.Sf().boundaryField()[interface.aPatchID()]
        );


}

