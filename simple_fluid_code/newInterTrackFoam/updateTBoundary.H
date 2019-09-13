///RG-TQ TEqn.H

forAll (T.boundaryField(), patchI)
{
    if
    (
	(
	    mesh.boundary()[patchI].name()
             == "leftWall"
	)
	&&
	(
	    T.boundaryField()[patchI].type()
	 == fixedValueFvPatchScalarField::typeName
	)

    )
    {
//	vectorField Cf = mesh.boundary()[patchI].Cf();
//	Info << "CFLeft: " <<  Cf << endl;

	scalarField TFluid = T.boundaryField()[patchI].patchInternalField();
//	Info << "TFluid:" << TFluid << endl;
	scalarField DnF = mesh.boundary()[patchI].deltaCoeffs();
//	Info << "DnF:" << DnF << endl;

	scalarField sigma = kFluid.boundaryField()[patchI]*DnF+interface.kWall().value()/interface.thWall().value();
//	Info << "sigma:" << sigma << endl;
//	Info << "kFluid.boundaryField()[patchI]:" << kFluid.boundaryField()[patchI] << endl;

	scalarField TiLeft = ((interface.kWall().value()*interface.TLeft().value()/interface.thWall().value())+kFluid.boundaryField()[patchI]*TFluid*DnF)/sigma;

//	Info << "TLeftEff on the left wall:" << TiLeft << endl;

	T.boundaryField()[patchI] == TiLeft;

    }


    if
    (
	(
	    mesh.boundary()[patchI].name()
             == "rightWall"
	)
	&&
	(
	    T.boundaryField()[patchI].type()
	 == fixedValueFvPatchScalarField::typeName
	)

    )
    {
	scalarField TFluid = T.boundaryField()[patchI].patchInternalField();
//        vectorField nF = mesh.boundary()[patchI].nf();
//area center
//	    vectorField Cf = mesh.boundary()[patchI].Cf();
//neighbour cell center
//	    vectorField Cn = mesh.boundary()[patchI].Cn();
//	    vectorField delta = mesh.boundary()[patchI].delta();

//	Info << "TFluid:" << TFluid << endl;

//	Info << "nF: " <<  nF << endl;
//	Info << "CFRight: " <<  Cf << endl;
//	Info << "Cn: " <<  Cn << endl;
//	Info << "delta: " <<  delta << endl;
//	Info << "delta2: " <<  Cf-Cn << endl;

	scalarField DnF = mesh.boundary()[patchI].deltaCoeffs();
//	Info << "DnF:" << DnF << endl;
//	DnF = 1.0/mag(delta);
//	Info << "DnF2:" << DnF << endl;
//	DnF = delta&nF;
//	Info << "DnFnew:" << DnF << endl;
//	Info << "DnFnew:" << (1.0/DnF) << endl;

	scalarField sigma = kFluid.boundaryField()[patchI]*DnF+interface.kWall().value()/interface.thWall().value();
//	Info << "sigma:" << sigma << endl;
//	Info << "kFluid.boundaryField()[patchI]:" << kFluid.boundaryField()[patchI] << endl;

	scalarField TiRight = ((interface.kWall().value()*interface.TRight().value()/interface.thWall().value())+kFluid.boundaryField()[patchI]*TFluid*DnF)/sigma;

//	Info << "TRightEff on the right wall:" << TiRight << endl;

	T.boundaryField()[patchI] == TiRight;

    }

}




