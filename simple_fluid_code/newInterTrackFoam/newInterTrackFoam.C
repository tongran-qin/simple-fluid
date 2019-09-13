/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Application
    interfaceTrackinFoam

Description
    Incompressible laminar CFD code for simulation of a single bubble rising 
    in a stil liquid. Interface between fluid phases is tracked using moving
    mesh.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "freeSurface.H"

#include "wallFvPatch.H"
#include "fixedGradientFvPatchFields.H"
#	include "scalarMatrices.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"

#   include "createDynamicFvMesh.H"

#   include "createFields.H"

#   include "initContinuityErrs.H"
//#	include "scalarSquareMatrix.H"
//#	include "scalarMatrices.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;

    ///RG-TQ Initial volume of liqiud and vapor
    const dimensionedScalar& VlInitial = fvc::domainIntegrate(fluidIndicator);
    const dimensionedScalar& VvInitial = fvc::domainIntegrate((1.0-fluidIndicator));
    Info<< "Initial Liquid Volume ----> " << VlInitial.value() << endl;
    Info<< "Initial Vapor Volume ---->" << VvInitial.value() << endl;

    if (runTime.value() < 1e-15)
    {
    interface.initializeConcentration();
    rho = fluidIndicator
   *(
        interface.rhoFluidA() 
      - interface.rhoFluidB()
    ) 
  + interface.rhoFluidB();

    }
    interface.updateProperties();
//    interface.initializeMass();

	dimensionedScalar MlOld = fvc::domainIntegrate(fluidIndicator*rho);
//	dimensionedScalar MvOld = interface.MB1();
	dimensionedScalar MOld = MlOld + interface.MB1();

	interface.moveWholeFreeSurface(MOld, interface.Mair0());
  	dimensionedScalar L0("L0", dimLength, 1.0);
  	dimensionedVector R0 = L0*interface.geometricMeanPosition();
	pAbs = p + rho*(interface.g() & mesh.C()) - interface.rhoFluidA()*(interface.g() & R0)*fluidIndicator;
	dimensionedScalar pCorr = interface.updatePressureDensity2();
  	pAbs += pCorr;

//#   include "pAbs.H"
//#   include "pAbs2.H"

//    nt = pAbs/interface.gasConstant1()/interface.molarMass1()/T;
	dimensionedScalar pdim("pdim",dimPressure,1.0);
    	nt = interface.pTotal().value()*pdim/interface.gasConstant1()/interface.molarMass1()/((interface.TLeft() + interface.TRight())/2.0);
    	nt = fluidIndicator*(fvc::domainIntegrate((1.0-fluidIndicator)*(nt))/fvc::domainIntegrate(1.0-fluidIndicator)) + (1.0-fluidIndicator)*(nt);

    Info<< "initial concentration in vapor [nt] ----> " << min(nt) <<max(nt) <<endl;

    Dfn = Df*nt;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info << "Time = " << runTime.value() << endl << endl;

#       include "readPISOControls.H"

#       include "readFreeSurfaceControls.H"

        //interface.updateProperties();

//	dimensionedScalar VvOld = fvc::domainIntegrate(fluidIndicator);
//	dimensionedScalar MvOld = fvc::domainIntegrate((1.0-fluidIndicator)*pAbs/(interface.gasConstant1()*T));

        interface.moveMeshPointsForOldFreeSurfDisplacement();

        interface.updateDisplacementDirections();

        interface.predictPoints();

//	    dimensionedScalar Mnew = fvc::domainIntegrate(fluidIndicator*rho) + fvc::domainIntegrate((1.0-fluidIndicator)*con*pAbs/(interface.gasConstant1()*T));
//	    interface.moveWholeFreeSurface(Mnew,MOld);
//#           include "pAbs.H"
//#           include "pAbs2.H"

	interface.moveWholeFreeSurface(MOld, interface.Mair0());
//  	L0("L0", dimLength, 1.0);
  	R0 = L0*interface.geometricMeanPosition();
	pAbs = p + rho*(interface.g() & mesh.C()) - interface.rhoFluidA()*(interface.g() & R0)*fluidIndicator;
	pCorr = interface.updatePressureDensity2();
  	pAbs += pCorr;

        Info<< "\nMax surface Courant Number = "
            << interface.maxCourantNumber() << endl << endl;

//0820 update all the boundary condition once first
            interface.updateBoundaryConditions();

        ///RG-TQ setting the criteria of the loop by dp/dT
        dT=1.0;
        drho1=1.0;
	dnc1=1.0;
        dp=1.0;
        int corr=0;
        while(dp.value()>1e-4||dT.value()>1e-4||dnc1.value()>1e-4)
//        while(dp.value()>1e-4||dT.value()>1e-4||drho1.value()>1e-4)
        {  
        if(corr<1e3)
        {
            /// RG-TQ  
            volScalarField pOld = p;

            // Update interface bc
//	0820
//            interface.updateBoundaryConditions();
            interface.updateVelocity();
            interface.updatePressure();
#               include "correctPAtInterface.H"

            // Make the fluxes relative
            phi -= fvc::meshPhi(rho, U);

#           include "CourantNo.H"

            if(dp.value()<1e-4)
            {
            	volScalarField TOld = T;
//            	volScalarField rho1VaporOld = rho1Vapor;
//            	volScalarField nc1Old = nc1;
            	volScalarField conOld = con;
//0820
            interface.updateMassFlux();
            interface.updateTemperature();

            	///RG-TQ Temperature
#           	include "TEqn.H"
#           	include "updateTBoundary.H"

            	dT = (max(T-TOld)-min(T-TOld))/(max(T)-min(T));

//            if(dT.value()<1e-4)
//            {
            	///RG-TQ water vppor concentration in vapor phase
//#           	include "DensityandConcentration.H"
//	    	dnc1 = (max(nc1-nc1Old) - min(nc1-nc1Old))/(max(nc1)- min(nc1));
//        	dnc1 = max(mag(nc1*(1.0-fluidIndicator)-nc1Old*(1.0-fluidIndicator)))/(min(nc1));

#           	include "conc.H"
	    	dnc1 = (max(con-conOld) - min(con-conOld))/(max(con)- min(con));

//	    }

//        	dnc1 = max(mag(nc1*(1.0-fluidIndicator)-nc1Old*(1.0-fluidIndicator)))/(min(nc1));

//	    	drho1 = (max(rho1-rho1Old) - min(rho1-rho1Old))*(1.0-fluidIndicator)/(max(rho1Vapor)- min(rho1Vapor));
//        	drho1 = max(mag(rho1*(1.0-fluidIndicator)-rho1Old*(1.0-fluidIndicator)))/(min(rho1)+rhoSMALL);
//        	drho1 = max(mag(rho1Vapor-rho1VaporOld))/(min(rho1Vapor)+rhoSMALL);

		///buoyancy force is updated, then rho is spatially averaged
#           	include "buoyancy.H"
            }

    		volScalarField rhoRef = fluidIndicator
   			*(fvc::domainIntegrate(fluidIndicator*rho)/fvc::domainIntegrate(fluidIndicator))
			+ (1.0-fluidIndicator)
			*(fvc::domainIntegrate((1.0-fluidIndicator)*rho)/fvc::domainIntegrate(1.0 - fluidIndicator)); 

//set rho as spatial avergaged value in two phases
            rho = rhoRef;
//Matrix <scalarField,scalar> A (3,1)
//#include "scalarMatrices.H"
//#include "tetFemMatrices.H"
//scalarSquareMatrix A(3, 1.0);
//volTensorField gradU = fvc::grad(U);
            fvVectorMatrix UEqn
            (
                fvm::ddt(rho, U)
              + fvm::div(fvc::interpolate(rho)*phi, U, "div(phi,U)")
              - fvm::laplacian(mu, U)
	      - buoForce	
//              + rho*betav*(T - interface.TRef())*interface.g()             
            );

            solve(UEqn == -fvc::grad(p));

            // --- PISO loop
            for (int i=0; i<nCorr; i++)
            {
/*
                volScalarField UA = UEqn.A();
                volVectorField UH = UEqn.H();

                U = UH/UA;

#               include "calcPhi.H"

#               include "scalePhi.H"
*/
                volScalarField UA = UEqn.A();

                U = UEqn.H()/UA;

                phi = (fvc::interpolate(U) & mesh.Sf());

#               include "scalePhi.H"

#               include "correctPhiAtInterface.H"
	
#           	include "correctBCsAtWalls.H"

                // Non-orthogonal pressure corrector loop
                for (label nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                {
                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(1.0/fvc::interpolate(UA), p)
                     == fvc::div(phi)
                    );

#                   include "setReference.H"

                    pEqn.solve();

                    if (nonOrth == nNonOrthCorr)
                    {
                        phi -= pEqn.flux();
                    }
                }

#               include "continuityErrs.H"

                // Momentum corrector
                U -= fvc::grad(p)/UA;
                U.correctBoundaryConditions();

            }

            interface.correctPoints();

#           include "freeSurfaceContinuityErrs.H"

#           include "updateSf.H"

            ///RG-TQ Update dp and dT
            dp = (max(p-pOld)-min(p-pOld))/(max(p)-min(p));

            Info << "dp = " << dp.value() << ", dT = " << dT.value()<< ", dnc1 = " << dnc1.value() << endl;

            ///RG-TQ Update absolute pressure
//#           include "pAbs.H"
//#           include "pAbs2.H"

            corr+=1;
        }

        else
        {
            FatalErrorIn("newInterTrackFoam.C")
                << "Iterations for p,T and concentration haven't converged after 1e3 steps, probably won't converge"
                    << abort(FatalError);
        }

        }

#       include "volContinuity.H"


//	    Mnew = fvc::domainIntegrate(fluidIndicator*rho) + fvc::domainIntegrate((1.0-fluidIndicator)*con*pAbs/(interface.gasConstant1()*T));
//	    interface.moveWholeFreeSurface(Mnew,MOld);
//#           include "pAbs.H"
//#           include "pAbs2.H"

	interface.moveWholeFreeSurface(MOld, interface.Mair0());
//  	L0("L0", dimLength, 1.0);
  	R0 = L0*interface.geometricMeanPosition();
	pAbs = p + rho*(interface.g() & mesh.C()) - interface.rhoFluidA()*(interface.g() & R0)*fluidIndicator;
	pCorr = interface.updatePressureDensity2();
  	pAbs += pCorr;

	nt = (1.0-fluidIndicator)*(pAbs/interface.gasConstant1()/interface.molarMass1()/T)+fluidIndicator*nt;
	dimensionedScalar ntAvg = fvc::domainIntegrate((1.0-fluidIndicator)*nt)/fvc::domainIntegrate((1.0-fluidIndicator));
        Info << "ntAvg = " << ntAvg.value() << endl;
	nt = (1.0-fluidIndicator)*ntAvg + fluidIndicator*nt;

	    dimensionedScalar Mnew = fvc::domainIntegrate(fluidIndicator*rho) + fvc::domainIntegrate((1.0-fluidIndicator)*con*pAbs/(interface.gasConstant1()*T));
	    dimensionedScalar MaNew = fvc::domainIntegrate((1.0-fluidIndicator)*(1.0-con)*pAbs/(interface.gasConstantD()*T));
	    Info<<"Mass decrease in liquid: " << (fvc::domainIntegrate(fluidIndicator*rho) - MlOld).value() <<endl;
	    Info<<"Mass increase in vapor: " << (fvc::domainIntegrate((1.0-fluidIndicator)*con*pAbs/(interface.gasConstant1()*T))-interface.MB1()).value() <<endl;

	    Info<<"Relative Mass increase in air: " << (MaNew.value()-interface.Mair0().value())/interface.Mair0().value() <<endl;
/*
	    dimensionedScalar MvNew = fvc::domainIntegrate((1.0-fluidIndicator)*con*pAbs/(interface.gasConstant1()*T));
	    dimensionedScalar MaNew = fvc::domainIntegrate((1.0-fluidIndicator)*(1.0-con)*pAbs/(interface.gasConstantD()*T));
	    dimensionedScalar Vvnew1 = fvc::domainIntegrate(fluidIndicator);

	    Info<<"Mass increase in vapor: " << (MvNew.value()-interface.MB1().value()) <<endl;
	    Info<<"Relative Mass increase in vapor: " << (MvNew.value()-interface.MB1().value())/interface.MB1().value() <<endl;
	    Info<<"Relative Mass increase in air: " << (MaNew.value()-interface.Mair0().value())/interface.Mair0().value() <<endl;
*/
        ///RG-TQ Moving whole plane for single phase volume conservation    
//        dimensionedScalar VvCorr = fvc::domainIntegrate(1.0-fluidIndicator);
//        interface.updateVolume(VvCorr,VvInitial);
//        interface.moveWholeFreeSurface();
/*
        ///RG-TQ Calculate relative change in volume of liquid and vapor
        dimensionedScalar Vl = fvc::domainIntegrate(fluidIndicator);
        dimensionedScalar Vv = VvCorr;//fvc::domainIntegrate((1.0-fluidIndicator));
        Info<< "Relative Liquid Volume change ---->" << (Vl.value()-VlInitial.value())/VlInitial.value() << endl;
        Info<< "Relative Vapor  Volume change ---->" << (Vv.value()-VvInitial.value())/VvInitial.value() << endl;
*/

        ///RG-TQ output dp and dT
        Info<<"Iteration# = " << corr <<"---------> dp = " << dp.value() << endl;
        Info<< "---------> dT = " << dT.value() << endl;
        Info<< "---------> dnc1 = " << dnc1.value() << endl;
//        Info<< "---------> drho1 = " << drho1.value() << endl;

        ///RG-TQ Uliquid
        Uliquid = fluidIndicator * U;

	conVapor = fluidIndicator*(fvc::domainIntegrate((1.0-fluidIndicator)*con)/fvc::domainIntegrate(1.0-fluidIndicator)) + (1.0-fluidIndicator)*con;

	rhoVapor = fluidIndicator*(fvc::domainIntegrate((1.0-fluidIndicator)*rho)/fvc::domainIntegrate(1.0-fluidIndicator)) + (1.0-fluidIndicator)*rho;

//rho1 in the vapor, for visualization purposes only
	rho1Vapor = fluidIndicator*(fvc::domainIntegrate((1.0-fluidIndicator)*rho1)/fvc::domainIntegrate(1.0-fluidIndicator)) + (1.0-fluidIndicator)*rho1;

        Info << "Total surface tension force: "
            << interface.totalSurfaceTensionForce() << endl;

        vector totalForce =
            interface.totalViscousForce() 
          + interface.totalPressureForce();

        Info << "Total force: " << totalForce << endl;

        runTime.write();

        Info << "ExecutionTime = "
            << scalar(runTime.elapsedCpuTime())
            << " s\n" << endl << endl;
    }

    Info << "End\n" << endl;

    return(0);
}

// ************************************************************************* //
