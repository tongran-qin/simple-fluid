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

Description

\*---------------------------------------------------------------------------*/

#include "freeSurface.H"
#include "volFields.H"
#include "primitivePatchInterpolation.H"
#include "fixedGradientFvPatchFields.H"
#include "zeroGradientCorrectedFvPatchFields.H"
#include "fixedGradientCorrectedFvPatchFields.H"
#include "fixedValueCorrectedFvPatchFields.H"
#include "emptyFaPatch.H"
#include "wedgeFaPatch.H"
//#include "demandDrivenData.H"
#include "EulerDdtScheme.H"
#include "CrankNicholsonDdtScheme.H"
#include "backwardDdtScheme.H"
#include "facLnGrad.H"

#include "wallFvPatch.H"

#include "fixedGradientFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "symmetryFvPatchFields.H"
#include "transformField.H"

#include "fixedValueFaPatchFields.H"

#include "tetFemMatrices.H"
#include "tetPointFields.H"
#include "faceTetPolyPatch.H"
#include "tetPolyPatchInterpolation.H"
#include "fixedValueTetPolyPatchFields.H"

#include "twoDPointCorrector.H"

#include "fixedValuePointPatchFields.H"
//z0412s
#include "coordinateSystem.H"
#include "scalarMatrices.H"
#include "fixedGradientFaPatchFields.H"

//#include "demandDrivenData.H"
//#include "facLnGrad.H"
//#include "fixedValueFaPatchFields.H"
//z0412e

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(freeSurface, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void freeSurface::clearOut()
{
    deleteDemandDrivenData(interpolatorABPtr_);
    deleteDemandDrivenData(interpolatorBAPtr_);
    deleteDemandDrivenData(controlPointsPtr_);

    deleteDemandDrivenData(TiPtr_);
    deleteDemandDrivenData(TiDiffPtr_);
    deleteDemandDrivenData(dTiPtr_);   
    deleteDemandDrivenData(gradTiPtr_);
    deleteDemandDrivenData(massFluxPtr_);
    deleteDemandDrivenData(heatFluxlPtr_);    
    deleteDemandDrivenData(heatFluxvPtr_);
    deleteDemandDrivenData(UiPtr_);
    deleteDemandDrivenData(conciPtr_);    

    deleteDemandDrivenData(motionPointsMaskPtr_);
    deleteDemandDrivenData(pointsDisplacementDirPtr_);
    deleteDemandDrivenData(facesDisplacementDirPtr_);
    deleteDemandDrivenData(totalDisplacementPtr_);
    deleteDemandDrivenData(aMeshPtr_);
    deleteDemandDrivenData(UsPtr_);
    deleteDemandDrivenData(phisPtr_);
    deleteDemandDrivenData(phiPtr_);
    deleteDemandDrivenData(ddtPhiPtr_);
    deleteDemandDrivenData(surfactConcPtr_);
    deleteDemandDrivenData(surfaceTensionPtr_);
    deleteDemandDrivenData(surfactantPtr_);
    deleteDemandDrivenData(fluidIndicatorPtr_);
    deleteDemandDrivenData(contactAnglePtr_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

freeSurface::freeSurface
(
    dynamicFvMesh& m,
    const volScalarField& rho,
    volScalarField& rho1b,
    volScalarField& ntb,
    volVectorField& Ub, 
    volScalarField& Pb, 
    ///RG-TQ
    volScalarField& Tb,
//    volScalarField& conLiquidb,
    volScalarField& conb,
    const surfaceScalarField& sfPhi
)
:
    IOdictionary
    (
        IOobject
        (
            "freeSurfaceProperties",
            Ub.mesh().time().constant(),
            Ub.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(m),
    rho_(rho),
    rho1_(rho1b),
    nt_(ntb),
    U_(Ub),
    p_(Pb),
    ///RG-TQ
    T_(Tb),
//    conLiquid_(conLiquidb),
    con_(conb),
    phi_(sfPhi),
    curTimeIndex_(Ub.mesh().time().timeIndex()),
    twoFluids_
    (
        this->lookup("twoFluids")
    ),
    normalMotionDir_
    (
        this->lookup("normalMotionDir")
    ),
    cleanInterface_
    (
        this->lookup("cleanInterface")
    ),
    aPatchID_(-1),
    bPatchID_(-1),
    muFluidA_
    (
        this->lookup("muFluidA")
    ),
    muFluidB_
    (
        this->lookup("muFluidB")
    ),
    muFluidB1_
    (
        this->lookup("muFluidB1")
    ),
    rhoFluidA_
    (
        this->lookup("rhoFluidA")
    ),
    rhoFluidB_
    (
        this->lookup("rhoFluidB")
    ),
    rhoFluidB1_
    (
        this->lookup("rhoFluidB1")
    ),

    ///RG-TQ conductivity,Cp,Latent heat,pressure_corr,Volume, gas constant, etc..
    kFluidA_
    (
        this->lookup("kFluidA")
    ),
    ///RG-TQ
    kFluidB_
    (
        this->lookup("kFluidB")
    ),
    kFluidB1_
    (
        this->lookup("kFluidB1")
    ),
    kWall_
    (
        this->lookup("kWall")
    ),
    thWall_
    (
        this->lookup("thWall")
    ),

    TLeft_
    (
        this->lookup("TLeft")
    ),
    TRight_
    (
        this->lookup("TRight")
    ),
    acmCoef_
    (
        this->lookup("acmCoef")
    ),

    ///RG-TQ
    CpFluidA_
    (
        this->lookup("CpFluidA")
    ),
    ///RG-TQ
    CpFluidB_
    (
        this->lookup("CpFluidB")
    ),
    CpFluidB1_
    (
        this->lookup("CpFluidB1")
    ),
    betaFluidA_
    (
        this->lookup("betaFluidA")
    ),
    /*betaFluidB_
    (
        this->lookup("betaFluidB")
    ),*/

    ///RG-TQ
    DfFluidA_
    (
        this->lookup("DfFluidA")
    ),
    ///RG-TQ
    DfFluidB_
    (
        this->lookup("DfFluidB")
    ),
    DfFluidB0_
    (
        this->lookup("DfFluidB0")
    ),
    antoineA_
    (
        this->lookup("antoineA")
    ),
    antoineB_
    (
        this->lookup("antoineB")
    ),
    antoineC_
    (
        this->lookup("antoineC")
    ),

    pTotal_
    (
        this->lookup("pTotal")
    ),

    c0_
    (
        this->lookup("c0")
    ),

    ///RG-TQ
    pCorr_("pCorr", dimPressure, 0.0),

    MB1_("MB1", dimMass, 0.0),
    MBD_("MBD", dimMass, 0.0),

    Mair0_
    (
        this->lookup("Mair0")
    ),

    ///RG-TQ
    VvInitial_("VvInitial", dimVolume, 0.0),
    MvInitial_("MvInitial", dimMass, 0.0),
    ///RG-TQ
    Vv_("Vv", dimVolume, 0.0),
    Mv_("Mv", dimMass, 0.0),

    molarMass1_
    (
        this->lookup("molarMass1")
    ),

    molarMassD_
    (
        this->lookup("molarMassD")
    ),

    ///RG-TQ
    latentHeat_
    (
        this->lookup("latentHeat")
    ),
    ///RG-TQ
    gasConstant1_
    (
        this->lookup("gasConstant1")
    ),
    ///RG-TQ
    gasConstantD_
    (
        this->lookup("gasConstantD")
    ),

    ///RG-TQ
    TRef_
    (
        this->lookup("referenceTemperature")
    ),
//    rhoVaporRef_
//    (
//        this->lookup("rhoVaporRef")
//    ),

    g_(this->lookup("g")),

    ///RG-TQ
    ctangle_(this->lookup("ctangle")),

    cleanInterfaceSurfTension_
    (
        this->lookup("surfaceTension")
    ),
    tempCoeffSurfTension_
    (
        this->lookup("tempCoeffSurfTension")
    ),
    fixedFreeSurfacePatches_
    (
        this->lookup("fixedFreeSurfacePatches")
    ),
    pointNormalsCorrectionPatches_
    (
        this->lookup("pointNormalsCorrectionPatches")
    ),
//z0312s
    nFreeSurfCorr_
    (
        readInt(this->lookup("nFreeSurfaceCorrectors"))
    ),
//z0312e
//    dTheta_(0.0),
    smoothing_(false),
//    z0412s
    correctPointNormals_(false),
    correctDisplacement_(false),
    correctCurvature_(false),
    curvExtrapOrder_(0),
//    z0412e
    fvcNGradUn_(false),

    interpolatorABPtr_(NULL),
    interpolatorBAPtr_(NULL),
    controlPointsPtr_(NULL),

    TiPtr_(NULL),
    TiDiffPtr_(NULL),
    dTiPtr_(NULL),    
    gradTiPtr_(NULL), 
    massFluxPtr_(NULL),
    heatFluxlPtr_(NULL),      
    heatFluxvPtr_(NULL),    
    UiPtr_(NULL),
    conciPtr_(NULL),   
            
    motionPointsMaskPtr_(NULL),
    pointsDisplacementDirPtr_(NULL),
    facesDisplacementDirPtr_(NULL),
    totalDisplacementPtr_(NULL),
    aMeshPtr_(NULL),
    UsPtr_(NULL),
    phisPtr_(NULL),
    phiPtr_(NULL),
    ddtPhiPtr_(NULL),
    surfactConcPtr_(NULL),
    surfaceTensionPtr_(NULL),
    surfactantPtr_(NULL),
    fluidIndicatorPtr_(NULL),
    contactAnglePtr_(NULL),
    nGradUnPtr_(NULL) //0917-2012

{
    // Set point normal correction patches
    boolList& correction = aMesh().correctPatchPointNormals();

    forAll(pointNormalsCorrectionPatches_, patchI)
    {
        word patchName = pointNormalsCorrectionPatches_[patchI];

        label patchID = aMesh().boundary().findPatchID(patchName);

        if(patchID == -1)
        {
            FatalErrorIn
            (
                "freeSurface::freeSurface(...)"
            )   << "patch name for point normals correction don't exist"
                << abort(FatalError);
        }

        correction[patchID] = true;
    }

    // Detect the free surface patch
    forAll (mesh().boundary(), patchI)
    {
        if(mesh().boundary()[patchI].name() == "freeSurface")
        {
            aPatchID_ = patchI;
                
            Info<< "Found free surface patch. ID: " << aPatchID_
                << endl;
        }
    }

    if(aPatchID() == -1)
    {
        FatalErrorIn("freeSurface::freeSurface(...)")
            << "Free surface patch not defined.  Please make sure that "
                << " the free surface patches is named as freeSurface"
                << abort(FatalError);
    }


    // Detect the free surface shadow patch
    if (twoFluids())
    {
        forAll (mesh().boundary(), patchI)
        {
            if(mesh().boundary()[patchI].name() == "freeSurfaceShadow")
            {
                bPatchID_ = patchI;
                    
                Info<< "Found free surface shadow patch. ID: "
                    << bPatchID_ << endl;
            }
        }

        if(bPatchID() == -1)
        {
            FatalErrorIn("freeSurface::freeSurface(...)")
                << "Free surface shadow patch not defined. "
                    << "Please make sure that the free surface shadow patch "
                    << "is named as freeSurfaceShadow."
                    << abort(FatalError);
        }
    }


    // Mark free surface boundary points 
    // which belonge to processor patches
    forAll(aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type()
         == processorFaPatch::typeName
        )
        {
            const labelList& patchPoints =
                aMesh().boundary()[patchI].pointLabels();

            forAll(patchPoints, pointI)
            {
                motionPointsMask()[patchPoints[pointI]] = -1;
            }                
        }
    }


    // Mark fixed free surface boundary points 
    forAll(fixedFreeSurfacePatches_, patchI)
    {
        label fixedPatchID = 
            aMesh().boundary().findPatchID
            (
                fixedFreeSurfacePatches_[patchI]
            );

        if(fixedPatchID == -1)
        {
            FatalErrorIn("freeSurface::freeSurface(...)")
                << "Wrong faPatch name in the fixedFreeSurfacePatches list"
                    << " defined in the freeSurfaceProperties dictionary"
                    << abort(FatalError);
        }

        const labelList& patchPoints =
            aMesh().boundary()[fixedPatchID].pointLabels();

        forAll(patchPoints, pointI)
        {
            motionPointsMask()[patchPoints[pointI]] = 0;
        }
    }


    // Mark free-surface boundary point 
    // at the axis of 2-D axisymmetic cases
    forAll(aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type()
         == wedgeFaPatch::typeName
        )
        {
            const wedgeFaPatch& wedgePatch =
                refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

            if(wedgePatch.axisPoint() > -1)
            {
                motionPointsMask()[wedgePatch.axisPoint()] = 0;

                Info << "Axis point: " 
                    << wedgePatch.axisPoint()
                    << "vector: " 
                    << aMesh().points()[wedgePatch.axisPoint()] << endl;
            }
        }
    }


    // Read free-surface points total displacement if present
    readTotalDisplacement();


    // Read control points positions if present
    controlPoints();


    // Check if smoothing switch is set
    if (this->found("smoothing"))
    {
        smoothing_ = Switch(this->lookup("smoothing"));
        Info<< "smoothing: "<< smoothing_ <<endl ;
    }

    ///RG-TQ Set interface flux
    phi();

    ///RG-TQ There should be a prettier way to initialize the mass flux
    J = 0.0*aMesh().S();
    Vn = 0.0*aMesh().S();
    ///RG-TQ Initialize the vapor density
    rhoB1 = 0.0*aMesh().S();
    ncB1b = 0.0*aMesh().S();
    ncBDb = 0.0*aMesh().S();
    ncB1 = 0.0*aMesh().S();
    ncBD = 0.0*aMesh().S();
    gradn1B = 0.0*aMesh().S();
//    rhoB1 += rhoFluidB1().value();

    makecontactAngle();

/*
    // Check if contactAngle is defined
    IOobject contactAngleHeader
    (
        "contactAngle",
        DB().timeName(),
        mesh(),
        IOobject::MUST_READ
    );
    if (contactAngleHeader.headerOk())
    {
        Info << "Reading contact angle" << endl;

        contactAnglePtr_ =
//            new edgeVectorField
            new edgeScalarField
            (
                IOobject
                (
                    "contactAngle",
                    DB().timeName(),
                    mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                aMesh()
            );
        Info << "Reading contact angle finished" << endl;
    }
//    z0412e
//    correctContactLinePointNormals();
*/
//    z0412s
    // Check if correctPointNormals switch is set
    if (this->found("correctPointNormals"))
    {
        correctPointNormals_ = Switch(this->lookup("correctPointNormals"));
    }

    // Check if correctDisplacement switch is set
    if (this->found("correctDisplacement"))
    {
        correctDisplacement_ = Switch(this->lookup("correctDisplacement"));
    }

    // Check if correctCurvature switch is set
    if (this->found("correctCurvature"))
    {
        correctCurvature_ = Switch(this->lookup("correctCurvature"));
    }

    // Check if curvExtrapOrder parameter is set
    if (this->found("curvExtrapOrder"))
    {
        curvExtrapOrder_ = Switch(this->lookup("curvExtrapOrder"));
    }

    correctContactLinePointNormals();
//    z0412e

    // 0915-2012:Check if fvcNGradUn switch is set
    if (this->found("fvcNGradUn"))
    {
        fvcNGradUn_ = Switch(this->lookup("fvcNGradUn"));
    }
}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

freeSurface::~freeSurface()
{
    clearOut();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void freeSurface::updateDisplacementDirections()









{
    if(normalMotionDir())
    {
        // Update point displacement correction
        pointsDisplacementDir() = aMesh().pointAreaNormals();

        // Correcte point displacement direction 
        // at the "centerline" symmetryPlane which represents the axis
        // of an axisymmetric case
        forAll(aMesh().boundary(), patchI)
        {
            if(aMesh().boundary()[patchI].type() == wedgeFaPatch::typeName)
            {
                const wedgeFaPatch& wedgePatch =
                    refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);
		
                vector axis = wedgePatch.axis();

                label centerLinePatchID = 
                    aMesh().boundary().findPatchID("centerline");
	
                if(centerLinePatchID != -1)
                {
                    const labelList& pointLabels = 
                        aMesh().boundary()[centerLinePatchID].pointLabels();
                    
                    forAll(pointLabels, pointI)
                    {
                        vector dir = 
                            pointsDisplacementDir()[pointLabels[pointI]];
                    
                        dir = (dir&axis)*axis;
                        dir /= mag(dir);
                
                        pointsDisplacementDir()[pointLabels[pointI]] = dir;
                    }
                }
                else
                {
                    Info << "Warning: centerline polyPatch does not exist. " 
                        << "Free surface points displacement directions "
                        << "will not be corrected at the axis (centerline)" 
                        << endl; 
                }
            
                break;   
            }
        }

        // Update face displacement direction
        facesDisplacementDir() =
            aMesh().faceAreaNormals().internalField();

        // Correction of control points postion
        const vectorField& Cf = aMesh().areaCentres().internalField();

        controlPoints() = 
            facesDisplacementDir()
           *(facesDisplacementDir()&(controlPoints() - Cf))
          + Cf;
    }
}
//z0312 freeSurface::predictPoints()
bool freeSurface::predictPoints()
{
    // Smooth interface
//Info << "before smoothing === "<< endl;
    if (smoothing_)
    {
        smoothing();
    }
//Info << "after smoothing === "<< endl;
    ///RG-TQ For incoporating phase change
    const scalarField& Sf = aMesh().S();

    for
    (
        int freeSurfCorr=0;
        freeSurfCorr<nFreeSurfCorr_;
        freeSurfCorr++
    )
    {
        movePoints(phi_.boundaryField()[aPatchID()]
               ///RG-TQ Considering phase change
               - J/rhoFluidA().value()*Sf
		);
//         moveMeshPoints();
    }

    return true;
}

//z0312 freeSurface::correctPoints()
bool freeSurface::correctPoints()
{
    ///RG-TQ For incoporating phase change
    const scalarField& Sf = aMesh().S();

    for
    (
        int freeSurfCorr=0;
        freeSurfCorr<nFreeSurfCorr_;
        freeSurfCorr++
    )
    {
        movePoints(phi_.boundaryField()[aPatchID()]
                  ///RG-TQ Considering phase change
                  - J/rhoFluidA().value()*Sf
                  );
//         moveMeshPoints();
    }

    return true;
}

//z0312 freeSurface::movePoints(const scalarField& interfacePhi)
bool freeSurface::movePoints(const scalarField& interfacePhi)
{
    pointField newMeshPoints = mesh().allPoints();

    scalarField sweptVolCorr = 
        interfacePhi
      - fvc::meshPhi(rho(),U())().boundaryField()[aPatchID()];

    word ddtScheme
    (
        mesh().ddtScheme("ddt(" + rho().name() + ',' + U().name()+')')
    );

    if 
    (
        ddtScheme
     == fv::CrankNicholsonDdtScheme<vector>::typeName
    )
    {
        sweptVolCorr *= (1.0/2.0)*DB().deltaT().value();
    }
    else if
    (
        ddtScheme
     == fv::EulerDdtScheme<vector>::typeName
    )
    {
        sweptVolCorr *= DB().deltaT().value();
    }
    else if
    (
        ddtScheme
     == fv::backwardDdtScheme<vector>::typeName
    )
    {
        if (DB().timeIndex() == 1)
        {
            sweptVolCorr *= DB().deltaT().value();
        }
        else
        {
            sweptVolCorr *= (2.0/3.0)*DB().deltaT().value();
        }
    }   
    else
    {
        FatalErrorIn("freeSurface::movePoints()")
            << "Unsupported temporal differencing scheme : "
                << ddtScheme
                << abort(FatalError);
    }

    const scalarField& Sf = aMesh().S();
    const vectorField& Nf = aMesh().faceAreaNormals().internalField();

    scalarField deltaH =
        sweptVolCorr/(Sf*(Nf & facesDisplacementDir()));

    forAll(fixedFreeSurfacePatches_, patchI)
    {
        label fixedPatchID = 
            aMesh().boundary().findPatchID
            (
                fixedFreeSurfacePatches_[patchI]
            );

        if(fixedPatchID == -1)
        {
            FatalErrorIn("freeSurface::freeSurface(...)")
                << "Wrong faPatch name in the fixedFreeSurfacePatches list"
                    << " defined in the freeSurfaceProperties dictionary"
                    << abort(FatalError);
        }
        
        const labelList& eFaces =
            aMesh().boundary()[fixedPatchID].edgeFaces();

        forAll(eFaces, edgeI)
        {
            deltaH[eFaces[edgeI]] *= 2.0;
        }
    }

    pointField displacement = pointDisplacement(deltaH);

    if (correctDisplacement_)
    {
        correctPointDisplacement(sweptVolCorr, displacement);
    }

    // Move only free-surface points

    const labelList& meshPointsA = 
        mesh().boundaryMesh()[aPatchID()].meshPoints();

    forAll (displacement, pointI)
    {
        newMeshPoints[meshPointsA[pointI]] += displacement[pointI];
    }

    if(twoFluids_)
    {
        const labelList& meshPointsB = 
            mesh().boundaryMesh()[bPatchID_].meshPoints();

        pointField displacementB =
            interpolatorAB().pointInterpolate
            (
                displacement
            );
        
        forAll (displacementB, pointI)
        {
            newMeshPoints[meshPointsB[pointI]] += displacementB[pointI]; 
        }
    }

    // Update total displacement field

    if(totalDisplacementPtr_ && (curTimeIndex_ < DB().timeIndex()))
    {
        FatalErrorIn("freeSurface::movePoints()")
            << "Total displacement of free surface points "
                << "from previous time step is not absorbed by the mesh."
                << abort(FatalError);
    }
    else if (curTimeIndex_ < DB().timeIndex())
    {
        totalDisplacement() = displacement;

        curTimeIndex_ = DB().timeIndex();
    }
    else
    {
        totalDisplacement() += displacement;
    }

    twoDPointCorrector twoDPointCorr(mesh());

    twoDPointCorr.correctPoints(newMeshPoints);

    mesh().movePoints(newMeshPoints);

    aMesh().movePoints(mesh().points());

    correctContactLinePointNormals();

    if (correctPointNormals_)
    {
        correctPointNormals();
    }

    if (correctCurvature_)
    {
//        correctCurvature();
        smoothCurvature();
    }


    // Move correctedFvPatchField fvSubMeshes

    forAll(U().boundaryField(), patchI)
    {
        if
        (
            (
                U().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<vector>::typeName
            )
            ||
            (
                U().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<vector>::typeName
            )
            ||
            (
                U().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<vector>::typeName
            )
        )
        {
            correctedFvPatchField<vector>& pU =
                refCast<correctedFvPatchField<vector> >
                (
                    U().boundaryField()[patchI]
                );

            pU.movePatchSubMesh();
        }
    }

    forAll(p().boundaryField(), patchI)
    {
        if
        (
            (
                p().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<scalar>::typeName
            )
            ||
            (
                p().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<scalar>::typeName
            )
            ||
            (
                p().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<scalar>::typeName
            )
        )
        {
            correctedFvPatchField<scalar>& pP =
                refCast<correctedFvPatchField<scalar> >
                (
                    p().boundaryField()[patchI]
                );

            pP.movePatchSubMesh();
        }
    }

    return true;
}

//z0312 freeSurface::moveMeshPointsForOldFreeSurfDisplacement(
bool freeSurface::moveMeshPointsForOldFreeSurfDisplacement()
{
    if(totalDisplacementPtr_)
    {
        pointField newPoints = mesh().allPoints();

        const labelList& meshPointsA = 
            mesh().boundaryMesh()[aPatchID()].meshPoints();

        forAll (totalDisplacement(), pointI)
        {
            newPoints[meshPointsA[pointI]] -= totalDisplacement()[pointI]; 
        }


        // Check mesh motion solver type 

        bool feMotionSolver = 
            mesh().objectRegistry::foundObject<tetPointVectorField>
            (
                "motionU"
            );
    
        bool fvMotionSolver =
            mesh().objectRegistry::foundObject<pointVectorField>
            (
                "pointMotionU"
            );

        if (feMotionSolver)
        {
            tetPointVectorField& motionU =
                const_cast<tetPointVectorField&>
                (
                    mesh().objectRegistry::
                    lookupObject<tetPointVectorField>
                    (
                        "motionU"
                    )
                );

            fixedValueTetPolyPatchVectorField& motionUaPatch =
                refCast<fixedValueTetPolyPatchVectorField>
                (
                    motionU.boundaryField()[aPatchID()]
                );

            tetPolyPatchInterpolation tppiAPatch
            (
                refCast<const faceTetPolyPatch>
                (
                    motionUaPatch.patch()
                )
            );

            motionUaPatch ==
                tppiAPatch.pointToPointInterpolate
                (
                    totalDisplacement()/DB().deltaT().value()
                );

            if(twoFluids_)
            {
                const labelList& meshPointsB = 
                    mesh().boundaryMesh()[bPatchID()].meshPoints();

                pointField totDisplacementB =
                    interpolatorAB().pointInterpolate
                    (
                        totalDisplacement()
                    );
        
                forAll (totDisplacementB, pointI)
                {
                    newPoints[meshPointsB[pointI]] -= 
                        totDisplacementB[pointI]; 
                }
            
                fixedValueTetPolyPatchVectorField& motionUbPatch =
                    refCast<fixedValueTetPolyPatchVectorField>
                    (
                        motionU.boundaryField()[bPatchID()]
                    );

                tetPolyPatchInterpolation tppiBPatch
                (
                    refCast<const faceTetPolyPatch>(motionUbPatch.patch())
                );

                motionUbPatch == 
                    tppiBPatch.pointToPointInterpolate
                    (
                        totDisplacementB/DB().deltaT().value()
                    );
            }
        }
        else if (fvMotionSolver)
        {
            pointVectorField& motionU =
                const_cast<pointVectorField&>
                (
                    mesh().objectRegistry::
                    lookupObject<pointVectorField>
                    (
                        "pointMotionU"
                    )
                );

            fixedValuePointPatchVectorField& motionUaPatch =
                refCast<fixedValuePointPatchVectorField>
                (
                    motionU.boundaryField()[aPatchID()]
                );

            motionUaPatch ==
                totalDisplacement()/DB().deltaT().value();

            if(twoFluids_)
            {
                const labelList& meshPointsB = 
                    mesh().boundaryMesh()[bPatchID()].meshPoints();

                pointField totDisplacementB =
                    interpolatorAB().pointInterpolate
                    (
                        totalDisplacement()
                    );
        
                forAll (totDisplacementB, pointI)
                {
                    newPoints[meshPointsB[pointI]] -= 
                        totDisplacementB[pointI]; 
                }
            
                fixedValuePointPatchVectorField& motionUbPatch =
                    refCast<fixedValuePointPatchVectorField>
                    (
                        motionU.boundaryField()[bPatchID()]
                    );

                motionUbPatch == 
                    totDisplacementB/DB().deltaT().value();
            }
        }

        twoDPointCorrector twoDPointCorr(mesh());

        twoDPointCorr.correctPoints(newPoints);

        mesh().movePoints(newPoints);

        deleteDemandDrivenData(totalDisplacementPtr_);

        mesh().update();

        aMesh().movePoints(mesh().points());

        correctContactLinePointNormals();

        if (correctPointNormals_)
        {
            correctPointNormals();
        }

        if (correctCurvature_)
        {
//            correctCurvature();
            smoothCurvature();
        }


        // Move correctedFvPatchField fvSubMeshes

        forAll(U().boundaryField(), patchI)
        {
            if
            (
                (
                    U().boundaryField()[patchI].type()
                 == fixedGradientCorrectedFvPatchField<vector>::typeName
                )
                ||
                (
                    U().boundaryField()[patchI].type()
                 == fixedValueCorrectedFvPatchField<vector>::typeName
                )
                ||
                (
                    U().boundaryField()[patchI].type()
                 == zeroGradientCorrectedFvPatchField<vector>::typeName
                )
            )
            {
                correctedFvPatchField<vector>& aU =
                    refCast<correctedFvPatchField<vector> >
                    (
                        U().boundaryField()[patchI]
                    );

                aU.movePatchSubMesh();
            }
        }

        forAll(p().boundaryField(), patchI)
        {
            if
            (
                (
                    p().boundaryField()[patchI].type()
                 == fixedGradientCorrectedFvPatchField<scalar>::typeName
                )
                ||
                (
                    p().boundaryField()[patchI].type()
                 == fixedValueCorrectedFvPatchField<scalar>::typeName
                )
                ||
                (
                    p().boundaryField()[patchI].type()
                 == zeroGradientCorrectedFvPatchField<scalar>::typeName
                )
            )
            {
                correctedFvPatchField<scalar>& aP =
                    refCast<correctedFvPatchField<scalar> >
                    (
                        p().boundaryField()[patchI]
                    );

                aP.movePatchSubMesh();
            }
        }
    }

    return true;
}

//z0312 freeSurface::moveMeshPoints()

bool freeSurface::moveMeshPoints()
{
    scalarField sweptVolCorr = 
        phi_.boundaryField()[aPatchID()]
      - fvc::meshPhi(rho(),U())().boundaryField()[aPatchID()];

    word ddtScheme
    (
        mesh().ddtScheme("ddt(" + rho().name() + ',' + U().name()+')')


    );

    if 
    (
        ddtScheme
     == fv::CrankNicholsonDdtScheme<vector>::typeName
    )
    {
        sweptVolCorr *= (1.0/2.0)*DB().deltaT().value();
    }
    else if
    (
        ddtScheme
     == fv::EulerDdtScheme<vector>::typeName
    )
    {
        sweptVolCorr *= DB().deltaT().value();
    }
    else if
    (
        ddtScheme
     == fv::backwardDdtScheme<vector>::typeName
    )
    {
        sweptVolCorr *= (2.0/3.0)*DB().deltaT().value();
    }   
    else
    {
        FatalErrorIn("freeSurface::movePoints()")
            << "Unsupported temporal differencing scheme : "
                << ddtScheme
                << abort(FatalError);
    }

    const scalarField& Sf = aMesh().S();
    const vectorField& Nf = aMesh().faceAreaNormals().internalField();

    scalarField deltaH =
        sweptVolCorr/(Sf*(Nf & facesDisplacementDir()));

    forAll(fixedFreeSurfacePatches_, patchI)
    {
        label fixedPatchID = 
            aMesh().boundary().findPatchID
            (
                fixedFreeSurfacePatches_[patchI]
            );

        if(fixedPatchID == -1)
        {
            FatalErrorIn("freeSurface::freeSurface(...)")
                << "Wrong faPatch name in the fixedFreeSurfacePatches list"
                    << " defined in the freeSurfaceProperties dictionary"
                    << abort(FatalError);
        }
        
        const labelList& eFaces =
            aMesh().boundary()[fixedPatchID].edgeFaces();

        forAll(eFaces, edgeI)
        {
            deltaH[eFaces[edgeI]] *= 2.0;
        }
    }

    pointField displacement = pointDisplacement(deltaH);

    //-- Set mesh motion boundary conditions

    // Check mesh motion solver type 

    bool feMotionSolver = 
        mesh().objectRegistry::foundObject<tetPointVectorField>
        (
            "motionU"
        );
    
    bool fvMotionSolver =
        mesh().objectRegistry::foundObject<pointVectorField>
        (
            "pointMotionU"
        );

    if (feMotionSolver)
    {
        tetPointVectorField& motionU =
            const_cast<tetPointVectorField&>
            (
                mesh().objectRegistry::
                lookupObject<tetPointVectorField>
                (
                    "motionU"
                )
            );

        fixedValueTetPolyPatchVectorField& motionUaPatch =
            refCast<fixedValueTetPolyPatchVectorField>
            (
                motionU.boundaryField()[aPatchID()]
            );
        
        tetPolyPatchInterpolation tppiAPatch
        (
            refCast<const faceTetPolyPatch>
            (
                motionUaPatch.patch()
            )
        );

        motionUaPatch ==
            tppiAPatch.pointToPointInterpolate
            (
                displacement/DB().deltaT().value()
            );

        if(twoFluids_)
        {
            pointField displacementB =
                interpolatorAB().pointInterpolate
                (
                    displacement
                );
        
            fixedValueTetPolyPatchVectorField& motionUbPatch =
                refCast<fixedValueTetPolyPatchVectorField>
                (
                    motionU.boundaryField()[bPatchID()]
                );

            tetPolyPatchInterpolation tppiBPatch
            (
                refCast<const faceTetPolyPatch>(motionUbPatch.patch())
            );

            motionUbPatch == 
                tppiBPatch.pointToPointInterpolate
                (
                    displacement/DB().deltaT().value()
                );
        }
    }
    else if (fvMotionSolver)
    {
        pointVectorField& motionU =
            const_cast<pointVectorField&>
            (
                mesh().objectRegistry::
                lookupObject<pointVectorField>
                (
                    "pointMotionU"
                )
            );

        fixedValuePointPatchVectorField& motionUaPatch =
            refCast<fixedValuePointPatchVectorField>
            (
                motionU.boundaryField()[aPatchID()]
            );
        
        motionUaPatch ==
            displacement/DB().deltaT().value();

        if(twoFluids_)
        {
            pointField displacementB =
                interpolatorAB().pointInterpolate
                (
                    displacement
                );
            
            fixedValuePointPatchVectorField& motionUbPatch =
                refCast<fixedValuePointPatchVectorField>
                (
                    motionU.boundaryField()[bPatchID()]
                );

            motionUbPatch == 
                displacementB/DB().deltaT().value();
        }
    }

    mesh().update();

    aMesh().movePoints(mesh().points());

    correctContactLinePointNormals();

    if (correctPointNormals_)
    {
        correctPointNormals();
    }

    if (correctCurvature_)
    {
//        correctCurvature();
        smoothCurvature();
    }

    // Move correctedFvPatchField fvSubMeshes

    forAll(U().boundaryField(), patchI)
    {
        if
        (
            (
                U().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<vector>::typeName
            )
         ||
            (
                U().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<vector>::typeName
            )
         ||
            (
                U().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<vector>::typeName
            )
        )
        {
            correctedFvPatchField<vector>& aU =
                refCast<correctedFvPatchField<vector> >
                (
                    U().boundaryField()[patchI]
                );

            aU.movePatchSubMesh();
        }
    }

    forAll(p().boundaryField(), patchI)
    {
        if
        (
            (
                p().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<scalar>::typeName
            )
         ||
            (
                p().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<scalar>::typeName
            )
         ||
            (
                p().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<scalar>::typeName
            )
        )
        {
            correctedFvPatchField<scalar>& aP =
                refCast<correctedFvPatchField<scalar> >
                (
                    p().boundaryField()[patchI]
                );

                aP.movePatchSubMesh();
        }
    }

    return true;
}
/*
bool freeSurface::predictPoints()
{
    // Smooth interface

    if (smoothing_)
    {
        Info<< "Perform smoothing " <<endl ;
        controlPoints() = aMesh().areaCentres().internalField();
        movePoints(scalarField(controlPoints().size(), 0));
        movePoints(-fvc::meshPhi(rho(),U())().boundaryField()[aPatchID()]);
    }

    // Predict points position

//     areaScalarField predictedPhi = phi().oldTime();

//     areaScalarField predictedPhi = 
//         phi().oldTime()
//       + (
//             (3.0/2.0)*ddtPhi().oldTime() 
//           - (1.0/2.0)*ddtPhi().oldTime().oldTime()
//         )
//        *DB().deltaT();

//     if (DB().timeIndex() == 1)
//     {


//         predictedPhi = 
//             phi().oldTime()
//           + ddtPhi().oldTime()*DB().deltaT();
//     }

    ///RG-TQ For incoporating phase change
    const scalarField& Sf = aMesh().S();

    movePoints(
               phi_.boundaryField()[aPatchID()]
               ///RG-TQ Considering phase change
               - J/rhoFluidA().value()*Sf
               );
    return true;
}


bool freeSurface::correctPoints()
{
    label nFreeSurfCorr_ = 1;

    ///RG-TQ For incoporating phase change
    const scalarField& Sf = aMesh().S();

    for
    (
        int freeSurfCorr=0;
        freeSurfCorr<nFreeSurfCorr_;
        freeSurfCorr++
    )

    {
        movePoints(
                  phi_.boundaryField()[aPatchID()]
                  ///RG-TQ Considering phase change
                  - J/rhoFluidA().value()*Sf
                  );
    }
    
//     ddtPhi() = fac::ddt(phi());

    return true;
}

//z0312 freeSurface::movePoints(const scalarField& interfacePhi)
bool freeSurface::movePoints(const scalarField& interfacePhi)
{
//    pointField newMeshPoints = mesh().points();
    pointField newMeshPoints = mesh().allPoints();

    scalarField sweptVolCorr = 
        interfacePhi
      - fvc::meshPhi(rho(),U())().boundaryField()[aPatchID()];

    word ddtScheme
    (
        mesh().ddtScheme("ddt(" + rho().name() + ',' + U().name()+')')
    );

    if 
    (
        ddtScheme
     == fv::CrankNicholsonDdtScheme<vector>::typeName
    )
    {
        sweptVolCorr *= (1.0/2.0)*DB().deltaT().value();
    }
    else if
    (
        ddtScheme
     == fv::EulerDdtScheme<vector>::typeName
    )
    {
        sweptVolCorr *= DB().deltaT().value();
    }
    else if
    (
        ddtScheme
     == fv::backwardDdtScheme<vector>::typeName
    )
    {
        if (DB().timeIndex() == 1)
        {
            sweptVolCorr *= DB().deltaT().value();
        }
        else
        {
            sweptVolCorr *= (2.0/3.0)*DB().deltaT().value();
        }
    }   
    else
    {
        FatalErrorIn("freeSurface::movePoints()")
            << "Unsupported temporal differencing scheme : "
                << ddtScheme
                << abort(FatalError);
    }

    const scalarField& Sf = aMesh().S();
    const vectorField& Nf = aMesh().faceAreaNormals().internalField();

    scalarField deltaH =
        sweptVolCorr/(Sf*(Nf & facesDisplacementDir()));

    forAll(fixedFreeSurfacePatches_, patchI)
    {
        label fixedPatchID = 
            aMesh().boundary().findPatchID
            (
                fixedFreeSurfacePatches_[patchI]
            );

        if(fixedPatchID == -1)
        {
            FatalErrorIn("freeSurface::freeSurface(...)")
                << "Wrong faPatch name in the fixedFreeSurfacePatches list"
                    << " defined in the freeSurfaceProperties dictionary"
                    << abort(FatalError);
        }
        
        const labelList& eFaces =
            aMesh().boundary()[fixedPatchID].edgeFaces();

        forAll(eFaces, edgeI)
        {
            deltaH[eFaces[edgeI]] *= 2.0;
        }
    }

    pointField displacement = pointDisplacement(deltaH);

//    if (correctDisplacement_)
//    {
//       correctPointDisplacement(sweptVolCorr, displacement);
//    }

    // Move only free-surface points

    const labelList& meshPointsA = 
        mesh().boundaryMesh()[aPatchID()].meshPoints();

    forAll (displacement, pointI)
    {
        newMeshPoints[meshPointsA[pointI]] += displacement[pointI];
    }

    if(twoFluids_)
    {
        const labelList& meshPointsB = 
            mesh().boundaryMesh()[bPatchID_].meshPoints();

        pointField displacementB =
            interpolatorAB().pointInterpolate
            (
                displacement
            );
        
        forAll (displacementB, pointI)
        {
            newMeshPoints[meshPointsB[pointI]] += displacementB[pointI]; 
        }
    }

    // Update total displacement field

    if(totalDisplacementPtr_ && (curTimeIndex_ < DB().timeIndex()))
    {
        FatalErrorIn("freeSurface::movePoints()")
            << "Total displacement of free surface points "
                << "from previous time step is not absorbed by the mesh."
                << abort(FatalError);
    }
    else if (curTimeIndex_ < DB().timeIndex())
    {
        totalDisplacement() = displacement;

        curTimeIndex_ = DB().timeIndex();
    }
    else
    {
        totalDisplacement() += displacement;
    }

    twoDPointCorrector twoDPointCorr(mesh());

    twoDPointCorr.correctPoints(newMeshPoints);

    mesh().movePoints(newMeshPoints);

    aMesh().movePoints(mesh().points());

    correctContactLinePointNormals();

//    if (correctPointNormals_)
//    {
//        correctPointNormals();
//    }

//    if (correctCurvature_)
//    {
//        correctCurvature();
//    }


    // Move correctedFvPatchField fvSubMeshes

    forAll(U().boundaryField(), patchI)
    {
        if
        (
            (
                U().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<vector>::typeName
            )
            ||
            (
                U().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<vector>::typeName
            )
            ||
            (
                U().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<vector>::typeName
            )
        )
        {
            correctedFvPatchField<vector>& pU =
                refCast<correctedFvPatchField<vector> >
                (
                    U().boundaryField()[patchI]
                );

            pU.movePatchSubMesh();
        }
    }

    forAll(p().boundaryField(), patchI)
    {
        if
        (
            (
                p().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<scalar>::typeName
            )
            ||
            (
                p().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<scalar>::typeName
            )
            ||
            (
                p().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<scalar>::typeName
            )
        )
        {
            correctedFvPatchField<scalar>& pP =
                refCast<correctedFvPatchField<scalar> >
                (
                    p().boundaryField()[patchI]
                );

            pP.movePatchSubMesh();
        }
    }

    return true;
}


bool freeSurface::movePoints(const scalarField& interfacePhi)
{
    pointField newMeshPoints = mesh().points();

    scalarField sweptVolCorr = 
        interfacePhi
      - fvc::meshPhi(rho(),U())().boundaryField()[aPatchID()];

//    Info<< "interfacePhi"<<endl<<interfacePhi<<endl ;
//    Info<< "sweptVolCorr"<<endl<<sweptVolCorr<<endl ;

    word ddtScheme
    (
        mesh().ddtScheme("ddt(" + rho().name() + ',' + U().name()+')')
    );

    if 
    (
        ddtScheme
     == fv::CrankNicholsonDdtScheme<vector>::typeName
    )
    {
        sweptVolCorr *= (1.0/2.0)*DB().deltaT().value();
    }
    else if
    (
        ddtScheme
     == fv::EulerDdtScheme<vector>::typeName
    )
    {
        sweptVolCorr *= DB().deltaT().value();
    }
    else if
    (
        ddtScheme
     == fv::backwardDdtScheme<vector>::typeName
    )
    {
        if (DB().timeIndex() == 1)
        {
            sweptVolCorr *= DB().deltaT().value();
        }
        else
        {
            sweptVolCorr *= (2.0/3.0)*DB().deltaT().value();
        }
    }   
    else
    {
        FatalErrorIn("freeSurface::movePoints()")
            << "Unsupported temporal differencing scheme : "
                << ddtScheme
                << abort(FatalError);
    }


    const scalarField& Sf = aMesh().S();
    const vectorField& Nf = aMesh().faceAreaNormals().internalField();

    scalarField deltaH =
        sweptVolCorr/(Sf*(Nf & facesDisplacementDir()));

    pointField displacement = pointDisplacement(deltaH);

    ///RG-TQ For updating displacement
    displacementInitial = 0.0*displacement;

    // Move only free-surface points

    const labelList& meshPointsA = 

        mesh().boundaryMesh()[aPatchID()].meshPoints();

    forAll (displacement, pointI)
    {
        newMeshPoints[meshPointsA[pointI]] += displacement[pointI];
    }

    if(twoFluids_)
    {
        const labelList& meshPointsB = 
            mesh().boundaryMesh()[bPatchID_].meshPoints();

        pointField displacementB =
            interpolatorAB().pointInterpolate
            (
                displacement
            );
        
        forAll (displacementB, pointI)
        {
            newMeshPoints[meshPointsB[pointI]] += displacementB[pointI]; 
        }
    }

    // Update total displacement field

    if(totalDisplacementPtr_ && (curTimeIndex_ < DB().timeIndex()))
    {
        FatalErrorIn("freeSurface::movePoints()")
            << "Total displacement of free surface points "
                << "from previous time step is not absorbed by the mesh."
                << abort(FatalError);
    }
    else if (curTimeIndex_ < DB().timeIndex())
    {
        totalDisplacement() = displacement;

        curTimeIndex_ = DB().timeIndex();
    }
    else
    {
        totalDisplacement() += displacement;
    }

    twoDPointCorrector twoDPointCorr(mesh());

    twoDPointCorr.correctPoints(newMeshPoints);
    
    mesh().movePoints(newMeshPoints);

//RG-TQ
//    Info << "pointNormals (0) " << aMesh().pointAreaNormals() <<endl;

    aMesh().movePoints(mesh().points());

//RG-TQ
    correctContactLinePointNormals();
//    Info << "pointNormals (1) " << aMesh().pointAreaNormals() <<endl;

    // Move correctedFvPatchField fvSubMeshes

    forAll(U().boundaryField(), patchI)
    {
        if
        (
            (
                U().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<vector>::typeName
            )
            ||
            (
                U().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<vector>::typeName
            )
            ||
            (
                U().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<vector>::typeName
            )
        )
        {
            correctedFvPatchField<vector>& pU =
                refCast<correctedFvPatchField<vector> >
                (
                    U().boundaryField()[patchI]
                );

            pU.movePatchSubMesh();
        }
    }

    forAll(p().boundaryField(), patchI)
    {
        if
        (
            (
                p().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<scalar>::typeName
            )
            ||
            (
                p().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<scalar>::typeName
            )
            ||
            (
                p().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<scalar>::typeName
            )
        )
        {
            correctedFvPatchField<scalar>& pP =
                refCast<correctedFvPatchField<scalar> >
                (
                    p().boundaryField()[patchI]
                );

            pP.movePatchSubMesh();
        }
    }

    return true;
}
*/

///RG-TQ Function to Move whole free surface for single-phase volume conservation
bool freeSurface::moveWholeFreeSurface
(
//  dimensionedScalar Vv,
//  dimensionedScalar VvInitial
  dimensionedScalar Mv,//mvapor
  dimensionedScalar MvInitial //mair
)
{
//    Vv_ = Vv;
//    VvInitial_ = VvInitial;
    Mv_ = Mv;
    MvInitial_ = MvInitial;

    scalarField rhoPA = rho().boundaryField()[aPatchID()].patchInternalField();

    scalarField rhoPB = interpolatorBA().faceInterpolate
        ( 
            rho().boundaryField()[bPatchID()].patchInternalField()
        );

    scalarField conPB = interpolatorBA().faceInterpolate
        ( 
            con().boundaryField()[bPatchID()].patchInternalField()
        );
//    pointField newMeshPoints = mesh().points();
    pointField newMeshPoints = mesh().allPoints();
        
    const scalarField& Sf = aMesh().S();
    const vectorField& Nf = aMesh().faceAreaNormals().internalField();
    vector Verticaly=vector(0, 1, 0);
	
	dimensionedScalar A1 = fvc::domainIntegrate((1.0-fluidIndicator())*p()*con()/T())/gasConstant1();
	dimensionedScalar B1 = fvc::domainIntegrate((1.0-fluidIndicator())*con()/T())/gasConstant1();

	dimensionedScalar A2 = fvc::domainIntegrate((1.0-fluidIndicator())*p()*(1.0-con())/T())/gasConstantD();
	dimensionedScalar B2 = fvc::domainIntegrate((1.0-fluidIndicator())*(1.0-con())/T())/gasConstantD();

  	dimensionedScalar D = fvc::domainIntegrate(fluidIndicator()*rho());

	dimensionedScalar C1 = gSum((rhoPA - rhoPB*conPB)*Sf*(Verticaly& Nf));
	dimensionedScalar C2 = gSum(rhoPB*(1.0-conPB)*Sf*(Verticaly& Nf));

//    	scalar E = gSum(J*Sf)*DB().deltaT().value();

	vector discorr = -Verticaly*(A1.value() * B2.value() - A2.value() * B1.value() + B1.value() * MvInitial_.value()  - B2.value() * Mv_.value() + B2.value() * D.value()) / (B1.value() * C2.value() + B2.value() * C1.value());//+ B2.value() * E
	pCorr_.value() = -(A1.value() * C2.value() + A2.value() * C1.value() - C1.value() * MvInitial_.value()  - C2.value() * Mv_.value() + C2.value() * D.value()) / (B1.value() * C2.value() + B2.value() * C1.value());//+ C2.value() * E

//    scalar VCorr = gSum(((J)/rhoFluidA().value())*Sf)*DB().deltaT().value();
//    scalar JCorr = gSum(J*Sf)*DB().deltaT().value();

//    	scalarField Sfproj = Sf*(Verticaly& Nf);
//    	scalar Sftotal=gSum(Sfproj);
//    	scalar Diffrho=gAverage(rhoPA - rhoPB);
//    vector discorr = Verticaly*((VvInitial_.value() - Vv_.value() - VCorr)/Sftotal);
//    vector discorr = Verticaly*((MvInitial_.value() - Mv_.value() - JCorr*sign(JCorr))/Sftotal/Diffrho);

//    Info<< "Sftotal " << Sftotal << endl;
//    Info<< "Vv().value() " << Vv_.value() << endl;
//    Info<< "VvInitial().value() " << VvInitial_.value() << endl;
//    Info<< "VCorr " << VCorr << endl;
    Info<< "discorr " << discorr << endl;

    scalarField deltaH1 = scalarField(aMesh().nFaces(), 0.0);
    pointField displacement1 = pointDisplacement(deltaH1);

    ///RG-TQ For updating displacement
    displacementInitial = 0.0*displacement1;

    pointField displacement = displacementInitial + discorr;

    Info<< "Displacement " << displacement<< endl;
    Info<< "Mass accumulated from fluxes: " << gSum(J*Sf)*DB().deltaT().value() << endl;
    // Move only free-surface points

    const labelList& meshPointsA = 
        mesh().boundaryMesh()[aPatchID()].meshPoints();

    forAll (displacement, pointI)
    {
        newMeshPoints[meshPointsA[pointI]] += displacement[pointI];
    }

    if(twoFluids_)
    {
        const labelList& meshPointsB = 
            mesh().boundaryMesh()[bPatchID_].meshPoints();

        pointField displacementB =
            interpolatorAB().pointInterpolate
            (
                displacement
            );
        
        forAll (displacementB, pointI)
        {
            newMeshPoints[meshPointsB[pointI]] += displacementB[pointI]; 
        }
    }

    // Update total displacement field

    if(totalDisplacementPtr_ && (curTimeIndex_ < DB().timeIndex()))
    {
        FatalErrorIn("freeSurface::movePoints()")
            << "Total displacement of free surface points "
                << "from previous time step is not absorbed by the mesh."
                << abort(FatalError);
    }
    else if (curTimeIndex_ < DB().timeIndex())
    {
        totalDisplacement() = displacement;

        curTimeIndex_ = DB().timeIndex();
    }
    else
    {
        totalDisplacement() += displacement;
    }

    twoDPointCorrector twoDPointCorr(mesh());

    twoDPointCorr.correctPoints(newMeshPoints);

    mesh().movePoints(newMeshPoints);

    aMesh().movePoints(mesh().points());

    correctContactLinePointNormals();

    if (correctPointNormals_)
    {
        correctPointNormals();
    }

    if (correctCurvature_)
    {
//        correctCurvature();
        smoothCurvature();
    }

    // Move correctedFvPatchField fvSubMeshes

    forAll(U().boundaryField(), patchI)
    {
        if
        (
            (
                U().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<vector>::typeName
            )
            ||
            (
                U().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<vector>::typeName
            )
            ||
            (
                U().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<vector>::typeName
            )
        )
        {
            correctedFvPatchField<vector>& pU =
                refCast<correctedFvPatchField<vector> >
                (
                    U().boundaryField()[patchI]
                );

            pU.movePatchSubMesh();
        }
    }

    forAll(p().boundaryField(), patchI)
    {
        if
        (
            (
                p().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<scalar>::typeName
            )
            ||
            (
                p().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<scalar>::typeName
            )
            ||
            (
                p().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<scalar>::typeName
            )
        )
        {
            correctedFvPatchField<scalar>& pP =
                refCast<correctedFvPatchField<scalar> >
                (
                    p().boundaryField()[patchI]
                );

            pP.movePatchSubMesh();
        }
    }
/*
    // Move only free-surface points

    const labelList& meshPointsA = 
        mesh().boundaryMesh()[aPatchID()].meshPoints();

    forAll (displacement, pointI)
    {
        newMeshPoints[meshPointsA[pointI]] += discorr;
    //        newMeshPoints[meshPointsA[pointI]] += displacement[pointI];
    }

    if(twoFluids_)
    {
        const labelList& meshPointsB = 
            mesh().boundaryMesh()[bPatchID_].meshPoints();

        pointField displacementB =
            interpolatorAB().pointInterpolate
            (
                displacement
            );
        
        forAll (displacementB, pointI)
        {
            newMeshPoints[meshPointsB[pointI]] += discorr; 
    //            newMeshPoints[meshPointsB[pointI]] += displacementB[pointI]; 
        }
    }

    // Update total displacement field

    if(totalDisplacementPtr_ && (curTimeIndex_ < DB().timeIndex()))
    {
        FatalErrorIn("freeSurface::movePoints()")
            << "Total displacement of free surface points "
                << "from previous time step is not absorbed by the mesh."
                << abort(FatalError);
    }
    else if (curTimeIndex_ < DB().timeIndex())
    {
 
        totalDisplacement() = displacement;
        curTimeIndex_ = DB().timeIndex();
    }
    else
    {
        totalDisplacement() += displacement;
    }

    twoDPointCorrector twoDPointCorr(mesh());

    twoDPointCorr.correctPoints(newMeshPoints);
  
    mesh().movePoints(newMeshPoints);


    aMesh().movePoints(mesh().points());

//RG-TQ
    correctContactLinePointNormals();
//    Info << "pointNormals (3) " << aMesh().pointAreaNormals() <<endl;

    forAll(U().boundaryField(), patchI)
    {
        if
        (
            (
                U().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<vector>::typeName
            )
            ||
            (
                U().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<vector>::typeName
            )
            ||
            (
                U().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<vector>::typeName
            )
        )
        {
            correctedFvPatchField<vector>& pU =
                refCast<correctedFvPatchField<vector> >
                (
                    U().boundaryField()[patchI]
                );

            pU.movePatchSubMesh();
        }
    }

    forAll(p().boundaryField(), patchI)
    {
        if
        (
            (
                p().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<scalar>::typeName
            )
            ||
            (
                p().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<scalar>::typeName
            )
            ||
            (
                p().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<scalar>::typeName
            )
        )
        {
            correctedFvPatchField<scalar>& pP =
                refCast<correctedFvPatchField<scalar> >
                (
                    p().boundaryField()[patchI]
                );

            pP.movePatchSubMesh();
        }
    }
*/
    return true;
}
/*
//z0312 freeSurface::moveMeshPointsForOldFreeSurfDisplacement(
bool freeSurface::moveMeshPointsForOldFreeSurfDisplacement()
{
    if(totalDisplacementPtr_)
    {
//        pointField newPoints = mesh().points();
        pointField newPoints = mesh().allPoints();

        const labelList& meshPointsA = 
            mesh().boundaryMesh()[aPatchID()].meshPoints();

        forAll (totalDisplacement(), pointI)
        {
            newPoints[meshPointsA[pointI]] -= totalDisplacement()[pointI]; 
        }


        // Check mesh motion solver type 

        bool feMotionSolver = 
            mesh().objectRegistry::foundObject<tetPointVectorField>
            (
                "motionU"
            );
    
        bool fvMotionSolver =
            mesh().objectRegistry::foundObject<pointVectorField>
            (
                "pointMotionU"
            );

        if (feMotionSolver)
        {
//        Info << "s1feMotionSolver " <<  endl;
            tetPointVectorField& motionU =
                const_cast<tetPointVectorField&>
                (
                    mesh().objectRegistry::
                    lookupObject<tetPointVectorField>
                    (
                        "motionU"
                    )
                );

            fixedValueTetPolyPatchVectorField& motionUaPatch =
                refCast<fixedValueTetPolyPatchVectorField>
                (
                    motionU.boundaryField()[aPatchID()]
                );

            tetPolyPatchInterpolation tppiAPatch
            (
                refCast<const faceTetPolyPatch>
                (
                    motionUaPatch.patch()
                )
            );

            motionUaPatch ==
                tppiAPatch.pointToPointInterpolate
                (
                    totalDisplacement()/DB().deltaT().value()
                );

            if(twoFluids_)
            {
                const labelList& meshPointsB = 
                    mesh().boundaryMesh()[bPatchID()].meshPoints();

                pointField totDisplacementB =
                    interpolatorAB().pointInterpolate
                    (
                        totalDisplacement()
                    );
        
                forAll (totDisplacementB, pointI)
                {
                    newPoints[meshPointsB[pointI]] -= 
                        totDisplacementB[pointI]; 
                }
            
                fixedValueTetPolyPatchVectorField& motionUbPatch =
                    refCast<fixedValueTetPolyPatchVectorField>
                    (
                        motionU.boundaryField()[bPatchID()]
                    );

                tetPolyPatchInterpolation tppiBPatch
                (
                    refCast<const faceTetPolyPatch>(motionUbPatch.patch())
                );

                motionUbPatch == 
                    tppiBPatch.pointToPointInterpolate
                    (
                        totDisplacementB/DB().deltaT().value()
                    );
            }
        }
        else if (fvMotionSolver)
        {
//        Info << "s1fvMotionSolver " <<  endl;
            pointVectorField& motionU =
                const_cast<pointVectorField&>
                (
                    mesh().objectRegistry::
                    lookupObject<pointVectorField>
                    (
                        "pointMotionU"
                    )
                );

            fixedValuePointPatchVectorField& motionUaPatch =
                refCast<fixedValuePointPatchVectorField>
                (
                    motionU.boundaryField()[aPatchID()]
                );

            motionUaPatch ==
                totalDisplacement()/DB().deltaT().value();

            if(twoFluids_)
            {
                const labelList& meshPointsB = 
                    mesh().boundaryMesh()[bPatchID()].meshPoints();

                pointField totDisplacementB =
                    interpolatorAB().pointInterpolate
                    (
                        totalDisplacement()
                    );
        
                forAll (totDisplacementB, pointI)
                {
                    newPoints[meshPointsB[pointI]] -= 
                        totDisplacementB[pointI]; 
                }
            
                fixedValuePointPatchVectorField& motionUbPatch =
                    refCast<fixedValuePointPatchVectorField>
                    (
                        motionU.boundaryField()[bPatchID()]
                    );
//        Info << "s13 " <<  endl;

                motionUbPatch == 
                    totDisplacementB/DB().deltaT().value();
            }
        }

        twoDPointCorrector twoDPointCorr(mesh());
        twoDPointCorr.correctPoints(newPoints);
        mesh().movePoints(newPoints);
        deleteDemandDrivenData(totalDisplacementPtr_);
//    }

        mesh().update();

        aMesh().movePoints(mesh().points());

        correctContactLinePointNormals();

//        if (correctPointNormals_)
//        {
//            correctPointNormals();
//        }

//        if (correctCurvature_)
//        {
//            correctCurvature();
//        }


        // Move correctedFvPatchField fvSubMeshes

        forAll(U().boundaryField(), patchI)
        {
            if
            (
                (
                    U().boundaryField()[patchI].type()
                 == fixedGradientCorrectedFvPatchField<vector>::typeName
                )
                ||
                (
                    U().boundaryField()[patchI].type()
                 == fixedValueCorrectedFvPatchField<vector>::typeName
                )
                ||
                (
                    U().boundaryField()[patchI].type()
                 == zeroGradientCorrectedFvPatchField<vector>::typeName
                )
            )
            {
                correctedFvPatchField<vector>& aU =
                    refCast<correctedFvPatchField<vector> >
                    (
                        U().boundaryField()[patchI]
                    );

                aU.movePatchSubMesh();
            }
        }

        forAll(p().boundaryField(), patchI)
        {
            if
            (
                (
                    p().boundaryField()[patchI].type()
                 == fixedGradientCorrectedFvPatchField<scalar>::typeName
                )
                ||
                (
                    p().boundaryField()[patchI].type()
                 == fixedValueCorrectedFvPatchField<scalar>::typeName
                )
                ||
                (
                    p().boundaryField()[patchI].type()
                 == zeroGradientCorrectedFvPatchField<scalar>::typeName
                )
            )
            {
                correctedFvPatchField<scalar>& aP =
                    refCast<correctedFvPatchField<scalar> >
                    (
                        p().boundaryField()[patchI]
                    );

                aP.movePatchSubMesh();
            }
        }
    }

    return true;
}
*/
/*
bool freeSurface::moveMeshPointsForOldFreeSurfDisplacement()
{
    if(totalDisplacementPtr_)
    {
        pointField newPoints = mesh().points();

        const labelList& meshPointsA = 
            mesh().boundaryMesh()[aPatchID()].meshPoints();

        forAll (totalDisplacement(), pointI)
        {
            newPoints[meshPointsA[pointI]] -= totalDisplacement()[pointI]; 
        }


        // Check mesh motion solver type 
        bool feMotionSolver = 
            mesh().objectRegistry::foundObject<tetPointVectorField>
            (
                "motionU"
            );
        bool fvMotionSolver =
            mesh().objectRegistry::foundObject<pointVectorField>
            (
                "pointMotionU"
            );

        if (feMotionSolver)
        {
            tetPointVectorField& motionU =
                const_cast<tetPointVectorField&>
                (
                    mesh().objectRegistry::
                    lookupObject<tetPointVectorField>
                    (
                        "motionU"
                    )
                );

            fixedValueTetPolyPatchVectorField& motionUaPatch =
                refCast<fixedValueTetPolyPatchVectorField>
                (
                    motionU.boundaryField()[aPatchID()]
                );

            tetPolyPatchInterpolation tppiAPatch
            (
                refCast<const faceTetPolyPatch>
                (
                    motionUaPatch.patch()
                )
            );

            motionUaPatch ==
                tppiAPatch.pointToPointInterpolate
                (
                    totalDisplacement()/DB().deltaT().value()
                );

            if(twoFluids_)
            {
                const labelList& meshPointsB = 
                    mesh().boundaryMesh()[bPatchID()].meshPoints();

                pointField totDisplacementB =
                    interpolatorAB().pointInterpolate
                    (
                        totalDisplacement()
                    );
        
                forAll (totDisplacementB, pointI)
                {
                    newPoints[meshPointsB[pointI]] -= 
                        totDisplacementB[pointI]; 
                }
            
                fixedValueTetPolyPatchVectorField& motionUbPatch =
                    refCast<fixedValueTetPolyPatchVectorField>
                    (
                        motionU.boundaryField()[bPatchID()]
                    );

                tetPolyPatchInterpolation tppiBPatch
                (
                    refCast<const faceTetPolyPatch>(motionUbPatch.patch())
                );

                motionUbPatch == 
                    tppiBPatch.pointToPointInterpolate
                    (
                        totDisplacementB/DB().deltaT().value()
                    );
            }
        }
        else if (fvMotionSolver)
        {
            pointVectorField& motionU =
                const_cast<pointVectorField&>
                (
                    mesh().objectRegistry::
                    lookupObject<pointVectorField>
                    (
                        "pointMotionU"
                    )
                );

            fixedValuePointPatchVectorField& motionUaPatch =
                refCast<fixedValuePointPatchVectorField>
                (
                    motionU.boundaryField()[aPatchID()]
                );

            motionUaPatch ==
                totalDisplacement()/DB().deltaT().value();

            if(twoFluids_)
            {
                const labelList& meshPointsB = 
                    mesh().boundaryMesh()[bPatchID()].meshPoints();

                pointField totDisplacementB =
                    interpolatorAB().pointInterpolate
                    (
                        totalDisplacement()
                    );
        
                forAll (totDisplacementB, pointI)
                {
                    newPoints[meshPointsB[pointI]] -= 
                        totDisplacementB[pointI]; 
                }
            
                fixedValuePointPatchVectorField& motionUbPatch =
                    refCast<fixedValuePointPatchVectorField>
                    (
                        motionU.boundaryField()[bPatchID()]
                    );

                motionUbPatch == 
                    totDisplacementB/DB().deltaT().value();
            }
        }

        twoDPointCorrector twoDPointCorr(mesh());

        twoDPointCorr.correctPoints(newPoints);

        mesh().movePoints(newPoints);

        deleteDemandDrivenData(totalDisplacementPtr_);
    }
        Info << "s15e " <<  endl;
    mesh().update();

    aMesh().movePoints(mesh().points());

//RG-TQ
    correctContactLinePointNormals();
//    Info << "pointNormals (2) " << aMesh().pointAreaNormals() <<endl;



    // Move correctedFvPatchField fvSubMeshes

    forAll(U().boundaryField(), patchI)
    {
        if
        (
            (
                U().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<vector>::typeName
            )
           ||
            (
                U().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<vector>::typeName
            )
           ||
            (
                U().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<vector>::typeName
            )
        )
        {
            correctedFvPatchField<vector>& aU =
                refCast<correctedFvPatchField<vector> >
                (
                    U().boundaryField()[patchI]
                );
            
            aU.movePatchSubMesh();
        }
    }

    forAll(p().boundaryField(), patchI)
    {
        if
        (
            (
                p().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<scalar>::typeName
            )
           ||
            (
                p().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<scalar>::typeName
            )
           ||
            (
                p().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<scalar>::typeName
            )
        )
        {
            correctedFvPatchField<scalar>& aP =
                refCast<correctedFvPatchField<scalar> >
                (
                    p().boundaryField()[patchI]
                );
            
            aP.movePatchSubMesh();
        }
    }

    return true;
}
*/

void freeSurface::updateBoundaryConditions()
{
    ///RG-TQ updateMassFlux
    updateMassFlux();
    updateVelocity();
    updateSurfactantConcentration();
    updatePressure();
    ///RG-TQ updateTemperature();
    updateTemperature();
    ///RG-TQ concentration();
//    updateCon();
}

///RG-TQ updateMassFlux
void freeSurface::updateMassFlux()
{
    if(twoFluids())
    {
        // Compute the absolute pressure on both sides of the interface

        vectorField nB = mesh().boundary()[bPatchID()].nf();

        scalarField DnA = mesh().boundary()[aPatchID()].deltaCoeffs();

        scalarField DnB = interpolatorBA().faceInterpolate
        (
            mesh().boundary()[bPatchID()].deltaCoeffs()
        );

        scalarField pB =
            interpolatorBA().faceInterpolate
            (
                p().boundaryField()[bPatchID()]
            );

//	scalarField rhoB = ntB*conB1*molarMass1().value()+ntB*(1-conB1)*molarMassD().value();

        pB += rhoFluidB().value()*
                  (g_.value() & mesh().C().boundaryField()[aPatchID()]);

//        pB += pCorr().value();
        pB += pCorr_.value();
      	Info<< min(pB) << " < pB < " << max(pB) << endl;
//      Info<< "total pAbs " << pB << endl;

        scalarField pA = pB;        

        const scalarField& K = aMesh().faceCurvatures().internalField();

        if(cleanInterface())
        {   
            pA -= cleanInterfaceSurfTension().value()*K;
        }
        else
        {
            scalarField surfTensionK =
                surfaceTension().internalField()*K;
                
// RG: We should not be subtracting off the average, should we?
// There should be a more elegant (and general) way to write this
            pA -= surfTensionK - gAverage(surfTensionK);
        }

        updateNGradUn();
        pA -= 2.0*(muFluidA().value() - muFluidB().value())*nGradUn();
//             *fac::div(Us())().internalField();

        scalarField TPA = 
            T().boundaryField()[aPatchID()].patchInternalField();
//          Info<< "T().boundaryField()[aPatchID()] " <<  TPA << endl;

        scalarField TPB = interpolatorBA().faceInterpolate
        (
            T().boundaryField()[bPatchID()].patchInternalField()
        );
        scalarField beta =
                 kFluidA().value()*DnA + kFluidB().value()*DnB + VSMALL;

        scalarField T0 = (kFluidA().value()*DnA*TPA 
                        + kFluidB().value()*DnB*TPB)/beta;

        TFs = T0;
        scalarField TSat = T0;

	scalarField KA = kFluidA().value()*DnA/latentHeat().value()/molarMass1().value();
	scalarField KB = kFluidB().value()*DnB/latentHeat().value()/molarMass1().value();
	scalarField KAB = KA + KB;
//	scalarField aX = DfFluidB().value()*DnB*pB/(gasConstantD().value()*TPB*TPB);

        // Compute saturation temperature at the interface
        scalarField conB1b = interpolatorBA().faceInterpolate
        (
            con().boundaryField()[bPatchID()].patchInternalField()
        );
		
	scalarField gradcon1n =interpolatorBA().faceInterpolate
        ( 
//      nB & gradcon1
	con().boundaryField()[bPatchID()].snGrad()
	);

        scalarField conB1 = conB1b + gradcon1n/DnB;
        Info<< min(conB1) << " < conB1 < " << max(conB1) << endl;
//        Info<< " < conB1: " << conB1 << endl;
//        scalarField TSat = antoineB().value()/(antoineA().value()
//                   -log10((pB*conB)/133.322368))-antoineC().value()+273.15;
				   
	if (min(conB1)< 0.0)
		{
		FatalErrorIn("freeSurface::updateMassflux()")
		<< "concentration field of vapor in the air has a minimum value that is negative."
		<< abort(FatalError);
		}

/*
       	scalarField ncB1b = interpolatorBA().faceInterpolate
        (
            nc1().boundaryField()[bPatchID()].patchInternalField()
        );
*/
//        scalarField ntB = pB/gasConstant1().value()/molarMass1().value()/TPB;
//        scalarField ntB = pB/gasConstant1().value()/molarMass1().value()/TFs;
//        ncB1b = ntB*conB;

	dimensionedScalar pdimi("pdim",dimPressure,1.0);
    	dimensionedScalar ntBv = pTotal().value()*pdimi/gasConstant1()/molarMass1()/((TLeft() + TRight())/2.0);

	volScalarField ntt = (1.0-fluidIndicator())*((p()+pCorr())/gasConstant1()/molarMass1()/T())+fluidIndicator()*nt();
	dimensionedScalar ntAvg = fvc::domainIntegrate((1.0-fluidIndicator())*ntt)/fvc::domainIntegrate((1.0-fluidIndicator()));

        scalarField ntB = interpolatorBA().faceInterpolate
        (
            nt().boundaryField()[bPatchID()].patchInternalField()
        );


    	ntB = scalarField(aMesh().nFaces(), ntAvg.value());
//	scalarField ntB = pB/gasConstant1().value()/molarMass1().value()/TFs;

//       	scalarField ncB1b = conB1*ntB;
//        Info<< min(ncB1b) << " < ncB1b < " << max(ncB1b) << endl;
//        ncB1 = ncB1b;	
//      scalarField ntB = ncB1/conB;
	Info<< min(ntB) << " < ntB < " << max(ntB) << endl;
//	ntB = gAverage(ntB);
//      ncBDb = ncB1b*(1.0-conB)/conB;

//        Info<< min(ncB1*molarMass1().value()) << " < rhoB1 = ncB1*molarMass1 < " << max(ncB1*molarMass1().value()) << endl;

        // Compute temperature at the interface

/*
        scalarField rhoB = interpolatorBA().faceInterpolate
        ( 
            rho().boundaryField()[bPatchID()].patchInternalField()
        );


        scalarField rhoBD = interpolatorBA().faceInterpolate
        ( 
            (rho() - rho1())().boundaryField()[bPatchID()].patchInternalField()
        );
*/

/*
        rhoB1 = interpolatorBA().faceInterpolate
        ( 
            rho1().boundaryField()[bPatchID()].patchInternalField()
        );

        scalarField rhoBD = pB*(1.0-conB)/(gasConstantD().value()*TFs);
        scalarField rhoB = rhoBD + rhoB1;
        Info<< min(rhoB) << " < rhoB < " << max(rhoB) << endl;
*/

//	07/04/2012 accomodation coefficient
    	dimensionedScalar preFactor = 2.0*(acmCoef()/(2.0-acmCoef()))/sqrt(6.2831852*gasConstant1());
	scalarField alpha = preFactor.value()*latentHeat().value()*sqrt(TFs)/(TSat*TSat);
	scalarField gamma = preFactor.value()*(pA-pB)*sqrt(TFs)/rhoFluidA().value()/TSat;

	dimensionedScalar alpha1 = preFactor*latentHeat();
	scalarField gamma1 = preFactor.value()*(pA-pB)/rhoFluidA().value();

    	dimensionedScalar gasC = gasConstant1()*molarMass1();
    	scalarField TgasC = sqrt(TFs)*TSat*TSat;

	scalarField aJ = DnB*DfFluidB().value()+Vn;	
//	scalarField aR = aB.value()/(aA.value()-log10(gCD*TFs/133.322368))/(aA.value()-log10(gCD*TFs/133.322368))/log(10.0);
	scalarField aR = antoineB().value()
	/(antoineA().value()-log10(ntB*conB1*gasC.value()*TFs/133.322368))
	/(antoineA().value()-log10(ntB*conB1*gasC.value()*TFs/133.322368))
	/log(10.0);

	scalarField a23 = 1.0-conB1;
	scalarField a32 = 0.5*conB1 * ntB * (3.0 * alpha1.value() * TFs - alpha1.value() * TSat + gamma1 * TSat)/TgasC;
	scalarField a34 = ((-alpha1.value() + gamma1) * TSat + alpha1.value() * TFs) * ntB *TFs / TgasC;
	scalarField a35 = -ntB * conB1 * (2 * alpha1.value() * TFs - alpha1.value() * TSat + gamma1 * TSat)*TFs/TSat/TgasC;
	scalarField a52 = aR/TFs;
	scalarField a54 = aR/conB1;

/*
    volVectorField gradnc1v = fvc::grad(nc1());
	vectorField gradnc1 = gradnc1v.boundaryField()[bPatchID()].patchInternalField();
	scalarField gradnc1n = interpolatorBA().faceInterpolate
        (
            (nB & gradnc1)/ntB
	);


    volVectorField gradnc1v = fvc::grad(con());
	vectorField gradnc1 = gradnc1v.boundaryField()[bPatchID()].patchInternalField();
	scalarField gradnc1n = interpolatorBA().faceInterpolate
        (
            (nB & gradnc1)
	);
		
        Info<< min(gradnc1n) << " < gradCon1N before newton iteration < " << max(gradnc1n) << endl;
*/
//	Vn = 0.0*Vn;
//	J = 0.0*J;

/*
        volVectorField gradrho1v = fvc::grad(rho1());
	vectorField gradrho1 = gradrho1v.boundaryField()[bPatchID()].patchInternalField();
	scalarField gradrho1n = interpolatorBA().faceInterpolate
        (
            nB & gradrho1
	);
		
        Info<< min(gradrho1n) << " < gradrho1N < " << max(gradrho1n) << endl;

        const volScalarField rhoBDv = p()*(1.0-con())/(gasConstantD()*T());
//        volVectorField gradrhobdv = fvc::grad(rhoBDv);

        volVectorField gradrhov = fvc::grad(rho1())+fvc::grad(rhoBDv);
	vectorField gradrho = gradrhov.boundaryField()[bPatchID()].patchInternalField();
	scalarField gradrhon = interpolatorBA().faceInterpolate
        (
            nB & gradrho
	);
        Info<< min(gradrhon) << " < gradrhoN < " << max(gradrhon) << endl;

        scalar mD1 = molarMassD().value()/molarMass1().value();
//	Vn = (mD1/rhoBD)*DfFluidB().value()*gradrho1n;
	Vn = (-DfFluidB().value() * gradrhon + J) / (rhoB);
	J = DfFluidB().value()*gradrho1n + Vn * rhoB1;
	const scalarField& Sf = aMesh().S();
	scalar JCorr = gSum(J*Sf)/gSum(Sf);

        Info<< min(Vn) << " < Vn < " << max(Vn) << endl;
*/
/*
        volVectorField gradconv = fvc::grad(con());
	vectorField gradCon = gradconv.boundaryField()[bPatchID()].patchInternalField();
	scalarField gradConN =interpolatorBA().faceInterpolate
        ( 
            nB & gradCon
	);
*/	
        dimensionedScalar dTFs("dTFs",dimless,1.0);
        dimensionedScalar dTSat("dTSat",dimless,1.0);
//        dimensionedScalar drhoB1("drhoB1",dimless,1.0);
        dimensionedScalar dVn("dVn",dimless,1.0);

    	int corr=0;
	double Normold = 1.0;
	double Normnew = 1.0;
	double r1 = 1.0;
	double r2 = 1.0;
	double r3 = 1.0;
	double r4 = 1.0;
	double r5 = 1.0;
//	double r6 = 1.0;

//        while(dTFs.value()>1e-6||dVn.value()>1e-6||drhoB1.value()>1e-6)
//        while(dTFs.value()>1e-6||drhoB1.value()>1e-6)
        while(Normnew>1e-8&&corr++<20)
        {
//            if(corr++<100)       
//            {
        	scalarField Jold = J;
        	scalarField TFsold = TFs;
        	scalarField Vnold = Vn;
//        	scalarField rhoBDold = rhoBD;
//        	scalarField rhoB1old = rhoB1;
//        	scalarField gradrho1nold = gradrho1n;
//        	scalarField gradnc1nold = gradnc1n;
        	scalarField conB1old = conB1;
        	scalarField TSatold = TSat;
/*
		alpha = preFactor.value()*latentHeat().value()*sqrt(TFs)/(TSat*TSat);
		gamma = preFactor.value()*(pA-pB)*sqrt(TFs)/rhoFluidA().value()/TSat;
		scalarField GbT = gamma - alpha*TSat + alpha*TFs;			
		ntB = pB/gasConstant1().value()/molarMass1().value()/TFs;
		dimensionedScalar aA = antoineA();
		dimensionedScalar aB = antoineB();
		scalarField aJ = (DnB*DfFluidB().value()+Vn)*ntB;			
//		scalarField aR = aB.value()/(aA.value()-log10(gCD*TFs/133.322368))/(aA.value()-log10(gCD*TFs/133.322368))/log(10.0);
		scalarField aR = aB.value()
		/(aA.value()-log10(ntB*conB1*gasConstant1().value()*molarMass1().value()*TFs/133.322368))
		/(aA.value()-log10(ntB*conB1*gasConstant1().value()*molarMass1().value()*TFs/133.322368))
		/log(10.0);
*/
//		ntB = pB/gasConstant1().value()/molarMass1().value()/TFs;
		aJ = DnB*DfFluidB().value()+Vn;			
		aR = antoineB().value()
	/(antoineA().value()-log10(ntB*conB1*gasC.value()*TFs/133.322368))
	/(antoineA().value()-log10(ntB*conB1*gasC.value()*TFs/133.322368))
	/log(10.0);

    		TgasC = sqrt(TFs)*TSat*TSat;

		alpha = alpha1.value()*sqrt(TFs)/(TSat*TSat);
		gamma = gamma1*sqrt(TFs)/TSat;

		a23 = 1.0-conB1;
		a32 = 0.5*conB1 * ntB * (3.0 * alpha1.value() * TFs - alpha1.value() * TSat + gamma1 * TSat)/TgasC;
		a34 = ((-alpha1.value() + gamma1) * TSat + alpha1.value() * TFs) * ntB *TFs / TgasC;
		a35 = -ntB * conB1 * (2 * alpha1.value() * TFs - alpha1.value() * TSat + gamma1 * TSat)*TFs/TSat/TgasC;
		a52 = aR/TFs;
		a54 = aR/conB1;		
/*
	       	scalarField g1 = -J + DfFluidB().value()*gradrho1n + Vn*rhoB1;
//        	scalarField g2 = -mD1*DfFluidB().value()*gradrho1n + Vn*rhoBD;
			scalarField g2 = aX*(TFs-TPB)-mD1*DfFluidB().value()*gradrho1n + Vn*rhoBD;

        	scalarField g3 = -J + alpha*(TFs - TSat)*rhoB1 + gamma*rhoB1;
        	scalarField g4 = KA*(TPA - TFs) - KB*(TFs-TPB) - J;

        	scalarField g5 = -pB + (rhoBD*gasConstantD().value() + rhoB1*gasConstant1().value())*TFs;

        	scalarField g6 = antoineB().value()/(antoineA().value()
		           -log10(rhoB1*gasConstant1().value()*TFs/133.322368))-(antoineC().value()-273.15) - TSat;
*/
	       	scalarField g1 = -J + ntB*DfFluidB().value()*(conB1-conB1b)*DnB + Vn*ntB*conB1;
		scalarField g2 = - DfFluidB().value()*(conB1-conB1b)*DnB + Vn*(1.0-conB1);

        	scalarField g3 = -J + alpha*(TFs - TSat)*ntB*conB1 + gamma*ntB*conB1;
        	scalarField g4 = KA*(TPA - TFs) - KB*(TFs-TPB) - J;

//        	scalarField g5 = -pB + (rhoBD*gasConstantD().value() + rhoB1*gasConstant1().value())*TFs;

        	scalarField g5 = antoineB().value()/(antoineA().value()
		           -log10(ntB*conB1*gasC.value()*TFs/133.322368))-(antoineC().value()-273.15) - TSat;
//        	Info<<"---------------g1-g5---------------" << endl;
/*
		r1 = max( sign(g1)*(g1))/(max(sign(J)*(J))+VSMALL );
//		r2 = max( sign(g2)*(g2))/(max(sign(Vn*rhoBD)*(Vn*rhoBD))+VSMALL );
		r2 = max( sign(g2)*(g2))/(max(sign(Vn*ncBD)*(Vn*ncBD))+VSMALL );
		r3 = max( sign(g3)*(g3))/(max(sign(J)*(J))+VSMALL );
		r4 = max( sign(g4)*(g4))/(max(sign(J*molarMass1().value())*(J*molarMass1().value()))+VSMALL );
//		r5 = max( max(g5/pB), -min(g5/pB));
		r5 = max( max(g5/TSat), -min(g5/TSat));

        	Normold = sqrt(r1*r1 + r2*r2 + r3*r3 + r4*r4+ r5*r5);
        	//Info<<"original_maxNorm," << Normold << endl;
*/

		scalarField h1 = (-KAB * a23 * (a35 * a54 + a34) * g1 + KAB * conB1 * ntB * (a35 * a54 + a34) * g2 + KAB * aJ * ntB * (a23 + conB1) * g3 + ntB * aJ * (a35 * a52 + a32) * (a23 + conB1) * g4 + KAB * a35 * aJ * ntB * (a23 + conB1) * g5) / ((a23 * a35 * a54 - a23 * aJ * ntB - aJ * conB1 * ntB + a23 * a34) * KAB - a23 * a35 * a52 * aJ * ntB - a35 * a52 * aJ * conB1 * ntB - a23 * a32 * aJ * ntB - a32 * aJ * conB1 * ntB);

		scalarField h2 = (a23 * (a35 * a54 + a34) * g1 - conB1 * ntB * (a35 * a54 + a34) * g2 - ntB * aJ * (a23 + conB1) * g3 + (-a23 * a35 * a54 + a23 * aJ * ntB + aJ * conB1 * ntB - a23 * a34) * g4 - a35 * aJ * ntB * (a23 + conB1) * g5) / ((a23 * a35 * a54 - a23 * aJ * ntB - aJ * conB1 * ntB + a23 * a34) * KAB - a23 * a35 * a52 * aJ * ntB - a35 * a52 * aJ * conB1 * ntB - a23 * a32 * aJ * ntB - a32 * aJ * conB1 * ntB);

		scalarField h3 = (-aJ * (a35 * a52 + KAB + a32) * g1 + (-a35 * a52 * aJ * ntB + KAB * a35 * a54 - KAB * aJ * ntB - a32 * aJ * ntB + KAB * a34) * g2 + KAB * aJ * g3 + aJ * (a35 * a52 + a32) * g4 + KAB * a35 * aJ * g5) / ((a23 * a35 * a54 - a23 * aJ * ntB - aJ * conB1 * ntB + a23 * a34) * KAB - a23 * a35 * a52 * aJ * ntB - a35 * a52 * aJ * conB1 * ntB - a23 * a32 * aJ * ntB - a32 * aJ * conB1 * ntB);

		scalarField h4 = (-a23 * (a35 * a52 + KAB + a32) * g1 + conB1 * ntB * (a35 * a52 + KAB + a32) * g2 + KAB * a23 * g3 + a23 * (a35 * a52 + a32) * g4 + KAB * a23 * a35 * g5) / ((a23 * a35 * a54 - a23 * aJ * ntB - aJ * conB1 * ntB + a23 * a34) * KAB - a23 * a35 * a52 * aJ * ntB - a35 * a52 * aJ * conB1 * ntB - a23 * a32 * aJ * ntB - a32 * aJ * conB1 * ntB);

		scalarField h5 = (-a23 * (KAB * a54 + a32 * a54 - a34 * a52) * g1 + conB1 * ntB * (KAB * a54 + a32 * a54 - a34 * a52) * g2 + (-a23 * a52 * aJ * ntB - a52 * aJ * conB1 * ntB + KAB * a23 * a54) * g3 + (a23 * a52 * aJ * ntB + a52 * aJ * conB1 * ntB + a23 * a32 * a54 - a23 * a34 * a52) * g4 + (KAB * a23 * aJ * ntB + KAB * aJ * conB1 * ntB + a23 * a32 * aJ * ntB + a32 * aJ * conB1 * ntB - KAB * a23 * a34) * g5) / ((a23 * a35 * a54 - a23 * aJ * ntB - aJ * conB1 * ntB + a23 * a34) * KAB - a23 * a35 * a52 * aJ * ntB - a35 * a52 * aJ * conB1 * ntB - a23 * a32 * aJ * ntB - a32 * aJ * conB1 * ntB);

/*
		scalarField h1 = (-TFs * KAB * ntB * (conB1 - 1.0) * (-aR * alpha + GbT) * g1 - KAB * TFs * conB1 * ntB * (-aR * alpha + GbT) * g2 - KAB * TFs * aJ * g3 - aJ * alpha * conB1 * ntB * (TFs - aR) * g4 + KAB * TFs * aJ * alpha * conB1 * g5 * ntB) / ((conB1 * ntB - ntB) * TFs * KAB * GbT + (-aR * alpha * conB1 * ntB + aR * alpha * ntB + aJ) * TFs * KAB + TFs * aJ * alpha * conB1 * ntB - aJ * aR * alpha * conB1 * ntB);

		scalarField h2 = (TFs * ntB * (conB1 - 1) * (-aR * alpha + GbT) * g1 + TFs * conB1 * ntB * (-aR * alpha + GbT) * g2 + TFs * aJ * g3 - (-aR * alpha * conB1 * ntB + GbT * conB1 * ntB + aR * alpha * ntB - GbT * ntB + aJ) * TFs * g4 - TFs * aJ * alpha * conB1 * ntB * g5) / ((-aR * alpha * conB1 * ntB + GbT * conB1 * ntB + aR * alpha * ntB - GbT * ntB + aJ) * TFs * KAB + TFs * aJ * alpha * conB1 * ntB - aJ * aR * alpha * conB1 * ntB);

		scalarField h3 = (aJ * (TFs * alpha * conB1 * ntB - aR * alpha * conB1 * ntB + KAB * TFs) * g1 + (KAB * TFs * aR * alpha * ntB + TFs * aJ * alpha * conB1 * ntB - aJ * aR * alpha * conB1 * ntB - GbT * KAB * TFs * ntB + KAB * TFs * aJ) * g2 - KAB * TFs * aJ * g3 - aJ * alpha * conB1 * ntB * (TFs - aR) * g4 + KAB * TFs * aJ * alpha * conB1 * g5 * ntB) / (ntB * (conB1 * ntB - ntB) * TFs * KAB * GbT + ntB * (-aR * alpha * conB1 * ntB + aR * alpha * ntB + aJ) * TFs * KAB + aJ * alpha * conB1 * ntB * ntB * TFs - ntB * ntB * aJ * aR * alpha * conB1);

		scalarField h4 = (-(conB1 - 1) * (TFs * alpha * conB1 * ntB - aR * alpha * conB1 * ntB + KAB * TFs) * g1 - conB1 * (TFs * alpha * conB1 * ntB - aR * alpha * conB1 * ntB + KAB * TFs) * g2 + KAB * TFs * (conB1 - 1) * g3 + alpha * ntB * conB1 * (conB1 - 1) * (TFs - aR) * g4 - KAB * TFs * alpha * conB1 * ntB * (conB1 - 1) * g5) / ((-aR * alpha * conB1 * ntB + GbT * conB1 * ntB + aR * alpha * ntB - GbT * ntB + aJ) * TFs * KAB + aJ * alpha * conB1 * ntB * TFs - aJ * aR * alpha * conB1 * ntB);

		scalarField h5 = (aR * (conB1 - 1) * (-TFs * alpha * conB1 * ntB + GbT * conB1 * ntB - KAB * TFs) * g1 + aR * conB1 * (-TFs * alpha * conB1 * ntB + GbT * conB1 * ntB - KAB * TFs) * g2 + aR * (KAB * TFs * conB1 - KAB * TFs + aJ * conB1) * g3 - aR * conB1 * (-TFs * alpha * conB1 * ntB + GbT * conB1 * ntB + TFs * alpha * ntB - GbT * ntB + aJ) * g4 - TFs * conB1 * (GbT * KAB * conB1 * ntB + aJ * alpha * conB1 * ntB - GbT * KAB * ntB + KAB * aJ) * g5) / ((-aR * alpha * conB1 * ntB + GbT * conB1 * ntB + aR * alpha * ntB - GbT * ntB + aJ) * conB1 * TFs * KAB + aJ * alpha * conB1 * conB1 * ntB * TFs - aJ * aR * alpha * conB1 * conB1 * ntB);
*/

/*
		scalarField h1 = (-KAB * TFs * g3 - alpha * ntB * conB1 * (TFs - aR) * g4 + KAB * TFs * alpha * conB1 * ntB * g5) / (alpha * ntB * conB1 * TFs - aR * alpha * conB1 * ntB + KAB * TFs);

		scalarField h2 = (-TFs * alpha * conB1 * g5 * ntB + TFs * g3 - TFs * g4) / (alpha * ntB * conB1 * TFs - aR * alpha * conB1 * ntB + KAB * TFs);
		
		scalarField h5 = (aR * g3 - aR * g4 - TFs * (alpha * ntB * conB1 + KAB) * g5) / (alpha * ntB * conB1 * TFs - aR * alpha * conB1 * ntB + TFs * KAB);
*/

/* OLD 1-5

		scalarField h1 = (-KAB * TFs * g3 - alpha * ncB1 * (TFs - aR) * g4 + KAB * TFs * alpha * g5 * ncB1) / (alpha * ncB1 * molarMass1().value() * TFs - aR * alpha * ncB1 * molarMass1().value() + KAB * TFs);

		scalarField h2 = (-alpha * ncB1 * molarMass1().value() * TFs * g5 + TFs * molarMass1().value() * g3 - TFs * g4) / (alpha * ncB1 * molarMass1().value() * TFs - aR * alpha * ncB1 * molarMass1().value() + KAB * TFs);

		scalarField h5 = (aR * molarMass1().value() * g3 - aR * g4 - TFs * (alpha * ncB1 * molarMass1().value() + KAB) * g5) / (alpha * ncB1 * molarMass1().value() * TFs - aR * alpha * ncB1 * molarMass1().value() + KAB * TFs);
//        	Info<<"---------------h1-h5---------------" << endl;
*/
		J = -h1 + Jold;
		TFs = -h2 + TFsold;
		Vn = -h3 + Vnold;
//		rhoBD = -h4 + rhoBDold;
//		gradrho1n = -h5 + gradrho1nold;
//		gradnc1n = -h4 + gradnc1nold;
		conB1 = -h4 + conB1old;
		TSat = -h5 + TSatold;

//		Vn = J/ntB;
//		conB1 = (DnB*conB1b*ntB*DfFluidB().value()+J)/(DnB*ntB*DfFluidB().value()+J);

/*
		Vn = J/ntB;
		gradnc1n = J * ncBD / ntB / DfFluidB().value() / ntB;
		scalarField gradntBn = pB*DnB*(TPB - TFs)/(gasConstant1().value()*molarMass1().value()*TFs*TFs);
		gradn1B = ntB*gradnc1n + conB*gradntBn;
*/
/*
		conB = conBb + gradnc1n/DnB;
		ntB = pB/gasConstant1().value()/molarMass1().value()/TFs;
		ncB1 = ntB*conB;
		ncBD = ntB - ncB1;
*/
/*
		ntB = pB/gasConstant1().value()/molarMass1().value()/TFs;
		ncB1 = ncB1b + gradn1B/DnB;
		ncBD = ntB - ncB1;
		conB = ncB1/ntB;

  		Info<< " ntB*gradnc1n + ncB1*gradntBn <<< " << (gradn1B) << endl;
  		Info<< " conB1 <<< " << (conB) << endl;
        	Info<<"---------------r1-r5---------------" << endl;
		Info<< min(ncB1) << " < ncB1 < " << max(ncB1) << endl;
		Info<< min(gradn1B) << " < gradn1B < " << max(gradn1B) << endl;
*/

/*
		alpha = preFactor.value()*latentHeat().value()*sqrt(TFs)/(TSat*TSat);
		gamma = preFactor.value()*(pA-pB)*sqrt(TFs)/rhoFluidA().value()/TSat;

	       	g1 = -J + DfFluidB().value()*gradrho1n + Vn*rhoB1;
//        	g2 = -mD1*DfFluidB().value()*gradrho1n + Vn*rhoBD;
//        	g2 = DfFluidB().value()*(gradrhon -gradrho1n) + Vn*rhoBD;
		g2 = aX*(TFs-TPB)-mD1*DfFluidB().value()*gradrho1n + Vn*rhoBD;
        	g3 = -J + alpha*(TFs - TSat)*rhoB1 + gamma*rhoB1;
        	g4 = KA*(TPA - TFs) - KB*(TFs-TPB) - J;
        	g5 = -pB + (rhoBD*gasConstantD().value() + rhoB1*gasConstant1().value())*TFs;

        	g6 = antoineB().value()/(antoineA().value()
		           -log10(rhoB1*gasConstant1().value()*TFs/133.322368))-(antoineC().value()-273.15) - TSat;

		r1 = max( sign(g1)*(g1))/(max(sign(J)*(J))+VSMALL );
		r2 = max( sign(g2)*(g2))/(max(sign(Vn*rhoBD)*(Vn*rhoBD))+VSMALL );
		r3 = max( sign(g3)*(g3))/(max(sign(J)*(J))+VSMALL );
		r4 = max( sign(g4)*(g4))/(max(sign(J)*(J))+VSMALL );
		r5 = max( max(g5/pB), -min(g5/pB));
		r6 = max( max(g6/TSat), -min(g6/TSat));
*/
		r1 = max( sign(g1)*(g1))/(max(sign(J)*(J))+VSMALL );
//		r2 = max( sign(g2)*(g2))/(max(sign(Vn*rhoBD)*(Vn*rhoBD))+VSMALL );
//		r2 = max( sign(g2)*(g2))/(max(sign(Vn*ntB*(1.0-conB1))*(Vn*ntB*(1.0-conB1)))+VSMALL );
		r2 = max( sign(g2)*(g2))/(max(sign(Vn*(1.0-conB1))*(Vn*(1.0-conB1)))+VSMALL );
		r3 = max( sign(g3)*(g3))/(max(sign(J)*(J))+VSMALL );
		r4 = max( sign(g4)*(g4))/(max(sign(J)*(J))+VSMALL );
//		r5 = max( max(g5/pB), -min(g5/pB));
		r5 = max( max(g5/TSat), -min(g5/TSat));

//        	Normnew = sqrt(r3*r3 + r4*r4+ r5*r5);
        	Normnew = sqrt(r1*r1 + r2*r2 + r3*r3 + r4*r4+ r5*r5);

        	Info<<"original_maxNorm," << Normold <<" ----- new_maxNorm," << Normnew << endl;

        	dTFs = max(mag(TFs-TFsold))/(min(TFs)+VSMALL);
//        	drhoB1 = max(mag(rhoB1-rhoB1old))/(min(rhoB1)+VSMALL);
        	dVn = max(mag(Vn-Vnold))/(min(Vn)+VSMALL);
			
                if(min(TFs)<0)
		{
	         	FatalErrorIn("freeSurface::updateMassflux()")
	                << "Negative interfacial temperature "
	                    << abort(FatalError);
		}
/*
            }
            else
            {
            	FatalErrorIn("freeSurface::updateMassflux()")
                << "Iterations for interfacial temperature failed to converge "

                    << abort(FatalError);
            }
*/

        }
//		JCorr = gSum(J*Sf)/gSum(Sf);
//		J = J - JCorr;

//*new
//	Vn = J/ntB;
//	conB1 = (DnB*conB1b*ntB*DfFluidB().value()+J)/(DnB*ntB*DfFluidB().value()+J);
//	gradcon1n = (conB1 - conB1b) * DnB;
//	scalarField gradntBn = pB*DnB*(TPB - TFs)/(gasConstant1().value()*molarMass1().value()*TFs*TFs);
//////////////////*/

//	gradn1B = ntB*gradnc1n;// + conB*gradntBn;

/*5.0
	Vn = J/ntB;
	gradnc1n = J * ncBD / ntB / DfFluidB().value() / ntB;
//////////////////*/

	ncB1 = ntB*conB1;
	ncBD = ntB*(1.0-conB1);
	J = J*molarMass1().value();
        massFlux().internalField() = J;

//        rhoB1 = ntBv.value()*conB1*molarMass1().value();
//        scalarField rhoBD = pB*(1.0-conB)/(gasConstantD().value()*TFs);

//        scalar mD1 = molarMassD().value()/molarMass1().value();

/*	1.0
	scalarField aX = 0.0*TPB;
	Vn = (J * mD1 - aX * TFs + aX * TPB) / (mD1 * rhoB1 + rhoBD);
	scalarField gradrho1n = (TFs * aX * rhoB1 - TPB * aX * rhoB1 + J * rhoBD) / DfFluidB().value() / (mD1 * rhoB1 + rhoBD);
	gradnc1n = (gradrho1n/molarMass1().value())/ntBv.value();

/////////////////*/
/*/	2.0
	scalarField aX = DfFluidB().value()*DnB*pB/(gasConstantD().value()*TPB*TPB);
	Vn = (J * mD1 - aX * TFs + aX * TPB) / (mD1 * rhoB1 + rhoBD);
	scalarField gradrho1n = (TFs * aX * rhoB1 - TPB * aX * rhoB1 + J * rhoBD) / DfFluidB().value() / (mD1 * rhoB1 + rhoBD);
	gradnc1n = (gradrho1n/molarMass1().value())/ntBv.value();
//////////////////
//	3.0
	scalarField aX = 0.0*TPB;
	Vn = (J * mD1 - aX * TFs + aX * TPB) / (mD1 * rhoB1 + rhoBD);
	scalarField gradrho1n = (TFs * aX * rhoB1 - TPB * aX * rhoB1 + J * rhoBD) / DfFluidB().value() / (mD1 * rhoB1 + rhoBD);
	gradnc1n = (gradrho1n/molarMass1().value()-conB*gradntBn)/ntBv.value();
/////////////////*/
/*/	4.0
	scalarField aX = DfFluidB().value()*DnB*pB/(gasConstantD().value()*TPB*TPB);
	Vn = (J * mD1 - aX * TFs + aX * TPB) / (mD1 * rhoB1 + rhoBD);
	scalarField gradrho1n = (TFs * aX * rhoB1 - TPB * aX * rhoB1 + J * rhoBD) / DfFluidB().value() / (mD1 * rhoB1 + rhoBD);
	gradnc1n = (gradrho1n/molarMass1().value()-conB*gradntBn)/ntBv.value();

/////////////////*/
//        Info<< "T().boundaryField()[aPatchID()] " <<  (TPA - TFs) << endl;
	heatFluxl().internalField() = kFluidA().value()*DnA*(TPA - TFs);
	heatFluxv().internalField() = kFluidB().value()*DnB*(TFs - TPB);
	conci().internalField() = conB1;
        dTi().internalField() = TFs - TSat;
        Info<<"interface iteration# = " << corr <<" dTFs = " << dTFs.value()<< " dVn = " << dVn.value()  <<" r1= " << r1 << "; r2= " << r2 <<" r3= " << r3 << "; r4= " << r4<< "; r5= " << r5<<endl;

//	Info<< min(TSat) << " < TSat < " << max(TSat) << endl;	
//  	Info<< min(TFs) << " < TFs < " << max(TFs) << endl;
        Info<< min(TFs-293) << " < TFs < " << max(TFs-293) <<endl;
        Info<< min(TSat-293) << " < TSat < " << max(TSat-293) <<endl;
	Info<< min(J) << " < J < " << max(J) << endl;
	Info<< min(Vn) << " < Vn < " << max(Vn) << endl;
  	Info<< min(conB1) << " < conB1 < " << max(conB1) << endl;
//        Info<< min((gradnc1n)) << " < gradnc1n before correction < " << max((gradnc1n)) <<endl;
/*
        Info<< " conB1--- " << conB <<endl;
        Info<< " J--- " << J <<endl;
        Info<< " TSat--- " << TSat-293 <<endl;
        Info<< " TFs--- " << TFs -293 <<endl;
*/
/*/	Info<< min(rhoB1) << " < rhoB1 < " << max(rhoB1) << endl;
//  	Info<< " TFs <<< " << (TFs) << endl;
  	Info<< " J <<< " << (J) << endl;
  	Info<< " Vn <<< " << (Vn) << endl;
  	Info<< " gradnc1n <<< " << (gradnc1n) << endl;
  	Info<< " ntB*gradnc1n + ncB1*gradntBn <<< " << (gradn1B) << endl;
*/
/*
        if
        (
	nc1().boundaryField()[bPatchID()].type()
         == fixedGradientFvPatchField<scalar>::typeName
        )
        {
            fixedGradientFvPatchField<scalar>& nc1B =
                refCast<fixedGradientFvPatchField<scalar> >
                (
                    nc1().boundaryField()[bPatchID()]
                );
//  	Info<< " gradnc1n <<< " << nc1B.gradient()<< endl;
//	nc1B.gradient() = interpolatorAB().faceInterpolate(gradnc1n);
        nc1B.gradient() = interpolatorAB().faceInterpolate(gradn1B);
//  	Info<< " gradnc1n2 <<< " << nc1B.gradient()<< endl;
        }

        else
        {
            FatalErrorIn("freeSurface::updatemass()")
                << "Bounary condition on " << nc1().name() 
                    <<  " for freeSurfaceShadow patch is " 
                    << nc1().boundaryField()[bPatchID()].type() 
                    << ", instead of " 
                    << fixedGradientFvPatchField<scalar>::typeName 
                    << abort(FatalError);
        }
*/
/*
        if
        (
	con().boundaryField()[bPatchID()].type()
         == fixedGradientFvPatchField<scalar>::typeName
        )
        {
            fixedGradientFvPatchField<scalar>& conB =
                refCast<fixedGradientFvPatchField<scalar> >
                (
                    con().boundaryField()[bPatchID()]
                );
	conB.gradient() = interpolatorAB().faceInterpolate(gradcon1n);
        }

        else
        {
            FatalErrorIn("freeSurface::updatemass()")
                << "Bounary condition on " << con().name() 
                    <<  " for freeSurfaceShadow patch is " 
                    << con().boundaryField()[bPatchID()].type() 
                    << ", instead of " 
                    << fixedGradientFvPatchField<scalar>::typeName 
                    << abort(FatalError);
        }

        if
        (
	rho1().boundaryField()[bPatchID()].type()
         == fixedGradientFvPatchField<scalar>::typeName
        )
        {
            fixedGradientFvPatchField<scalar>& rho1B =
                refCast<fixedGradientFvPatchField<scalar> >
                (
                    rho1().boundaryField()[bPatchID()]
                );

            rho1B.gradient() = interpolatorAB().faceInterpolate(gradrho1n);

        }

        else
        {
            FatalErrorIn("freeSurface::updatemass()")
                << "Bounary condition on " << rho1().name() 
                    <<  " for freeSurfaceShadow patch is " 
                    << rho1().boundaryField()[bPatchID()].type() 
                    << ", instead of " 
                    << fixedGradientFvPatchField<scalar>::typeName 
                    << abort(FatalError);
        }
*/

//	fixed value BC on air side for component 1
	if
        (
            con().boundaryField()[bPatchID()].type()
         == fixedValueFvPatchField<scalar>::typeName
        )
        {
            con().boundaryField()[bPatchID()] == 
		interpolatorAB().faceInterpolate(conB1);           
        }

        else
        {
            FatalErrorIn("freeSurface::updatemass()")
                << "Bounary condition on " << con().name() 
                    <<  " for freeSurfaceShadow patch is " 
                    << con().boundaryField()[bPatchID()].type() 
                    << ", instead " 
                    << fixedValueFvPatchField<scalar>::typeName 
                    << abort(FatalError);
        }

    }
}


void freeSurface::updateVelocity()
{
    if(twoFluids())
    {  	
        vectorField nA = mesh().boundary()[aPatchID()].nf();

        vectorField nB = mesh().boundary()[bPatchID()].nf();

        scalarField DnB = interpolatorBA().faceInterpolate
        (
            mesh().boundary()[bPatchID()].deltaCoeffs()
        );

        scalarField DnA = mesh().boundary()[aPatchID()].deltaCoeffs();

        vectorField UtPA = 

            U().boundaryField()[aPatchID()].patchInternalField();
//        gradTi().internalField() = UtPA;
        Ui().internalField() = UtPA;
        scalarField rhoB = rho().boundaryField()[bPatchID()].patchInternalField();

        if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientCorrectedFvPatchField<vector>::typeName
        )
        {
            fixedGradientCorrectedFvPatchField<vector>& aU =
                refCast<fixedGradientCorrectedFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );

            UtPA += aU.corrVecGrad();

        }

        UtPA -= nA*(nA & UtPA);
        Ui().internalField() = UtPA;

        vectorField UtPB = interpolatorBA().faceInterpolate
        (
            U().boundaryField()[bPatchID()].patchInternalField()
        );

        vectorField UPB = interpolatorBA().faceInterpolate
        (
            U().boundaryField()[bPatchID()].patchInternalField()
        );

        if
        (
            U().boundaryField()[bPatchID()].type()
         == fixedValueCorrectedFvPatchField<vector>::typeName
        )
        {
            fixedValueCorrectedFvPatchField<vector>& bU =
                refCast<fixedValueCorrectedFvPatchField<vector> >
                (
                    U().boundaryField()[bPatchID()]
                );

            UtPB += interpolatorBA().faceInterpolate(bU.corrVecGrad());
        }

        UtPB -= nA*(nA & UtPB);

        vectorField UtFs = 
            muFluidA().value()*DnA*UtPA 
          + muFluidB().value()*DnB*UtPB;

        vectorField UnFs = 
            nA*(phi_.boundaryField()[aPatchID()]
           /mesh().boundary()[aPatchID()].magSf());
//		-J/rhoFluidA().value());

        Us().internalField() += UnFs - nA*(nA&Us().internalField());
        correctUsBoundaryConditions();

        UtFs -= (muFluidA().value() - muFluidB().value())*
            (fac::grad(Us())&aMesh().faceAreaNormals())().internalField();

        vectorField tangentialSurfaceTensionForce(nA.size(), vector::zero);
/*
        wordList patchFieldTypes
        (
            aMesh().boundary().size(),
            zeroGradientFaPatchVectorField::typeName
        );

        forAll(aMesh().boundary(), patchI)
        {
            label ngbPolyPatchID = 
                aMesh().boundary()[patchI].ngbPolyPatchIndex();

            if (ngbPolyPatchID != -1)
            {
                if
                (
                    T().boundaryField()[ngbPolyPatchID].type()
                 == fixedValueFvPatchScalarField::typeName
                )
                {
                    patchFieldTypes[patchI] =
                        fixedValueFaPatchScalarField::typeName;
                }
            }
        }

        areaScalarField Ti
        (
            IOobject
            (
                "Ti",
                DB().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            aMesh(),
            dimensioned<scalar>("Ti", dimTemperature, 0),
            patchFieldTypes
        );
*/
        Ti().internalField() = TFs;

        // Set fixedValue boundary conditions
/*
        forAll(Ti().boundaryField(), patchI)
        {

//	Info << "Ti().boundaryField()[patchI]: " << Ti().boundaryField()[patchI] << endl;

            if
            (
                Ti().boundaryField()[patchI].type()
             == fixedValueFaPatchScalarField::typeName
            )
            {
                label ngbPolyPatchID =
                    aMesh().boundary()[patchI].ngbPolyPatchIndex();
                label ngbPolyPatchStart = 
                    mesh().boundaryMesh()[ngbPolyPatchID].start();
                labelList ngbPolyPatchFaces =
                    aMesh().boundary()[patchI].ngbPolyPatchFaces() 
                  - ngbPolyPatchStart;

		Info << "ngbPolyPatchID: " << ngbPolyPatchID << endl;
		Info << "ngbPolyPatchStart: " << ngbPolyPatchStart << endl;
		Info << "ngbPolyPatchFaces: " << ngbPolyPatchFaces << endl;
		Info << "ngbPolyPatchName: " << mesh().boundaryMesh()[ngbPolyPatchID].name() << endl;

                Ti.boundaryField()[patchI] ==
                    scalarField
                    (
                        T().boundaryField()[ngbPolyPatchID],
                        ngbPolyPatchFaces
                    );


		if
		(
		    mesh().boundaryMesh()[ngbPolyPatchID].name()
		=="leftWall"
		)
		{
		
		    scalarField TiFluid = Ti().boundaryField()[patchI].patchInternalField();

		    Info <<"-----------------------------"<<endl<< "TiFluid on the left wall: " << TiFluid << endl;

		    scalarField DnF = aMesh().boundary()[patchI].deltaCoeffs();

		    Info << "DnFi on the left wall: " << DnF << endl;

//		    dimensionedScalar theta = gAverage (contactAngle().boundaryField()[patchI]);

		    Info << "theta on the left wall: " << ctangle().value() << endl;

		    dimensionedScalar kFluidEff = (ctangle()*kFluidA() + (180.0-ctangle())*kFluidB())/180.0;

		    Info << "conductivityEff on the left wall: " << kFluidEff.value() << endl;

		    scalarField sigmai = kFluidEff.value()*DnF+kWall().value()/thWall().value();

		    Info << "sigmai on the left: " << sigmai << endl;
		    
		    scalarField TiLeftEff = ((kWall().value()*TLeft().value()/thWall().value())+kFluidEff.value()*TiFluid*DnF)/sigmai;

		    Info << "TiLeftEff on the left wall: " << TiLeftEff << endl<<"-----------------------------"<<endl;
		    Ti().boundaryField()[patchI] == TiLeftEff;

		}

		if
		(
		    mesh().boundaryMesh()[ngbPolyPatchID].name()
		=="rightWall"
		)
		{
		
		    scalarField TiFluid = Ti().boundaryField()[patchI].patchInternalField();

		    Info << "TiFluid on the right wall: " << TiFluid << endl;

		    scalarField DnF = aMesh().boundary()[patchI].deltaCoeffs();
		    Info << "DnFi on the right wall: " << DnF << endl;

//		    dimensionedScalar theta = gAverage (contactAngle().boundaryField()[patchI]);

		    Info << "theta on the right wall: " << ctangle().value() << endl;

		    dimensionedScalar kFluidEff = (ctangle()*kFluidA() + (180.0-ctangle())*kFluidB())/180.0;

		    Info << "conductivityEff on the right wall: " << kFluidEff.value() << endl;

		    scalarField sigmai = kFluidEff.value()*DnF+kWall().value()/thWall().value();

		    Info << "sigmai on the right: " << sigmai << endl;
		    
		    scalarField TiRightEff = ((kWall().value()*TRight().value()/thWall().value())+kFluidEff.value()*TiFluid*DnF)/sigmai;

		    Info << "TiRightEff on the right wall: " << TiRightEff << endl<<"-----------------------------"<<endl;

		    Ti().boundaryField()[patchI] == TiRightEff;

		}

            }
        }
*/
//	Info<< "Ti boundary: " << Ti().boundaryField() << endl;


        Ti().correctBoundaryConditions();

	dimensionedScalar TFsAvg("TFsAvg", dimTemperature, gAverage(TFs));
	TiDiff()=Ti()-TFsAvg;
//	Info<< "TiDiff: " << TiDiff() << endl;

//        vectorField gradTFs = fac::grad(Ti())().internalField();
	vectorField gradTFs = grad(Ti())().internalField();
        gradTi().internalField() = gradTFs;

//        vectorField tangradTFs = gradTFs-nA*(nA & gradTFs);
//	Info<< "controlpoints: " << aMesh().areaCentres().internalField() << endl;
//	Info<< "gradT " << gradTFs << endl;
//	Info<< "tangradT " << tangradTFs << endl;
//	Info<< "magtangradTFs " << mag(tangradTFs) << endl;
//        gradTFs.write();
//        tangradTFs.write();

/*
	tangentialSurfaceTensionForce +=
                          tempCoeffSurfTension().value()*gradTFs;
          Info<< "gradT " << gradTFs << endl;
          Info<< "tangentialSurfaceTensionForce " << tangentialSurfaceTensionForce << endl;

*/

//        volVectorField gradT = fvc::grad(T());
//	vectorField gradTFs = 
//            gradT.boundaryField()[aPatchID()].patchInternalField();
//	gradTFs -= nA*(nA & gradTFs);

        ///RG-TQ  Add thermocapillary force


	tangentialSurfaceTensionForce +=
                          tempCoeffSurfTension().value()*gradTFs;

	tangentialSurfaceTensionForce -= nA*(nA&tangentialSurfaceTensionForce);
//	Info<< "thermocapillarity " << tangentialSurfaceTensionForce << endl;
	
//wrong	tangentialSurfaceTensionForce += (J*J/rhoB1)*(UPB - (UPB&nB)*nB)/mag(UPB - (UPB&nB)*nB);
//	tangentialSurfaceTensionForce += (J*J/rhoB1)*UtPB/(mag(UPB)+VSMALL);
	tangentialSurfaceTensionForce += (J*J/(ncB1*molarMass1().value()))*UtPB/(mag(UPB)+VSMALL);
//	Info<< "tangentialSurfaceTension " << tangentialSurfaceTensionForce << endl;	

/*
        if(!cleanInterface())
        {
            tangentialSurfaceTensionForce = 
                surfaceTensionGrad()().internalField();
        }
        else
        {
            vectorField surfaceTensionForce =
                cleanInterfaceSurfTension().value()
               *fac::edgeIntegrate
                (
                    aMesh().Le()*aMesh().edgeLengthCorrection()
                )().internalField();

            tangentialSurfaceTensionForce =
                surfaceTensionForce 
              - cleanInterfaceSurfTension().value()
               *aMesh().faceCurvatures().internalField()*nA;
        }


        ///RG-TQ Use fvc::grad to incorporate proper boundary conditions for T
        volVectorField gradT = fvc::grad(T());
	vectorField gradTFs = 
            gradT.boundaryField()[aPatchID()].patchInternalField();
	gradTFs -= nA*(nA & gradTFs);

        //  Info<< "gradT " << gradTFs << endl;

        ///RG-TQ  Add thermocapillary force
//	tangentialSurfaceTensionForce +=
//                          tempCoeffSurfTension().value()*gradTFs;
*/

        UtFs += tangentialSurfaceTensionForce;

        UtFs /= muFluidA().value()*DnA + muFluidB().value()*DnB + VSMALL;

//	scalarField qadv = rho().boundaryField()[aPatchID()].patchInternalField() * mag(UtFs) * CpFluidA().value() * TFs;
//	Info<< "Heat flux advective by UtFs: " << qadv << endl;
//	scalarField qadve = mag(qadv - (qadv&nA) );
//        TiDiff().internalField() = qadv;
        Us().internalField() = UnFs + UtFs;
        correctUsBoundaryConditions();

        // Store old-time velocity field U()
        U().oldTime();

        U().boundaryField()[bPatchID()] == 
            interpolatorAB().faceInterpolate(UtFs)
          + nB*fvc::meshPhi(rho(),U())().boundaryField()[bPatchID()]/
            mesh().boundary()[bPatchID()].magSf()
          + interpolatorAB().faceInterpolate(nA*Vn);

//        Info<< "UtFs " << UtFs << endl;
//        Info<< "U().boundaryField()[bPatchID()] " << U().boundaryField()[bPatchID()]<< endl;

        if
        (
            p().boundaryField()[bPatchID()].type()
         == fixedGradientFvPatchField<scalar>::typeName
        )
        {
            fixedGradientFvPatchField<scalar>& pB =
                refCast<fixedGradientFvPatchField<scalar> >
                (
                    p().boundaryField()[bPatchID()]
                );

            pB.gradient() = 
               - rhoFluidB().value()
                *(
                     nB&fvc::ddt(U())().boundaryField()[bPatchID()]
                 );

//	    dimensionedScalar TRefB = fvc::domainIntegrate((1.0-fluidIndicator())*T())/fvc::domainIntegrate(1.0 - fluidIndicator());
//	    vectorField buoForceB = rhoFluidB().value()*(interpolatorAB().faceInterpolate(TFs)-TRefB.value())*g_.value()/interpolatorAB().faceInterpolate(TFs);
//          pB.gradient() -= nB & buoForceB;

//	    dimensionedScalar rhoRefB = fvc::domainIntegrate((1.0-fluidIndicator())*rho())/fvc::domainIntegrate(1.0 - fluidIndicator());
//	    vectorField buoForceB = (rhoB - rhoRefB.value())*g().value();
//            pB.gradient() += nB & buoForceB;

        }

        // Update fixedGradient boundary condition on patch A
	//0917-2012 update normal gradient of velocity
        updateNGradUn();

        vectorField nGradU =
            muFluidB().value()*(UtPB - UtFs)*DnB
          + tangentialSurfaceTensionForce
          + muFluidA().value()*nA*nGradUn()
//          - muFluidA().value()*nA*fac::div(Us())().internalField()
          + (muFluidB().value() - muFluidA().value())
           *(fac::grad(Us())().internalField()&nA);

        nGradU /= muFluidA().value() + VSMALL;

        
        if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientCorrectedFvPatchField<vector>::typeName
        )
        {
            fixedGradientCorrectedFvPatchField<vector>& aU =
                refCast<fixedGradientCorrectedFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );

            aU.gradient() = nGradU;
        }
        else if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientFvPatchField<vector>::typeName
        )
        {
            fixedGradientFvPatchField<vector>& aU =
                refCast<fixedGradientFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );
            
            aU.gradient() = nGradU;
        }
        else
        {
            FatalErrorIn("freeSurface::updateVelocity()")
                << "Bounary condition on " << U().name() 
                    <<  " for freeSurface patch is " 
                    << U().boundaryField()[aPatchID()].type() 
                    << ", instead " 
                    << fixedGradientCorrectedFvPatchField<vector>::typeName 
                    << " or "
                    << fixedGradientFvPatchField<vector>::typeName 
                    << abort(FatalError);
        }
    }
    else
    {
        vectorField nA = aMesh().faceAreaNormals().internalField();

        vectorField UnFs =
            nA*phi_.boundaryField()[aPatchID()]
           /mesh().boundary()[aPatchID()].magSf();

        // Correct normal component of free-surface velocity
        Us().internalField() += UnFs - nA*(nA&Us().internalField());
        correctUsBoundaryConditions();

        vectorField tangentialSurfaceTensionForce(nA.size(), vector::zero);

        if(!cleanInterface())
        {
             tangentialSurfaceTensionForce = 
                 surfaceTensionGrad()().internalField();
        }
        else
        {
            vectorField surfaceTensionForce =
                cleanInterfaceSurfTension().value()
               *fac::edgeIntegrate
                (
                    aMesh().Le()*aMesh().edgeLengthCorrection()
                )().internalField();

            tangentialSurfaceTensionForce =
                surfaceTensionForce 
              - cleanInterfaceSurfTension().value()
               *aMesh().faceCurvatures().internalField()*nA;

            if (muFluidA().value() < SMALL)
            {
                tangentialSurfaceTensionForce = vector::zero;
            }
        }

        vectorField tnGradU =
            tangentialSurfaceTensionForce/(muFluidA().value() + VSMALL)
          - (fac::grad(Us())&aMesh().faceAreaNormals())().internalField();

        vectorField UtPA =
            U().boundaryField()[aPatchID()].patchInternalField();
        UtPA -= nA*(nA & UtPA);

        scalarField DnA = mesh().boundary()[aPatchID()].deltaCoeffs();

        vectorField UtFs = UtPA + tnGradU/DnA;
        
        Us().internalField() = UtFs + UnFs;
        correctUsBoundaryConditions();

        vectorField nGradU =
            tangentialSurfaceTensionForce/(muFluidA().value() + VSMALL)
          - nA*fac::div(Us())().internalField()
          - (fac::grad(Us())().internalField()&nA);

        if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientCorrectedFvPatchField<vector>::typeName
        )
        {
            fixedGradientCorrectedFvPatchField<vector>& aU =
                refCast<fixedGradientCorrectedFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );

            aU.gradient() = nGradU;
        }
        else if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientFvPatchField<vector>::typeName
        )
        {
            fixedGradientFvPatchField<vector>& aU =
                refCast<fixedGradientFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );

            aU.gradient() = nGradU;
        }
        else
        {
            FatalErrorIn("freeSurface::updateVelocity()")
                << "Bounary condition on " << U().name() 
                    <<  " for freeSurface patch is " 
                    << U().boundaryField()[aPatchID()].type() 
                    << ", instead " 
                    << fixedGradientCorrectedFvPatchField<vector>::typeName 
                    << " or "
                    << fixedGradientFvPatchField<vector>::typeName 
                    << abort(FatalError);
        }
    }
}


void freeSurface::updatePressure()
{
    // Correct pressure boundary condition at the free-surface
        
    vectorField nA = mesh().boundary()[aPatchID()].nf();

    if(twoFluids())
    {
        scalarField pA =
            interpolatorBA().faceInterpolate
            (
                p().boundaryField()[bPatchID()]
            );

        const scalarField& K = aMesh().faceCurvatures().internalField();

        Info << "Free surface curvature: min = " << gMin(K)
            << ", max = " << gMax(K)
            << ", average = " << gAverage(K) << endl << flush;
//        TiDiff().internalField() = K;
        if(cleanInterface())
        {   
//             pA -= cleanInterfaceSurfTension().value()*(K - gAverage(K));
            pA -= cleanInterfaceSurfTension().value()*K;
        }
        else
        {
            scalarField surfTensionK =
                surfaceTension().internalField()*K;
                
            pA -= surfTensionK - gAverage(surfTensionK);
        }

        pA -= 2.0*(muFluidA().value() - muFluidB().value())*nGradUn();
//0917-2012             *fac::div(Us())().internalField();

        vector R0 = vector::zero;
//         vector R0 = gAverage(mesh().C().boundaryField()[aPatchID()]);

        pA -= (rhoFluidA().value() - rhoFluidB().value())*
            (
                g_.value()
              & (
                    mesh().C().boundaryField()[aPatchID()]
                  - R0
                )
            );

        p().boundaryField()[aPatchID()] == pA;
    }
    else
    {
        vector R0 = vector::zero;
//         vector R0 = gAverage(mesh().C().boundaryField()[aPatchID()]);

        scalarField pA =
          - rhoFluidA().value()*
            (
                g_.value()
              & (
                    mesh().C().boundaryField()[aPatchID()]
                  - R0
                )
            );

        const scalarField& K = aMesh().faceCurvatures().internalField();
        
        Info << "Free surface curvature: min = " << gMin(K)
            << ", max = " << gMax(K) << ", average = " << gAverage(K) 
            << endl;
        

        if(cleanInterface())
        {
//             pA -= cleanInterfaceSurfTension().value()*(K - gAverage(K));
            pA -= cleanInterfaceSurfTension().value()*K;
        }
        else
        {
            scalarField surfTensionK =
                surfaceTension().internalField()*K;
            
            pA -= surfTensionK - gAverage(surfTensionK);
        }

        pA -= 2.0*muFluidA().value()*fac::div(Us())().internalField();

        p().boundaryField()[aPatchID()] == pA;
    }


    // Set modified pressure at patches with fixed apsolute
    // pressure

    vector R0 = vector::zero;
//     vector R0 = gAverage(mesh().C().boundaryField()[aPatchID()]);

    for (int patchI=0; patchI < p().boundaryField().size(); patchI++)
    {
        if 
        (
            p().boundaryField()[patchI].type()
         == fixedValueFvPatchScalarField::typeName
        )
        {
            if (patchI != aPatchID())
            {
                p().boundaryField()[patchI] ==
//                     p().boundaryField()[patchI]
                  - rho().boundaryField()[patchI]*
                    (g_.value()&(mesh().C().boundaryField()[patchI] - R0));
            }
        }
    }
}

///RG-TQ update Temperature
void freeSurface::updateTemperature()
{
    if(twoFluids())
    {
	// Update fixedValue boundary condition on patch B
        if
        (
            T().boundaryField()[bPatchID()].type()
         == fixedValueFvPatchField<scalar>::typeName
        )
        {
            T().boundaryField()[bPatchID()] == 
              interpolatorAB().faceInterpolate(TFs);
        }
        else
        {
            FatalErrorIn("freeSurface::updateTemperature()")
                << "Bounary condition on " << T().name() 
                    <<  " for freeSurfaceShadow patch is " 
                    << T().boundaryField()[bPatchID()].type() 
                    << ", instead " 
                    << fixedValueFvPatchField<scalar>::typeName 
                    << abort(FatalError);
        }

        // Update fixedGradient boundary condition on patch A

        scalarField DnB = interpolatorBA().faceInterpolate
        (
            mesh().boundary()[bPatchID()].deltaCoeffs()
        );

        scalarField TPB = interpolatorBA().faceInterpolate
        (
            T().boundaryField()[bPatchID()].patchInternalField()
        );

        scalarField nGradT = kFluidB().value()*(TPB - TFs)*DnB
                           - latentHeat().value()*J;

        nGradT /= kFluidA().value() + VSMALL;

        if
        (
            T().boundaryField()[aPatchID()].type()
         == fixedGradientFvPatchField<scalar>::typeName
        )
        {
            fixedGradientFvPatchField<scalar>& aT =
                refCast<fixedGradientFvPatchField<scalar> >
                (
                    T().boundaryField()[aPatchID()]
                );            
            aT.gradient() = nGradT;
        }


        else
        {
            FatalErrorIn("freeSurface::updateTemperature()")

                << "Bounary condition on " << T().name() 
                    <<  " for freeSurface patch is " 
                    << T().boundaryField()[aPatchID()].type() 
                    << ", instead " 
                    << fixedGradientFvPatchField<scalar>::typeName 

                    << abort(FatalError);
        }
	


	/*
        // Update fixedValue boundary condition on patch A
        if
        (
            T().boundaryField()[aPatchID()].type()
         == fixedValueFvPatchField<scalar>::typeName
        )
        {
            T().boundaryField()[aPatchID()] == TFs;
//              interpolatorAB().faceInterpolate(TFs);
        }
        else
        {
            FatalErrorIn("freeSurface::updateTemperature()")
                << "Bounary condition on " << T().name() 
                    <<  " for freeSurfaceShadow patch is " 
                    << T().boundaryField()[aPatchID()].type() 
                    << ", instead " 
                    << fixedValueFvPatchField<scalar>::typeName 
                    << abort(FatalError);
        }

        // Update fixedGradient boundary condition on patch B

        scalarField DnA = mesh().boundary()[aPatchID()].deltaCoeffs();

        scalarField TPA = T().boundaryField()[aPatchID()].patchInternalField();

        scalarField TPB = interpolatorBA().faceInterpolate
        (
            T().boundaryField()[bPatchID()].patchInternalField()
        );

//        scalarField nGradT = kFluidB().value()*(TPB - TFs)*DnA
//                           - latentHeat().value()*J;

//        nGradT /= kFluidA().value() + VSMALL;

        scalarField nGradT =   kFluidA().value()*(TPA - TFs)*DnA
                           - latentHeat().value()*J;

        nGradT /= kFluidB().value() + VSMALL;
        Info << "bT.gradient() " << nGradT<<endl;
        if
        (
            T().boundaryField()[bPatchID()].type()
         == fixedGradientFvPatchField<scalar>::typeName
        )
        {
            fixedGradientFvPatchField<scalar>& bT =
                refCast<fixedGradientFvPatchField<scalar> >
                (
                    T().boundaryField()[bPatchID()]
                );
            
            bT.gradient() = interpolatorAB().faceInterpolate(nGradT);
        Info << "bT.gradient() " << bT.gradient()<<endl;
        Info << "bT.gradient() calculated " << (TPB-TFs)*DnA<<endl;
        }
        else
        {
            FatalErrorIn("freeSurface::updateTemperature()")
                << "Bounary condition on " << T().name() 
                    <<  " for freeSurface patch is " 
                    << T().boundaryField()[bPatchID()].type() 
                    << ", instead " 
                    << fixedGradientFvPatchField<scalar>::typeName 
                    << abort(FatalError);
        }
	*/

    }

    /*
    else
    {
//        scalarField nGradT = 0;
//          - (fvc::grad(T())().internalField()&nA);


        if
        (
            T().boundaryField()[aPatchID()].type()
         != fixedGradientFvPatchField<scalar>::typeName
        )
        {
            fixedGradientFvPatchField<scalar>& aT =
                refCast<fixedGradientFvPatchField<scalar> >
                (
                    T().boundaryField()[aPatchID()]
                );

            aT.gradient() = nGradT;
        }
        else 
        {
            FatalErrorIn("freeSurface::updateTemperature()")
                << "Bounary condition on " << T().name() 
                    <<  " for freeSurface patch is " 
                    << T().boundaryField()[aPatchID()].type() 
                    << ", instead " 
                    << fixedGradientFvPatchField<scalar>::typeName 
                    << abort(FatalError);
        }
    }
    */

}

/*
///RG-TQ update Concentration
void freeSurface::updateCon()
{
    if(twoFluids())
    {
        // Update fixedValue boundary condition on patch A

        scalarField conB = interpolatorBA().faceInterpolate
        (
            con().boundaryField()[bPatchID()].patchInternalField()
        );

        if
        (
            con().boundaryField()[aPatchID()].type()
         == fixedValueFvPatchField<scalar>::typeName
        )
        {
            con().boundaryField()[aPatchID()] == conB;           
        }
        else
        {
            FatalErrorIn("freeSurface::updatecon()")
                << "Bounary condition on " << con().name() 
                    <<  " for freeSurfaceShadow patch is " 
                    << con().boundaryField()[aPatchID()].type() 
                    << ", instead " 
                    << fixedValueFvPatchField<scalar>::typeName 
                    << abort(FatalError);
        }

        // Update fixedGradient boundary condition on patch B

        scalarField nGradcon = interpolatorAB().faceInterpolate(J/rhoVapor)/DfFluidB().value();
//        Info << "nGradcon" << nGradcon << endl;
        
        if
        (
            con().boundaryField()[bPatchID()].type()
         == fixedGradientFvPatchField<scalar>::typeName
        )
        {
            fixedGradientFvPatchField<scalar>& bcon =
                refCast<fixedGradientFvPatchField<scalar> >
                (
                    con().boundaryField()[bPatchID()]
                );
            bcon.gradient() = nGradcon;
        }
        else
        {
            FatalErrorIn("freeSurface::updatecon()")
                << "Bounary condition on " << con().name() 
                    <<  " for freeSurface patch is " 
                    << con().boundaryField()[bPatchID()].type() 
                    << ", instead " 
                    << fixedGradientFvPatchField<scalar>::typeName 
                    << abort(FatalError);
        }
    }
    else
    {
 
//        scalarField nGradT = 0;
//          - (fvc::grad(T())().internalField()&nA);

        if
        (
            con().boundaryField()[bPatchID()].type()
         != fixedGradientFvPatchField<scalar>::typeName
        )
       {
            fixedGradientFvPatchField<scalar>& aT =
                refCast<fixedGradientFvPatchField<scalar> >
                (
                    T().boundaryField()[aPatchID()]
                );

            aT.gradient() = nGradT;

        }
        else 
        {
            FatalErrorIn("freeSurface::updatecon()")
                << "Bounary condition on " << con().name() 
                    <<  " for freeSurface patch is " 
                    << con().boundaryField()[aPatchID()].type() 
                    << ", instead " 
                    << fixedGradientFvPatchField<scalar>::typeName 
                    << abort(FatalError);
        }

    }
}
*/

void freeSurface::updateSurfaceFlux()
{
    Phis() = fac::interpolate(Us()) & aMesh().Le();
}


void freeSurface::updateSurfactantConcentration()
{
    if(!cleanInterface())
    {
        Info << "Correct surfactant concentration" << endl << flush;
        
        updateSurfaceFlux();        

        // Crate and solve the surfactanta transport equation
        faScalarMatrix CsEqn
        (
            fam::ddt(surfactantConcentration())
          + fam::div(Phis(), surfactantConcentration())
          - fam::laplacian
            (
                surfactant().surfactDiffusion(),
                surfactantConcentration()
            )
        );


        if(surfactant().soluble())
        {
            const scalarField& C =
                mesh().boundary()[aPatchID()]
               .lookupPatchField<volScalarField, scalar>("C");

            areaScalarField Cb
            (
                IOobject
                (
                    "Cb",
                    DB().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                aMesh(),
                dimensioned<scalar>("Cb", dimMoles/dimVolume, 0),
                zeroGradientFaPatchScalarField::typeName
            );

            Cb.internalField() = C;
            Cb.correctBoundaryConditions();

            CsEqn += 
                fam::Sp
                (
                    surfactant().surfactAdsorptionCoeff()*Cb
                  + surfactant().surfactAdsorptionCoeff()
                   *surfactant().surfactDesorptionCoeff(),
                    surfactantConcentration()
                )
              - surfactant().surfactAdsorptionCoeff()
               *Cb*surfactant().surfactSaturatedConc();
        }

        CsEqn.solve();

        Info << "Correct surface tension" << endl;

        surfaceTension() =
            cleanInterfaceSurfTension()
          + surfactant().surfactR()
           *surfactant().surfactT()
           *surfactant().surfactSaturatedConc()
           *log(1.0 - surfactantConcentration()
           /surfactant().surfactSaturatedConc());
        
        if(neg(min(surfaceTension().internalField())))
        {
            FatalErrorIn
            (
                "void freeSurface::correctSurfactantConcentration()"
            ) 
                << "Surface tension has negative value" 
                    << abort(FatalError);
        }
    }
}


void freeSurface::correctUsBoundaryConditions()
{
    forAll(Us().boundaryField(), patchI)
    {
        if
        (
            Us().boundaryField()[patchI].type()
         == calculatedFaPatchVectorField::typeName
        )
        {
            vectorField& pUs = Us().boundaryField()[patchI];

            pUs = Us().boundaryField()[patchI].patchInternalField();

            label ngbPolyPatchID = 
                aMesh().boundary()[patchI].ngbPolyPatchIndex();

            if(ngbPolyPatchID != -1)
            {
                if
                (
                    (
                        U().boundaryField()[ngbPolyPatchID].type() 
                     == slipFvPatchVectorField::typeName
                    )
                 ||
                    (
                        U().boundaryField()[ngbPolyPatchID].type() 
                     == symmetryFvPatchVectorField::typeName
                    )
                )
                {
                    vectorField N = 
                        aMesh().boundary()[patchI].ngbPolyPatchFaceNormals();

                    pUs -= N*(N&pUs);

                }
            }
        }
    }

    Us().correctBoundaryConditions();
}


///Transfered from old ZT_version
dimensionedVector freeSurface::geometricMeanPosition() const
{
    return gAverage(mesh().C().boundaryField()[aPatchID()]);
}


vector freeSurface::totalPressureForce() const
{
    const scalarField& S = aMesh().S();

    const vectorField& n = aMesh().faceAreaNormals().internalField();

    const scalarField& P = p().boundaryField()[aPatchID()];

    vectorField pressureForces = S*P*n;

    return gSum(pressureForces);
}


vector freeSurface::totalViscousForce() const
{
    const scalarField& S = aMesh().S();
    const vectorField& n = aMesh().faceAreaNormals().internalField();

    vectorField nGradU =
        U().boundaryField()[aPatchID()].snGrad();

    vectorField viscousForces = 
      - muFluidA().value()*S
       *(
            nGradU
          + (fac::grad(Us())().internalField()&n)
          - (n*fac::div(Us())().internalField())
        );

    return gSum(viscousForces);
}


vector freeSurface::totalSurfaceTensionForce() const
{
    const scalarField& S = aMesh().S();

    const vectorField& n = aMesh().faceAreaNormals().internalField();

    const scalarField& K = aMesh().faceCurvatures().internalField();
//    Info<< "curvature " << K << endl;

    vectorField surfTensionForces(n.size(), vector::zero);

    if(cleanInterface())
    {
        surfTensionForces =
            S*cleanInterfaceSurfTension().value()
           *fac::edgeIntegrate
            (
                aMesh().Le()*aMesh().edgeLengthCorrection()
            )().internalField();
    }
    else
    {
        surfTensionForces *= surfaceTension().internalField()*K;
    }

    return gSum(surfTensionForces);
}

//z0312 freeSurface::initializeControlPointsPo
void freeSurface::initializeControlPointsPosition()
{
    scalarField deltaH = scalarField(aMesh().nFaces(), 0.0);

    pointField displacement = pointDisplacement(deltaH);

    const faceList& faces = aMesh().faces();
    const pointField& points = aMesh().points();


    pointField newPoints = points + displacement;

    scalarField sweptVol(faces.size(), 0.0);

    forAll(faces, faceI)

    {
        sweptVol[faceI] = -faces[faceI].sweptVol(points, newPoints);
    }

    vectorField faceArea(faces.size(), vector::zero);

    forAll (faceArea, faceI)
    {
        faceArea[faceI] = faces[faceI].normal(newPoints);
    }

    forAll(deltaH, faceI)
    {
        deltaH[faceI] = sweptVol[faceI]/
            (faceArea[faceI] & facesDisplacementDir()[faceI]);
    }

    forAll(fixedFreeSurfacePatches_, patchI)
    {
        label fixedPatchID = 
            aMesh().boundary().findPatchID
            (

                fixedFreeSurfacePatches_[patchI]
            );

        if(fixedPatchID == -1)
        {
            FatalErrorIn("freeSurface::freeSurface(...)")
                << "Wrong faPatch name in the fixedFreeSurfacePatches list"
                    << " defined in the freeSurfaceProperties dictionary"
                    << abort(FatalError);
        }

        const labelList& eFaces =
            aMesh().boundary()[fixedPatchID].edgeFaces();

        forAll(eFaces, edgeI)
        {
            deltaH[eFaces[edgeI]] *= 2.0;
        }
    }

    displacement = pointDisplacement(deltaH);
}


//z092812 freeSurface::smoothing()
void freeSurface::smoothing()
{
    controlPoints() = aMesh().areaCentres().internalField();

    scalarField deltaH = scalarField(aMesh().nFaces(), 0.0);

    pointField displacement = pointDisplacement(deltaH);

//    vectorField newMeshPoints = mesh().points();
    vectorField newMeshPoints = mesh().allPoints();

    const labelList& meshPointsA = 
        mesh().boundaryMesh()[aPatchID()].meshPoints();

    forAll (displacement, pointI)
    {
        newMeshPoints[meshPointsA[pointI]] += displacement[pointI];
    }

    if(twoFluids_)
    {
        const labelList& meshPointsB = 
            mesh().boundaryMesh()[bPatchID_].meshPoints();

        pointField displacementB =
            interpolatorAB().pointInterpolate
            (
                displacement
            );
        
        forAll (displacementB, pointI)
        {
            newMeshPoints[meshPointsB[pointI]] += displacementB[pointI]; 
        }
    }

    // Update total displacement field

    if(totalDisplacementPtr_ && (curTimeIndex_ < DB().timeIndex()))
    {
        FatalErrorIn("freeSurface::movePoints()")
            << "Total displacement of free surface points "
                << "from previous time step is not absorbed by the mesh."
                << abort(FatalError);
    }
    else if (curTimeIndex_ < DB().timeIndex())
    {
        totalDisplacement() = displacement;

        curTimeIndex_ = DB().timeIndex();
    }
    else
    {
        totalDisplacement() += displacement;
    }

    twoDPointCorrector twoDPointCorr(mesh());
    twoDPointCorr.correctPoints(newMeshPoints);

    mesh().movePoints(newMeshPoints);

    aMesh().movePoints(mesh().points());

    const scalarField& Sf = aMesh().S();
    const vectorField& Nf = aMesh().faceAreaNormals().internalField();

    deltaH =
       -mesh().phi().boundaryField()[aPatchID()]*DB().deltaT().value()
       /(Sf*(Nf & facesDisplacementDir()));
    
    forAll(fixedFreeSurfacePatches_, patchI)
    {
        label fixedPatchID = 
            aMesh().boundary().findPatchID
            (
                fixedFreeSurfacePatches_[patchI]
            );

        if(fixedPatchID == -1)
        {
            FatalErrorIn("freeSurface::freeSurface(...)")
                << "Wrong faPatch name in the fixedFreeSurfacePatches list"
                    << " defined in the freeSurfaceProperties dictionary"
                    << abort(FatalError);
        }

        const labelList& eFaces =
            aMesh().boundary()[fixedPatchID].edgeFaces();

        forAll(eFaces, edgeI)
        {
            deltaH[eFaces[edgeI]] *= 2.0;
        }
    }

    displacement = pointDisplacement(deltaH);

//    newMeshPoints = mesh().points();
    newMeshPoints = mesh().allPoints();
    forAll (displacement, pointI)
    {
        newMeshPoints[meshPointsA[pointI]] += displacement[pointI];
    }

    if(twoFluids_)
    {
        const labelList& meshPointsB = 
            mesh().boundaryMesh()[bPatchID_].meshPoints();

        pointField displacementB =
            interpolatorAB().pointInterpolate
            (
                displacement
            );
        
        forAll (displacementB, pointI)
        {
            newMeshPoints[meshPointsB[pointI]] += displacementB[pointI]; 
        }
    }

    // Update total displacement field

    if(totalDisplacementPtr_ && (curTimeIndex_ < DB().timeIndex()))
    {
        FatalErrorIn("freeSurface::movePoints()")
            << "Total displacement of free surface points "
                << "from previous time step is not absorbed by the mesh."
                << abort(FatalError);
    }
    else if (curTimeIndex_ < DB().timeIndex())
    {
        totalDisplacement() = displacement;

        curTimeIndex_ = DB().timeIndex();
    }
    else
    {
        totalDisplacement() += displacement;
    }

    twoDPointCorr.correctPoints(newMeshPoints);

    mesh().movePoints(newMeshPoints);

    aMesh().movePoints(mesh().points());

    if (correctPointNormals_)
    {
        correctPointNormals();
    }

    if (correctCurvature_)
    {
//         correctCurvature();
        smoothCurvature();
    }
}
/*
void freeSurface::initializeControlPointsPosition()
{
    scalarField deltaH = scalarField(aMesh().nFaces(), 0.0);

    pointField displacement = pointDisplacement(deltaH);

    const faceList& faces = aMesh().faces();
    const pointField& points = aMesh().points();


    pointField newPoints = points + displacement;

    scalarField sweptVol(faces.size(), 0.0);

    forAll(faces, faceI)
    {
        sweptVol[faceI] = -faces[faceI].sweptVol(points, newPoints);
    }

    vectorField faceArea(faces.size(), vector::zero);

    forAll (faceArea, faceI)
    {
        faceArea[faceI] = faces[faceI].normal(newPoints);
    }

    forAll(deltaH, faceI)
    {
        deltaH[faceI] = sweptVol[faceI]/
            (faceArea[faceI] & facesDisplacementDir()[faceI]);
    }

    displacement = pointDisplacement(deltaH);
}
*/

scalar freeSurface::maxCourantNumber()
{
    scalar CoNum = 0;

    if(cleanInterface())
    {
        const scalarField& dE =aMesh().lPN();

        CoNum = gMax
        (
            DB().deltaT().value()/
            sqrt
            (
                rhoFluidA().value()*dE*dE*dE/
                2.0/M_PI/(cleanInterfaceSurfTension().value() + SMALL)
            )
        );
    }
    else
    {
        scalarField sigmaE = 
            linearEdgeInterpolate(surfaceTension())().internalField()
          + SMALL;

        const scalarField& dE =aMesh().lPN();
        
        CoNum = gMax
        (
            DB().deltaT().value()/
            sqrt
            (
                rhoFluidA().value()*dE*dE*dE/
                2.0/M_PI/sigmaE
            )
        );
    }

    return CoNum;
}


///RG-TQ updatePressureDensity
void freeSurface::updatePressureDensity
(
//  dimensionedScalar rhoB,
//  dimensionedScalar rhoV,
  dimensionedScalar pCorr  
)
{
//    rhoFluidB_ = rhoB;
//    rhoFluidV_ = rhoV;
    pCorr_ = pCorr;
}

///RG-TQ updatePressureDensity
dimensionedScalar freeSurface::updatePressureDensity2
(
//  dimensionedScalar rhoB,
//  dimensionedScalar rhoV,
)
{
//    rhoFluidB_ = rhoB;
//    rhoFluidV_ = rhoV;
    dimensionedScalar pCorr = pCorr_;
return pCorr;
}

///RG-TQ update Volume
void freeSurface::updateVolume
(
  dimensionedScalar Vv,
  dimensionedScalar VvInitial
)
{
    Vv_ = Vv;

    VvInitial_ = VvInitial;
}

void freeSurface::initializeConcentration()

{
        // Compute temperature at the interface
/*
        scalarField DnA = mesh().boundary()[aPatchID()].deltaCoeffs();

        scalarField DnB = interpolatorBA().faceInterpolate
        (
            mesh().boundary()[bPatchID()].deltaCoeffs()
        );

        scalarField TPA = 
            T().boundaryField()[aPatchID()].patchInternalField();

        scalarField TPB = interpolatorBA().faceInterpolate
        (
            T().boundaryField()[bPatchID()].patchInternalField()
        );

        scalarField beta =
                 kFluidA().value()*DnA + kFluidB().value()*DnB + VSMALL;

        scalarField T0 = (kFluidA().value()*DnA*TPA 
                        + kFluidB().value()*DnB*TPB)/beta;

        dimensionedScalar TSat = gAverage(T0);
*/
//    	Info<< "TSat1 " << TSat1.value() << endl;
        dimensionedScalar TSat = (TLeft().value() + TRight().value())/2.0;
//	calculate euqlibrium pressure for pure substance through Antoine euqation
	dimensionedScalar p1 = 133.322368*pow(10.0,antoineA().value()-antoineB().value()/(antoineC().value()+TSat-273.15));

//	dimensionedScalar pAtm = 101325.0;
//	dimensionedScalar pAir = pAtm - pSat;

	dimensionedScalar pAir = pTotal() - p1;
	Info<<endl<<"==================freeSurface::initializeConcentration()========================"<<endl;
    	Info<< "p1 " << p1.value() << " pAir " << pAir.value()<< endl;

	DfFluidB_.value() = DfFluidB0().value()*(101325.0/(pTotal().value()));
	DfFluidA_.value() = DfFluidB().value();

//	calculate density through equation of state for ideal gas
	dimensionedScalar rhoVapor1("rhoVapor1", dimDensity, p1.value()/(gasConstant1().value()*TSat.value()));
	dimensionedScalar rhoVaporD("rhoVaporD", dimDensity, pAir.value()/(gasConstantD().value()*TSat.value()));
//	dimensionedScalar rhoVapor1 = p1/(gasConstant1()*TSat);
//	dimensionedScalar rhoVaporD = pAir/(gasConstantD()*TSat);

	rhoB1 += rhoVapor1.value();
	ncB1 += rhoVapor1.value()/molarMass1().value();
//	Info << min(rhoB1) << " < rhoB1 < " << max(rhoB1) << endl;

	rhoFluidB1_.value() = rhoVapor1.value();
	rhoFluidB_.value() = (rhoVapor1 + rhoVaporD).value(); 
//	rhoFluidB_.value() = rhoVaporD.value();
//	rhoVaporRef_.value() = (rhoVapor1 + rhoVaporD).value(); 

	dimensionedScalar conini1 = p1/pTotal();

	con()=fluidIndicator()*(1.0 - conini1.value())+ conini1.value();
	rho1()=(1.0-fluidIndicator())*rhoVapor1 + fluidIndicator()*rhoFluidA();
	volScalarField nc1=(1.0-fluidIndicator())*(rhoVapor1/molarMass1()) + fluidIndicator()*(rhoFluidA()/molarMass1());
	nt() =(1.0-fluidIndicator())*(rhoVapor1/molarMass1()+rhoVaporD/molarMassD()) + fluidIndicator()*(rhoFluidA()/molarMass1());
	MB1_.value() = (fvc::domainIntegrate((1.0-fluidIndicator())*rhoVapor1)).value();
	Mair0_.value() = (fvc::domainIntegrate((1.0-fluidIndicator())*rhoVaporD)).value();
	c0_.value() = conini1.value();
//	MBD_.value() = (fvc::domainIntegrate((1.0-fluidIndicator())*rhoVaporD)).value();
//	MB1_.value() = (fvc::domainIntegrate((1.0-fluidIndicator())*rho*con)).value();
//	MBD_.value() = (fvc::domainIntegrate((1.0-fluidIndicator())*rho*(1.0 - con))).value();

	Info << "Initial vapor mass [MB1]: " << MB1().value() << endl;
	Info << "Initial air mass [Mair0]: " << Mair0().value() << endl;
    	Info << "initial rhoVapor1 [rhoFluidB1]: " << rhoVapor1.value() << endl;
    	Info << "initial rhoVapor+Air [rhoFluidB]: " << rhoFluidB_.value()<< endl;
//    	Info << "(initial) rhoVaporRef----> " << rhoVaporRef_.value()<< endl;
	Info << "initial concentration in vapor [c0]: " << conini1.value() << endl;
	Info << "initial molar density in vapor [n1] ----> " << min(nc1) <<max(nc1) <<endl;
	Info << "initial molar density in vapor [nt] ----> " << min(nt()) <<max(nt()) <<endl;
	Info << "DfFluidB_.value() ----> " << DfFluidB_.value() << endl<<endl;

	muFluidB_ = muFluidB_*(1.0-conini1.value())+muFluidB1_*conini1.value();
	CpFluidB_ = CpFluidB_*(1.0-conini1.value())+CpFluidB1_*conini1.value();
	kFluidB_ = kFluidB_*(1.0-conini1.value())+kFluidB1_*conini1.value();
}


void freeSurface::initializeMass()

{
        // Compute temperature at the interface
/*
        scalarField DnA = mesh().boundary()[aPatchID()].deltaCoeffs();

        scalarField DnB = interpolatorBA().faceInterpolate
        (
            mesh().boundary()[bPatchID()].deltaCoeffs()
        );

        scalarField TPA = 
            T().boundaryField()[aPatchID()].patchInternalField();

        scalarField TPB = interpolatorBA().faceInterpolate
        (
            T().boundaryField()[bPatchID()].patchInternalField()
        );

        scalarField beta =
                 kFluidA().value()*DnA + kFluidB().value()*DnB + VSMALL;

        scalarField T0 = (kFluidA().value()*DnA*TPA 
                        + kFluidB().value()*DnB*TPB)/beta;

        dimensionedScalar TSat = gAverage(T0);

	TSat.value() = (TLeft().value() + TRight().value())/2.0;
*/
	dimensionedScalar TSat = (TLeft() + TRight())/2.0;

//	calculate euqlibrium pressure for pure substance through Antoine euqation
	dimensionedScalar p1 = 133.322368*pow(10.0,antoineA().value()-antoineB().value()/(antoineC().value()+TSat.value()-273.15));

//	dimensionedScalar pAtm = 101325.0;
//	dimensionedScalar pAir = pAtm - pSat;

	dimensionedScalar pAir = pTotal() - p1;

	DfFluidB_.value() = DfFluidB0().value()*(101325.0/(pTotal().value()));
	DfFluidA_.value() = DfFluidB().value();

//	calculate density through equation of state for ideal gas
	dimensionedScalar rhoVapor1 = p1/(gasConstant1()*TSat);
	dimensionedScalar rhoVaporD = pAir/(gasConstantD()*TSat);
//	rhoVaporRef_.value() = (rhoVapor1 + rhoVaporD).value(); 

	MB1_.value() = (fvc::domainIntegrate((1.0-fluidIndicator())*rhoVapor1)).value();

	MBD_.value() = (fvc::domainIntegrate((1.0-fluidIndicator())*rhoVaporD)).value();

//	MB1_.value() = (fvc::domainIntegrate((1.0-fluidIndicator())*rho*con)).value();
//	MBD_.value() = (fvc::domainIntegrate((1.0-fluidIndicator())*rho*(1.0 - con))).value();
	Info<<"==================freeSurface::initializeMass()========================"<<endl;
	Info<< "freeSurface::initializeMass::Initial vapor mass: " << MB1().value() << endl;
	Info<< "freeSurface::initializeMass::Initial air mass: " << MBD().value() << endl;
	Info<< "freeSurface::initializeMass::DfFluidB_.value() ----> " << DfFluidB_.value() << endl<<endl;
//    	Info << "(initial) rhoVaporRef----> " << rhoVaporRef_.value()<< endl;

}

void freeSurface::updateProperties()
{
    muFluidA_ = dimensionedScalar(this->lookup("muFluidA"));

    muFluidB_ = dimensionedScalar(this->lookup("muFluidB"));
    muFluidB1_ = dimensionedScalar(this->lookup("muFluidB1"));

    rhoFluidA_ = dimensionedScalar(this->lookup("rhoFluidA"));

    //rhoFluidB_ = dimensionedScalar(this->lookup("rhoFluidB"));

//    rhoFluidV_ = dimensionedScalar(this->lookup("rhoFluidV"));

    ///RG-TQ
    kFluidA_ = dimensionedScalar(this->lookup("kFluidA"));

    kFluidB_ = dimensionedScalar(this->lookup("kFluidB"));
    kFluidB1_ = dimensionedScalar(this->lookup("kFluidB1"));

    CpFluidA_ = dimensionedScalar(this->lookup("CpFluidA"));

    CpFluidB_ = dimensionedScalar(this->lookup("CpFluidB")); 
    CpFluidB1_ = dimensionedScalar(this->lookup("CpFluidB1")); 

    betaFluidA_ = dimensionedScalar(this->lookup("betaFluidA"));

    //betaFluidB_ = dimensionedScalar(this->lookup("betaFluidB"));  

    DfFluidA_ = dimensionedScalar(this->lookup("DfFluidA"));

    DfFluidB_ = dimensionedScalar(this->lookup("DfFluidB"));
    DfFluidB0_ = dimensionedScalar(this->lookup("DfFluidB0"));

    DfFluidB_.value() = DfFluidB0().value()*(101325.0/(pTotal().value()));
    DfFluidA_.value() = DfFluidB().value();

    latentHeat_ = dimensionedScalar(this->lookup("latentHeat"));   

//    gasConstantWater_ = dimensionedScalar(this->lookup("gasConstantWater")); 

//    gasConstantAir_ = dimensionedScalar(this->lookup("gasConstantAir")); 

    TRef_ = dimensionedScalar(this->lookup("referenceTemperature"));  

    g_ = dimensionedVector(this->lookup("g"));

    cleanInterfaceSurfTension_ = 
        dimensionedScalar(this->lookup("surfaceTension"));

    ///RG-TQ
    tempCoeffSurfTension_ = 
        dimensionedScalar(this->lookup("tempCoeffSurfTension"));

	muFluidB_ = muFluidB_*(1.0-c0_.value())+muFluidB1_*c0_.value();
	CpFluidB_ = CpFluidB_*(1.0-c0_.value())+CpFluidB1_*c0_.value();
	kFluidB_ = kFluidB_*(1.0-c0_.value())+kFluidB1_*c0_.value();
	Mair0_ = dimensionedScalar(this->lookup("Mair0"));
	MB1_.value() = (fvc::domainIntegrate((1.0-fluidIndicator())*rho1())).value();
	Info<<"==================freeSurface::updateProperties()========================"<<endl;
	Info<< "Initial vapor mass [MB1] after updateProperties(): " << MB1().value() << endl;
	Info<< "Initial air mass [Mair0] after updateProperties(): " << Mair0().value() << endl;	
	Info<< "DfFluidB_.value() after updateProperties(): " << DfFluidB_.value() << endl <<endl<<endl;
}


void freeSurface::writeVTK() const
{
    aMesh().patch().writeVTK
    (
        DB().timePath()/"freeSurface",
        aMesh().patch(),
        aMesh().patch().points()
    );
}


void freeSurface::writeVTKControlPoints()
{
    // Write patch and points into VTK
    fileName name(DB().timePath()/"freeSurfaceControlPoints");
    OFstream mps(name + ".vtk");

    mps << "# vtk DataFile Version 2.0" << nl
        << name << ".vtk" << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl
        << "POINTS " << controlPoints().size() << " float" << nl;

    forAll(controlPoints(), pointI)
    {
        mps << controlPoints()[pointI].x() << ' '
            << controlPoints()[pointI].y() << ' '
            << controlPoints()[pointI].z() << nl;

    }

    // Write vertices
    mps << "VERTICES " << controlPoints().size() << ' ' 
        << controlPoints().size()*2 << nl;

    forAll(controlPoints(), pointI)
    {
        mps << 1 << ' ' << pointI << nl;
    }
}

tmp<areaVectorField> freeSurface::grad(areaScalarField& Ts)
{
    tmp<areaVectorField> tGradTs
    (
        new areaVectorField
        (
            IOobject
            (
                "grad(" + Ts.name() + ")",
                DB().timeName(),
                mesh(),
                IOobject::NO_READ
            ),
            fac::grad(Ts)
        )
    );
    areaVectorField& gradTs = tGradTs();

    bool calculatedBC = false;

    // Zero gradient extrapolation for calculated boundary
    forAll(Ts.boundaryField(), patchI)
    {
        if 
        (
            Ts.boundaryField()[patchI].type()
         == calculatedFaPatchScalarField::typeName
        )
        {
            calculatedBC = true;

            Ts.boundaryField()[patchI] =
                Ts.boundaryField()[patchI].patchInternalField();
        }
    }


    if (calculatedBC)
    {
        gradTs = fac::grad(Ts);

        label iCorr = 1;
        do
        {
            iCorr++;
        
            // Extraplate value at calculated boundary
            forAll(Ts.boundaryField(), patchI)
            {
                if 
                (
                    Ts.boundaryField()[patchI].type()
                 == calculatedFaPatchScalarField::typeName
                )
                {
                    vectorField d = aMesh().boundary()[patchI].delta();
                
                    scalarField dTs =
                        (
                            d
                          & gradTs.boundaryField()[patchI]
                           .patchInternalField()
                        );

                    Ts.boundaryField()[patchI] =
                        Ts.boundaryField()[patchI].patchInternalField()
                      + dTs;
                }
            }

            gradTs = fac::grad(Ts);
        }
        while(iCorr < 5);
    }   

    return tGradTs;
}

//z0928 freeSurface::correctCurvature()
void freeSurface::correctCurvature()
{
    // Correct curvature next to fixed patches

    areaScalarField& K =
        const_cast<areaScalarField&>
        (
            aMesh().faceCurvatures()
        );

    scalarField& KI = K.internalField();

    if (curvExtrapOrder_ == 0)
    {
        forAll(fixedFreeSurfacePatches_, patchI)
        {
            label fixedPatchID = 
                aMesh().boundary().findPatchID
                (
                    fixedFreeSurfacePatches_[patchI]
                );

            if(fixedPatchID == -1)
            {
                FatalErrorIn("freeSurface::freeSurface(...)")
                    << "Wrong faPatch name in the fixedFreeSurfacePatches list"
                        << " defined in the freeSurfaceProperties dictionary"
                        << abort(FatalError);
            }

            const labelList& eFaces =
                aMesh().boundary()[fixedPatchID].edgeFaces();

            const labelListList& fFaces = aMesh().patch().faceFaces();

            forAll(eFaces, edgeI)
            {
                const label& curFace = eFaces[edgeI];
                const labelList& curFaceFaces = fFaces[curFace];

                scalar avrK = 0.0;
                label counter = 0;

                forAll(curFaceFaces, faceI)
                {
                    label index = findIndex(eFaces, curFaceFaces[faceI]);
                    
                    if (index == -1)
                    {
                        avrK += K[curFaceFaces[faceI]];
                        counter++;
                    }
                }

                avrK /= counter;
                
                KI[curFace] = avrK;
            }
        }
    }
    else if (curvExtrapOrder_ == 1)
    {
        label counter = 0;
        do
        {
            counter++;
            
            K.correctBoundaryConditions();
            areaVectorField gradK = fac::grad(K);
            vectorField& gradKI = gradK.internalField();

            forAll(fixedFreeSurfacePatches_, patchI)
            {
                label fixedPatchID = 
                    aMesh().boundary().findPatchID
                    (
                        fixedFreeSurfacePatches_[patchI]
                    );

                if(fixedPatchID == -1)
                {
                    FatalErrorIn("freeSurface::freeSurface(...)")
                        << "Wrong faPatch name in the fixedFreeSurfacePatches "
                            << "list defined in the freeSurfaceProperties "
                            << "dictionary"
                            << abort(FatalError);
                }
            
                const labelList& eFaces =
                    aMesh().boundary()[fixedPatchID].edgeFaces();

                const labelListList& fFaces = aMesh().patch().faceFaces();

                const vectorField& fCentres = aMesh().areaCentres();

                forAll(eFaces, edgeI)
                {
                    const label& curFace = eFaces[edgeI];
                    const labelList& curFaceFaces = fFaces[curFace];

                    scalar avrK = 0.0;
                    label counter = 0;

                    forAll(curFaceFaces, faceI)
                    {
                        label index = findIndex(eFaces, curFaceFaces[faceI]);

                        if (index == -1)
                        {
                            vector dr = 
                                fCentres[curFace] 
                              - fCentres[curFaceFaces[faceI]];

                            avrK += KI[curFaceFaces[faceI]]
                                + (dr&gradKI[curFaceFaces[faceI]]);
                            counter++;
                        }
                    }
            
                    avrK /= counter;
                    
                    KI[curFace] = avrK;
                }
            }
        }
        while(counter<10);
    }
    else
    {
        FatalErrorIn("freeSurface::freeSurface(...)")
            << "Curvature extrapolation order is not set correctly"
                << abort(FatalError);
    }

    K.correctBoundaryConditions();
}

//z0312 freeSurface::correctPointNormals()
void freeSurface::correctPointNormals()
{
    // Correct normals for fixed patches points

    vectorField& N =
        const_cast<vectorField&>
        (
            aMesh().pointAreaNormals()
        ); 

    const labelListList& pFaces =
        aMesh().patch().pointFaces();

    const labelListList& fFaces =

        aMesh().patch().faceFaces();

    const faceList& faces = 
        aMesh().patch().localFaces();

    const pointField& points = 
        aMesh().patch().localPoints();


    // Wedge points
    forAll(aMesh().boundary(), patchI)
    {
        if (aMesh().boundary()[patchI].type() == wedgeFaPatch::typeName)
        {
            const wedgeFaPatch& wedgePatch =
                refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

            labelList patchPoints = wedgePatch.pointLabels();

            forAll(patchPoints, pointI)
            {
                label curPoint = patchPoints[pointI];

                labelHashSet faceSet;
                forAll(pFaces[curPoint], faceI)
                {
                    faceSet.insert(pFaces[curPoint][faceI]);

                }
                labelList curFaces = faceSet.toc();
                
                labelHashSet pointSet;

                pointSet.insert(curPoint);
                for(label i=0; i<curFaces.size(); i++)
                {
                    const labelList& facePoints = faces[curFaces[i]];
                    for(label j=0; j<facePoints.size(); j++)
                    {
                        if(!pointSet.found(facePoints[j]))
                        {
                            pointSet.insert(facePoints[j]);
                        }
                    }
                }
                pointSet.erase(curPoint);
                labelList curPoints = pointSet.toc();

                
                labelHashSet addPointsSet;
                forAll(curPoints, pointI)
                {
                    label index = 
                        findIndex(patchPoints, curPoints[pointI]);

                    if (index != -1)
                    {
                        addPointsSet.insert(curPoints[pointI]);
                    }  
                }
                addPointsSet.insert(curPoint);
                labelList curAddPoints = addPointsSet.toc();


                if (curPoints.size() + curAddPoints.size() >= 5)
                {
                    vectorField allPoints
                    (
                        curPoints.size()+curAddPoints.size()
                    );
                    scalarField W(curPoints.size()+curAddPoints.size(), 1.0);
                    label I = -1;
                    for(label i=0; i<curPoints.size(); i++)
                    {
                        I++;
                        allPoints[I] = points[curPoints[i]];
                        W[I] = 1.0/magSqr(allPoints[I] - points[curPoint]);
                    }
                    for(label i=0; i<curAddPoints.size(); i++)
                    {
                        I++;
                        allPoints[I] = 
                            transform
                            (
                                wedgePatch.faceT(),
                                points[curAddPoints[i]]
                            );
                        W[I] = 1.0/magSqr(allPoints[I] - points[curPoint]);
                    }

                    // Transforme points
                    vector origin = points[curPoint];
                    vector axis = N[curPoint]/mag(N[curPoint]);
                    vector dir = (allPoints[0] - points[curPoint]);
                    dir -= axis*(axis&dir);
                    dir /= mag(dir);
                    coordinateSystem cs("cs", origin, axis, dir);
                    
                    forAll(allPoints, pI)
                    {
                        allPoints[pI] = cs.localPosition(allPoints[pI]);
                    }
                    
                    scalarRectangularMatrix M
                    (
                        allPoints.size(),
                        5,
                        0.0
                    );

                    for(label i = 0; i < allPoints.size(); i++)
                    {
                        M[i][0] = sqr(allPoints[i].x());
                        M[i][1] = sqr(allPoints[i].y());
                        M[i][2] = allPoints[i].x()*allPoints[i].y();
                        M[i][3] = allPoints[i].x();
                        M[i][4] = allPoints[i].y();
                    }
                    
                    scalarSquareMatrix MtM(5, 0.0);

                    for (label i = 0; i < MtM.n(); i++)
                    {

                        for (label j = 0; j < MtM.m(); j++)
                        {
                            for (label k = 0; k < M.n(); k++)
                            {
                                MtM[i][j] += M[k][i]*M[k][j]*W[k];
                            }
                        }
                    }
                    
                    scalarField MtR(5, 0);

                    for (label i=0; i<MtR.size(); i++)
                    {
                        for (label j=0; j<M.n(); j++)
                        {
                            MtR[i] += M[j][i]*allPoints[j].z()*W[j];
                        }
                    }
            
                    scalarSquareMatrix::LUsolve(MtM, MtR);

                    vector curNormal = vector(MtR[3], MtR[4], -1);

                    
                    curNormal = cs.globalVector(curNormal);
                    
                    curNormal *= sign(curNormal&N[curPoint]);
                    
                    N[curPoint] = curNormal;
                }
            }
        }
    }
    
    // Fixed boundary points

    forAll(fixedFreeSurfacePatches_, patchI)
    {
        label fixedPatchID = 
            aMesh().boundary().findPatchID
            (
                fixedFreeSurfacePatches_[patchI]
            );

        if(fixedPatchID == -1)
        {
            FatalErrorIn("freeSurface::freeSurface(...)")
                << "Wrong faPatch name in the fixedFreeSurfacePatches list"
                    << " defined in the freeSurfaceProperties dictionary"
                    << abort(FatalError);
        }

        const labelList& pLabels = 
            aMesh().boundary()[fixedPatchID].pointLabels();

        const labelList& eFaces =
            aMesh().boundary()[fixedPatchID].edgeFaces();
 
        forAll(pLabels, pointI)
        {
            label curPoint = pLabels[pointI];

            labelHashSet faceSet;
            forAll(pFaces[curPoint], faceI)
            {
                faceSet.insert(pFaces[curPoint][faceI]);
            }

            labelList curFaces = faceSet.toc();

            forAll(curFaces, faceI)
            {
                const labelList& curFaceFaces =
                    fFaces[curFaces[faceI]];
                
                forAll(curFaceFaces, fI)
                {
                    label curFaceFace = curFaceFaces[fI];
                        
                    label index = findIndex(eFaces, curFaceFace);

                    if( (index==-1) && !faceSet.found(curFaceFace) )
                    {
                        faceSet.insert(curFaceFace);
                    }
                }
            }
            curFaces = faceSet.toc();

            labelHashSet pointSet;

            pointSet.insert(curPoint);
            for(label i=0; i<curFaces.size(); i++)
            {
                const labelList& fPoints = faces[curFaces[i]];
                for(label j=0; j<fPoints.size(); j++)
                {
                    if(!pointSet.found(fPoints[j]))
                    {
                        pointSet.insert(fPoints[j]);
                    }
                }
            }

            pointSet.erase(curPoint);

            labelList curPoints = pointSet.toc();

            // LS quadric fit
            vectorField allPoints(curPoints.size());
            scalarField W(curPoints.size(), 1.0);
            for(label i=0; i<curPoints.size(); i++)
            {
                allPoints[i] = points[curPoints[i]];
                W[i] = 1.0/magSqr(allPoints[i] - points[curPoint]);
            }

            // Transforme points
            vector origin = points[curPoint];
            vector axis = N[curPoint]/mag(N[curPoint]);
            vector dir = (allPoints[0] - points[curPoint]);
            dir -= axis*(axis&dir);
            dir /= mag(dir);
            coordinateSystem cs("cs", origin, axis, dir);

            forAll(allPoints, pI)
            {
                allPoints[pI] = cs.localPosition(allPoints[pI]);
            }

            scalarRectangularMatrix M
            (
                allPoints.size(),
                5,
                0.0
            );

            for(label i = 0; i < allPoints.size(); i++)
            {
                M[i][0] = sqr(allPoints[i].x());
                M[i][1] = sqr(allPoints[i].y());
                M[i][2] = allPoints[i].x()*allPoints[i].y();
                M[i][3] = allPoints[i].x();
                M[i][4] = allPoints[i].y();
            }

            scalarSquareMatrix MtM(5, 0.0);

            for (label i = 0; i < MtM.n(); i++)
            {
                for (label j = 0; j < MtM.m(); j++)
                {
                    for (label k = 0; k < M.n(); k++)
                    {
                        MtM[i][j] += M[k][i]*M[k][j]*W[k];
                    }
                }
            }

            scalarField MtR(5, 0);

            for (label i=0; i<MtR.size(); i++)
            {
                for (label j=0; j<M.n(); j++)
                {
                    MtR[i] += M[j][i]*allPoints[j].z()*W[j];
                }
            }

            scalarSquareMatrix::LUsolve(MtM, MtR);

            vector curNormal = vector(MtR[3], MtR[4], -1);

            curNormal = cs.globalVector(curNormal);

            curNormal *= sign(curNormal&N[curPoint]);

            N[curPoint] = curNormal;
        }
    }

    // Correcte wedge points
    forAll (aMesh().boundary(), patchI)
    {
        if (aMesh().boundary()[patchI].type() == wedgeFaPatch::typeName)
        {
            const wedgeFaPatch& wedgePatch =
                refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

            labelList patchPoints = wedgePatch.pointLabels();

            vector n =
                transform
                (
                    wedgePatch.edgeT(),
                    wedgePatch.centreNormal()
                );

            n /= mag(n);

            forAll (patchPoints, pointI)
            {
                N[patchPoints[pointI]]
                    -= n*(n&N[patchPoints[pointI]]);
            }
        }
    }


    // Boundary points correction
    forAll (aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().correctPatchPointNormals(patchI) 
        && !aMesh().boundary()[patchI].coupled()
        )
        {
            if (aMesh().boundary()[patchI].ngbPolyPatchIndex() == -1)
            {
                FatalErrorIn
                    (
                        "void freeSurface::correctPointNormals const"
                    )   << "Neighbour polyPatch index is not defined "
                        << "for faPatch " << aMesh().boundary()[patchI].name()
                        << abort(FatalError);
            }

            labelList patchPoints = aMesh().boundary()[patchI].pointLabels();

            vectorField n = 
                aMesh().boundary()[patchI].ngbPolyPatchPointNormals();

            forAll (patchPoints, pointI)
            {
                N[patchPoints[pointI]]
                    -= n[pointI]*(n[pointI]&N[patchPoints[pointI]]);
            }
        }
    }


    N /= mag(N);
}

//z0312 freeSurface::correctPointDisplacement
void freeSurface::correctPointDisplacement
(
    const scalarField& sweptVolCorr, 
    vectorField& displacement
)
{
    const labelListList& pFaces =
        aMesh().patch().pointFaces();

    const faceList& faces = 
        aMesh().patch().localFaces();

    const pointField& points = 
        aMesh().patch().localPoints();

    forAll(fixedFreeSurfacePatches_, patchI)
    {
        label fixedPatchID =
            aMesh().boundary().findPatchID
            (

                fixedFreeSurfacePatches_[patchI]
            );

        if(fixedPatchID == -1)
        {
            FatalErrorIn("freeSurface::freeSurface(...)")
                 << "Wrong faPatch name in the fixedFreeSurfacePatches list"
                    << " defined in the freeSurfaceProperties dictionary"
                    << abort(FatalError);
        }

        const labelList& pLabels = 
            aMesh().boundary()[fixedPatchID].pointLabels();

        const labelList& eFaces =
            aMesh().boundary()[fixedPatchID].edgeFaces();

        labelHashSet pointSet;

        forAll(eFaces, edgeI)
        {
            label curFace = eFaces[edgeI];

            const labelList& curPoints = faces[curFace];
            
            forAll(curPoints, pointI)
            {
                label curPoint = curPoints[pointI];
                label index = findIndex(pLabels, curPoint);

                if (index == -1)
                {
                    if (!pointSet.found(curPoint))
                    {
                        pointSet.insert(curPoint);
                    }
                }
            }
        }


        labelList corrPoints = pointSet.toc();

        labelListList corrPointFaces(corrPoints.size());

        forAll(corrPoints, pointI)
        {
            label curPoint = corrPoints[pointI];

            labelHashSet faceSet;
            
            forAll(pFaces[curPoint], faceI)
            {
                label curFace = pFaces[curPoint][faceI];

                label index = findIndex(eFaces, curFace);

                if (index != -1)
                {
                    if (!faceSet.found(curFace))
                    {
                        faceSet.insert(curFace);
                    }
                }
            }

            corrPointFaces[pointI] = faceSet.toc();
        }

        forAll(corrPoints, pointI)
        {
            label curPoint = corrPoints[pointI];

            scalar curDisp = 0;

            const labelList& curPointFaces = corrPointFaces[pointI];

            forAll(curPointFaces, faceI)
            {
                const face& curFace = faces[curPointFaces[faceI]];

                label ptInFace = curFace.which(curPoint);
                label next = curFace.nextLabel(ptInFace);
                label prev = curFace.prevLabel(ptInFace);
                
                vector a = points[next] - points[curPoint];
                vector b = points[prev] - points[curPoint];
                const vector& c = pointsDisplacementDir()[curPoint];

                curDisp += 2*sweptVolCorr[curPointFaces[faceI]]/((a^b)&c);
            }

            curDisp /= curPointFaces.size();

            displacement[curPoint] = 
                curDisp*pointsDisplacementDir()[curPoint];
        }
    }
}

//z092812 freeSurface::smoothCurvature()
void freeSurface::smoothCurvature()
{
    areaScalarField& oldK =
        const_cast<areaScalarField&>
        (
            aMesh().faceCurvatures()
        ); 

    areaScalarField K
    (
        IOobject
        (
            "K",
            DB().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aMesh(),
        dimless/dimLength,
        fixedGradientFaPatchScalarField::typeName
    );

    K.internalField() = oldK.internalField();

    forAll(K.boundaryField(), patchI)
    {
        if
        (
            K.boundaryField()[patchI].type()
         == fixedGradientFaPatchScalarField::typeName
        )
        {
            fixedGradientFaPatchScalarField& Kp =
                refCast<fixedGradientFaPatchScalarField>
                (
                    K.boundaryField()[patchI]
                );
            
            Kp.gradient() = 0;
        }
    }

    K.correctBoundaryConditions();


    label counter = 0;

    // Set gradient
    do
    {
        counter++;

        areaVectorField gradK = fac::grad(K);

        forAll(K.boundaryField(), patchI)
        {
            if
            (
                K.boundaryField()[patchI].type()
             == fixedGradientFaPatchScalarField::typeName
            )
            {
                fixedGradientFaPatchScalarField& Kp =
                    refCast<fixedGradientFaPatchScalarField>
                    (
                        K.boundaryField()[patchI]
                    );

                Kp.gradient() =
                (
                    aMesh().boundary()[patchI].edgeNormals()
                   &gradK.boundaryField()[patchI].patchInternalField()
                );
            }
        }

        K.correctBoundaryConditions();
    }
    while(counter<5);


    areaScalarField indicator =
        fac::div
        (
            fac::lnGrad(K)*aMesh().magLe()
          - (aMesh().Le()&fac::interpolate(fac::grad(K)))
        );

    scalar minIndicator = GREAT;
    label refFace = -1;

    forAll(indicator, faceI)
    {
        if (mag(indicator[faceI]) < minIndicator)

        {
            minIndicator = mag(indicator[faceI]);
            refFace = faceI;
        }
    }

    scalar refK = K[refFace];

    counter = 0;

    do
    {
        counter++;

        faScalarMatrix KEqn
        (
            fam::laplacian(K) 
         == fac::div(aMesh().Le()&fac::interpolate(fac::grad(K)))
        );

        KEqn.setReference(refFace, refK);
        KEqn.solve();
    }
    while(counter<2);

    oldK = K;
}

void freeSurface::updateNGradUn()
{
    if (fvcNGradUn_)
    {
        Info << "Update normal derivative of normal velocity using fvc" 
            << endl;

        volVectorField phiU = fvc::reconstruct(phi_);

        vectorField nA = mesh().boundary()[aPatchID()].nf();
    
        scalarField UnP = 
            (nA&phiU.boundaryField()[aPatchID()].patchInternalField());

        scalarField UnFs = 
            phi_.boundaryField()[aPatchID()]
           /mesh().magSf().boundaryField()[aPatchID()];

        nGradUn() = 
            (UnFs - UnP)*mesh().deltaCoeffs().boundaryField()[aPatchID()];

//	0928update
        bool secondOrderCorrection = true;

        if (secondOrderCorrection)
        {
            // Correct normal component of phiU 
            // befor gradient calculation
            forAll(phiU.boundaryField(), patchI)
            {
                vectorField n = 
                    mesh().Sf().boundaryField()[patchI]
                   /mesh().magSf().boundaryField()[patchI];

                phiU.boundaryField()[patchI] +=
                    n
                   *(
                       (
                           phi_.boundaryField()[patchI]
                          /mesh().magSf().boundaryField()[patchI]
                       )
                     - (n&phiU.boundaryField()[patchI])
                    );
            }

            // Calc gradient
            tensorField gradPhiUp = 
                fvc::grad(phiU)().boundaryField()[aPatchID()]
               .patchInternalField();

            nGradUn() = 2*nGradUn() - (nA&(gradPhiUp&nA));
        }
    }
    else
    {
        Info << "Update normal derivative of normal velocity using fac"
            << endl;

        nGradUn() = -fac::div(Us())().internalField();
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

