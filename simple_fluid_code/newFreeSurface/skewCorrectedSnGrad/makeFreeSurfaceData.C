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

\*---------------------------------------------------------------------------*/

#include "freeSurface.H"
#include "primitivePatchInterpolation.H"
#include "wedgeFaPatch.H"
#include "wedgeFaPatchFields.H"
#include "wallFvPatch.H"
#include "slipFaPatchFields.H"

#include "fixedValueFaPatchFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void freeSurface::makeInterpolators()
{
    if (debug)
    {
        Info<< "freeSurface::makeInterpolators() : "
            << "making pathc to patch interpolator"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if 
    (
        interpolatorBAPtr_ ||  
        interpolatorABPtr_
    )
    {
        FatalErrorIn("freeSurface::makeInterpolators()")
            << "patch to patch interpolators already exists"
            << abort(FatalError);
    }


    if(aPatchID() == -1)
    {
        FatalErrorIn("freeSurface::makeInterpolators()")
            << "Free surface patch A not defined."
            << abort(FatalError);
    }


    if(bPatchID() == -1)
    {
        FatalErrorIn("freeSurface::makeInterpolators()")
            << "Free surface patch B not defined."
            << abort(FatalError);
    }

//     patchToPatchInterpolation::setDirectHitTol(1e-2);

    interpolatorBAPtr_ = new IOpatchToPatchInterpolation
    (
        IOobject
        (
            "baInterpolator",
            DB().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh().boundaryMesh()[bPatchID()],
        mesh().boundaryMesh()[aPatchID()],
        intersection::VISIBLE
        // intersection::HALF_RAY
    );

    
    const scalarField& faceDistBA = 
        interpolatorBAPtr_->faceDistanceToIntersection();

    forAll(faceDistBA, faceI)
    {
        if(mag(faceDistBA[faceI] - GREAT) < SMALL)
        {
            FatalErrorIn("freeSurface::makeInterpolators()")
                << "Error in B-to-A face patchToPatchInterpolation."
                << abort(FatalError);            
        }
    }

    const scalarField& pointDistBA = 
        interpolatorBAPtr_->pointDistanceToIntersection();

    forAll(pointDistBA, pointI)
    {
        if(mag(pointDistBA[pointI] - GREAT) < SMALL)
        {
            FatalErrorIn("freeSurface::makeInterpolators()")
                << "Error in B-to-A point patchToPatchInterpolation."
                << abort(FatalError);            
        }
    }


    interpolatorABPtr_ = new IOpatchToPatchInterpolation
    (
        IOobject
        (
            "abInterpolator",
            DB().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh().boundaryMesh()[aPatchID()],
        mesh().boundaryMesh()[bPatchID()],
        intersection::VISIBLE
        // intersection::HALF_RAY
    );


    const scalarField& faceDistAB = 
        interpolatorABPtr_->faceDistanceToIntersection();

    forAll(faceDistAB, faceI)
    {
        if(mag(faceDistAB[faceI] - GREAT) < SMALL)
        {
            FatalErrorIn("freeSurface::makeInterpolators()")
                << "Error in A-to-B face patchToPatchInterpolation."
                << abort(FatalError);            
        }
    }

    const scalarField& pointDistAB = 
        interpolatorABPtr_->pointDistanceToIntersection();

    forAll(pointDistAB, pointI)
    {
        if(mag(pointDistAB[pointI] - GREAT)<SMALL)
        {
            FatalErrorIn("freeSurface::makeInterpolators()")
                << "Error in A-to-B point patchToPatchInterpolation."
                << abort(FatalError);            
        }
    }


    Info << "\nCheck A-to-B and B-to-A interpolators" << endl;

    scalar maxDist = max
    (
        mag
        (
            interpolatorABPtr_->faceInterpolate
            (
                vectorField(mesh().boundaryMesh()[aPatchID()]
               .faceCentres())
            )
          - mesh().boundaryMesh()[bPatchID()].faceCentres()
        )
    );

    scalar maxDistPt = max
    (
        mag
        (
            interpolatorABPtr_->pointInterpolate
            (
                vectorField(mesh().boundaryMesh()[aPatchID()]
               .localPoints())
            )
          - mesh().boundaryMesh()[bPatchID()].localPoints()
        )
    );

    Info << "A-to-B interpolation error, face: " << maxDist
        << ", point: " << maxDistPt << endl;


    maxDist = max
    (
        mag
        (
            interpolatorBAPtr_->faceInterpolate
            (
                vectorField
                (
                    mesh().boundaryMesh()[bPatchID()].faceCentres()
                )
            )
          - mesh().boundaryMesh()[aPatchID()].faceCentres()
        )
    );

    maxDistPt = max
    (
        mag
        (
            interpolatorBAPtr_->pointInterpolate
            (
                vectorField
                (
                    mesh().boundaryMesh()[bPatchID()].localPoints()
                )
            )
          - mesh().boundaryMesh()[aPatchID()].localPoints()
        )
    );

    Info << "B-to-A interpolation error, face: " << maxDist
        << ", point: " << maxDistPt << endl;
}


void freeSurface::makeControlPoints()
{
    if (debug)
    {
        Info<< "freeSurface::makeControlPoints() : "
            << "making control points"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (controlPointsPtr_)
    {
        FatalErrorIn("freeSurface::makeInterpolators()")
            << "patch to patch interpolators already exists"
            << abort(FatalError);
    }

    IOobject controlPointsHeader
    (
        "controlPoints",
        DB().timeName(),
        mesh(),
        IOobject::MUST_READ
    );

    if (controlPointsHeader.headerOk())
    {
        controlPointsPtr_ =
            new vectorIOField
            (
                IOobject
                (
                    "controlPoints",
                    DB().timeName(),
                    mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            );
    }
    else
    {
        controlPointsPtr_ =
            new vectorIOField
            (
                IOobject
                (
                    "controlPoints",
                    DB().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                aMesh().areaCentres().internalField()
            );

//         initializeControlPointsPosition();
    }
}

void freeSurface::makeTi()
{
    if (debug)
    {
        Info<< "freeSurface::makeTi() : "
            << "making interfacial temperature"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (TiPtr_)
    {
        FatalErrorIn("freeSurface::makeTi()")
            << "Ti already exists"
            << abort(FatalError);
    }

    IOobject TiHeader
    (
        "Ti",
        DB().timeName(),
        mesh(),
        IOobject::MUST_READ
    );


    if (TiHeader.headerOk())
    {
/*
        TiPtr_ =
            new areaScalarField
            (
                IOobject
                (
                    "Ti",
                    DB().timeName(),
                    mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            );
*/
    }

    else
    {

        wordList patchFieldTypes
        (
            aMesh().boundary().size(),
//            zeroGradientFaPatchVectorField::typeName
           calculatedFaPatchScalarField::typeName
        );
/*
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
*/
        TiPtr_ =
            new areaScalarField
            (
                IOobject
                (
                    "Ti",
                    DB().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
            aMesh(),
            dimensioned<scalar>("Ti", dimTemperature, 0),
            patchFieldTypes
//                aMesh().areaCentres().internalField()
            );
    }
}

void freeSurface::makeTiDiff()
{
    if (debug)
    {
        Info<< "freeSurface::makeTiDiff() : "
            << "making interfacial temperature"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (TiDiffPtr_)
    {
        FatalErrorIn("freeSurface::makeTiDiff()")
            << "TiDiff already exists"
            << abort(FatalError);
    }
        wordList patchFieldTypes
        (
            aMesh().boundary().size(),
//            zeroGradientFaPatchVectorField::typeName
	calculatedFaPatchScalarField::typeName
        );
/*
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
*/
        TiDiffPtr_ =
            new areaScalarField
            (
                IOobject
                (
                    "TiDiff",
                    DB().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
            aMesh(),
            dimensioned<scalar>("TiDiff", dimTemperature, 0),
            patchFieldTypes
//                aMesh().areaCentres().internalField()
            );
    
}
void freeSurface::makedTi()
{
    if (debug)
    {
        Info<< "freeSurface::makedTi() : "
            << "making temperature jump across the interface (Tl-Tv)"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (dTiPtr_)
    {
        FatalErrorIn("freeSurface::makedTi()")
            << "dTi already exists"
            << abort(FatalError);
    }

    IOobject dTiHeader
    (
        "dTi",
        DB().timeName(),
        mesh(),
        IOobject::MUST_READ
    );


    if (dTiHeader.headerOk())
    {
/*
        dTiPtr_ =
            new areaScalarField
            (
                IOobject
                (
                    "dTi",
                    DB().timeName(),
                    mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            );
*/
    }

    else
    {
        wordList patchFieldTypes
        (
            aMesh().boundary().size(),
//             zeroGradientFaPatchVectorField::typeName
            calculatedFaPatchScalarField::typeName
        );
/*
        forAll(aMesh().boundary(), patchI)
        {
            label ngbPolyPatchID = 
                aMesh().boundary()[patchI].ngbPolyPatchIndex();

            if (ngbPolyPatchID != -1)
            {
                if
                (
                    dTi().boundaryField()[ngbPolyPatchID].type()
                 == fixedValueFvPatchScalarField::typeName
                )
                {
                    patchFieldTypes[patchI] =
                        fixedValueFaPatchScalarField::typeName;
                }
            }
        }
*/
        dTiPtr_ =
            new areaScalarField
            (
                IOobject
                (
                    "dTi",
                    DB().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
            aMesh(),
            dimensioned<scalar>("dTi", dimTemperature, 0),
            patchFieldTypes
//                aMesh().areaCentres().internalField()
            );
    }
}
void freeSurface::makegradTi()
{
    if (debug)
    {
        Info<< "freeSurface::makegradTi() : "
            << "making free-surface velocity field"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (gradTiPtr_)
    {
        FatalErrorIn("freeSurface::makegradTi()")
            << "free-surface velocity field already exists"
            << abort(FatalError);
    }

    wordList patchFieldTypes
    (
        aMesh().boundary().size(),
        zeroGradientFaPatchVectorField::typeName
//	fixedValueFaPatchVectorField::typeName
    );

    gradTiPtr_ = new areaVectorField
    (
        IOobject
        (
            "gradTi",
            DB().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh(),
        dimensioned<vector>("gradTi", dimless, vector::zero),
        patchFieldTypes
    );


}



void freeSurface::makemassFlux()
{
    if (debug)
    {
        Info<< "freeSurface::makemassFlux() : "
            << "making interfacial temperature"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (massFluxPtr_)
    {
        FatalErrorIn("freeSurface::makemassFlux()")
            << "massFlux already exists"
            << abort(FatalError);
    }

    IOobject massFluxHeader
    (
        "massFlux",
        DB().timeName(),
        mesh(),
        IOobject::MUST_READ
    );


    if (massFluxHeader.headerOk())
    {

    }

    else
    {
        wordList patchFieldTypes
        (
            aMesh().boundary().size(),
            fixedValueFaPatchScalarField::typeName
        );

        massFluxPtr_ =
            new areaScalarField
            (
                IOobject
                (
                    "massFlux",
                    DB().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
            aMesh(),
            dimensioned<scalar>("massFlux", dimless, 0),
            patchFieldTypes
            );
    }
}
void freeSurface::makeheatFluxl()
{
    if (debug)
    {
        Info<< "freeSurface::makeheatFluxl() : "
            << "making interfacial temperature"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (heatFluxlPtr_)
    {
        FatalErrorIn("freeSurface::makeheatFluxl()")
            << "heatFluxl already exists"
            << abort(FatalError);
    }

    IOobject heatFluxlHeader
    (
        "heatFluxl",
        DB().timeName(),
        mesh(),
        IOobject::MUST_READ
    );


    if (heatFluxlHeader.headerOk())
    {

    }

    else
    {
        wordList patchFieldTypes
        (
            aMesh().boundary().size(),
            fixedValueFaPatchScalarField::typeName
        );

        heatFluxlPtr_ =
            new areaScalarField
            (
                IOobject
                (
                    "heatFluxl",
                    DB().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
            aMesh(),
            dimensioned<scalar>("heatFluxl", dimless, 0),
            patchFieldTypes
            );
    }
}

void freeSurface::makeheatFluxv()
{
    if (debug)
    {
        Info<< "freeSurface::makeheatFluxv() : "
            << "making interfacial temperature"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (heatFluxvPtr_)
    {
        FatalErrorIn("freeSurface::makeheatFluxv()")
            << "heatFluxv already exists"
            << abort(FatalError);
    }

    IOobject heatFluxvHeader
    (
        "heatFluxv",
        DB().timeName(),
        mesh(),
        IOobject::MUST_READ
    );


    if (heatFluxvHeader.headerOk())
    {

    }

    else
    {
        wordList patchFieldTypes
        (
            aMesh().boundary().size(),
            fixedValueFaPatchScalarField::typeName
        );

        heatFluxvPtr_ =
            new areaScalarField
            (
                IOobject
                (
                    "heatFluxv",
                    DB().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
            aMesh(),
            dimensioned<scalar>("heatFluxv", dimless, 0),
            patchFieldTypes
            );
    }
}

void freeSurface::makeUi()
{
    if (debug)
    {
        Info<< "freeSurface::makeUi() : "
            << "making free-surface velocity field"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (UiPtr_)
    {
        FatalErrorIn("freeSurface::makeUi()")
            << "free-surface velocity field already exists"
            << abort(FatalError);
    }

    wordList patchFieldTypes
    (
        aMesh().boundary().size(),
        zeroGradientFaPatchVectorField::typeName
//	fixedValueFaPatchVectorField::typeName
    );

    UiPtr_ = new areaVectorField
    (
        IOobject
        (
            "Ui",
            DB().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh(),
        dimensioned<vector>("Ui", dimless, vector::zero),
        patchFieldTypes
    );

}


void freeSurface::makeconci()
{
    if (debug)
    {
        Info<< "freeSurface::makeconci() : "
            << "making interfacial temperature"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (conciPtr_)
    {
        FatalErrorIn("freeSurface::makeconci()")
            << "conci already exists"
            << abort(FatalError);
    }

    IOobject conciHeader
    (
        "conci",
        DB().timeName(),
        mesh(),
        IOobject::MUST_READ
    );


    if (conciHeader.headerOk())
    {
/*
        conciPtr_ =
            new areaScalarField
            (
                IOobject
                (
                    "Ti",
                    DB().timeName(),
                    mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            );
*/
    }

    else
    {

        wordList patchFieldTypes
        (
            aMesh().boundary().size(),
//            zeroGradientFaPatchVectorField::typeName
           calculatedFaPatchVectorField::typeName
        );
/*
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
*/
        conciPtr_ =
            new areaScalarField
            (
                IOobject
                (
                    "conci",
                    DB().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
            aMesh(),
            dimensioned<scalar>("conci", dimless, 0),
            patchFieldTypes
//                aMesh().areaCentres().internalField()
            );
    }
}

void freeSurface::makeMotionPointsMask()
{
    if (debug)
    {
        Info<< "freeSurface::makeMotionPointsMask() : "
            << "making motion points mask"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (motionPointsMaskPtr_)
    {
        FatalErrorIn("freeSurface::motionPointsMask()")
            << "motion points mask already exists"
            << abort(FatalError);
    }


    if(aPatchID() == -1)
    {
        FatalErrorIn("freeSurface::makeMotionPointsMask()")
            << "Free surface patch A not defined."
            << abort(FatalError);
    }


    motionPointsMaskPtr_ = new labelList
    (
        mesh().boundaryMesh()[aPatchID()].nPoints(),
        1
    );
}


void freeSurface::makeDirections()
{
    if (debug)
    {
        Info<< "freeSurface::makeDirections() : "
            << "making displacement directions for points and "
            << "control points"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if 
    (
        pointsDisplacementDirPtr_ ||  
        facesDisplacementDirPtr_
    )
    {
        FatalErrorIn("freeSurface::makeDirections()")
            << "points and control points displacement directions "
            << "already exists"
            << abort(FatalError);
    }


    if(aPatchID() == -1)
    {
        FatalErrorIn("freeSurface::makeDirections()")
            << "Free surface patch A not defined."
            << abort(FatalError);
    }


    pointsDisplacementDirPtr_ = 
        new vectorField
        (
            mesh().boundaryMesh()[aPatchID()].nPoints(),
            vector::zero
        );

    facesDisplacementDirPtr_ = 
        new vectorField
        (
            mesh().boundaryMesh()[aPatchID()].size(),
            vector::zero
        );

    if(!normalMotionDir())
    {
        if(mag(g_).value() <= SMALL)
        {
            FatalErrorIn("freeSurface::makeDirections()")
                << "Zero gravity"
                    << abort(FatalError);
        }

        facesDisplacementDir() = -(g_/mag(g_)).value();
        pointsDisplacementDir() = -(g_/mag(g_)).value();
    }

    updateDisplacementDirections();
}


void freeSurface::makeTotalDisplacement()
{
    if (debug)
    {
        Info<< "freeSurface::makeTotalDisplacement() : "
            << "making zero total points displacement"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (totalDisplacementPtr_)
    {
        FatalErrorIn("freeSurface::makeTotalDisplacement()")
            << "total points displacement already exists"
            << abort(FatalError);
    }

    totalDisplacementPtr_ =
        new vectorIOField
        (
            IOobject
            (
                "totalDisplacement",
                DB().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            vectorField
            (
                mesh().boundaryMesh()[aPatchID()].nPoints(), 
                vector::zero
            )
        );
}
 

void freeSurface::readTotalDisplacement()
{
    if (debug)
    {
        Info<< "freeSurface::readTotalDisplacement() : "
            << "reading total points displacement if present"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (totalDisplacementPtr_)
    {
        FatalErrorIn("freeSurface::makeTotalDisplacement()")
            << "total points displacement already exists"
            << abort(FatalError);
    }

    if
    (
        IOobject
        (
            "totalDisplacement",
            DB().timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ).headerOk()
    )
    {
        totalDisplacementPtr_ =
            new vectorIOField
            (
                IOobject
                (
                    "totalDisplacement",
                    DB().timeName(),
                    mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )     
            );
    }           
}


void freeSurface::makeFaMesh() const
{
    if (debug)
    {
        Info<< "freeSurface::makeFaMesh() : "
            << "making finite area mesh"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (aMeshPtr_)
    {
        FatalErrorIn("freeSurface::makeFaMesh()")
            << "finite area mesh already exists"
            << abort(FatalError);
    }


    aMeshPtr_ = new faMesh(mesh());
}


// void freeSurface::makeMeshMotionSolver()
// {
//     if (debug)
//     {
//         Info<< "freeSurface::makeMeshMotionSolver() : "
//             << "making mesh motion solver"
//             << endl;
//     }


//     // It is an error to attempt to recalculate
//     // if the pointer is already set
//     if (mSolverPtr_)
//     {
//         FatalErrorIn("freeSurface::makeMeshMotionSolver()")
//             << "mesh motion solver already exists"
//             << abort(FatalError);
//     }


//     mSolverPtr_ = dynamic_cast<tetDecompositionMotionSolver*>
//     (
//         motionSolver::New(mesh()).ptr()
//     );
// }


void freeSurface::makeUs() const
{
    if (debug)
    {
        Info<< "freeSurface::makeUs() : "
            << "making free-surface velocity field"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (UsPtr_)
    {
        FatalErrorIn("freeSurface::makeUs()")
            << "free-surface velocity field already exists"
            << abort(FatalError);
    }


    wordList patchFieldTypes
    (
        aMesh().boundary().size(),
        zeroGradientFaPatchVectorField::typeName
    );

    forAll(aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type()
         == wedgeFaPatch::typeName
        )
        {
            patchFieldTypes[patchI] = 
                wedgeFaPatchVectorField::typeName;
        }
        else
        {
            label ngbPolyPatchID = 
                aMesh().boundary()[patchI].ngbPolyPatchIndex();

            if (ngbPolyPatchID != -1)
            {
                if
                (
                    mesh().boundary()[ngbPolyPatchID].type()
                 == wallFvPatch::typeName
                )
                {
                    patchFieldTypes[patchI] =
                        slipFaPatchVectorField::typeName;
                }
            }
        }
    }


    UsPtr_ = new areaVectorField
    (
        IOobject
        (
            "Us",
            DB().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aMesh(),
        dimensioned<vector>("Us", dimVelocity, vector::zero),
        patchFieldTypes
    );
}


void freeSurface::makePhis()
{
    if (debug)
    {
        Info<< "freeSurface::makePhis() : "
            << "making free-surface fluid flux"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (phisPtr_)
    {
        FatalErrorIn("freeSurface::makePhis()")
            << "free-surface fluid flux already exists"
            << abort(FatalError);
    }


    phisPtr_ = new edgeScalarField
    (
        IOobject
        (
            "phis",
            DB().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearEdgeInterpolate(Us()) & aMesh().Le()
    );
}


void freeSurface::makePhi()
{
    if (debug)
    {
        Info<< "freeSurface::makePhi() : "
            << "making free-surface flux"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (phiPtr_)
    {
        FatalErrorIn("freeSurface::makePhis()")
            << "free-surface flux already exists"
            << abort(FatalError);
    }

    phiPtr_ = new areaScalarField
    (
        IOobject
        (
            "faPhi",
            DB().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aMesh(),
        dimensionedScalar("faPhi", dimVolume/dimTime, 0),
        zeroGradientFaPatchScalarField::typeName
    );
    phiPtr_->correctBoundaryConditions();
}


void freeSurface::makeDdtPhi()
{
    if (debug)
    {
        Info<< "freeSurface::makeDdtPhi() : "
            << "making free-surface flux time derivative"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ddtPhiPtr_)
    {
        FatalErrorIn("freeSurface::makePhis()")
            << "free-surface flux time derivative already exists"
            << abort(FatalError);
    }

    ddtPhiPtr_ = new areaScalarField
    (
        IOobject
        (
            "ddtPhi",
            DB().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fac::ddt(phi())
    );
}

void freeSurface::makeSurfactConc() const
{
    if (debug)
    {
        Info<< "freeSurface::makeSurfactConc() : "
            << "making free-surface surfactant concentration field"
            << endl;
    }


    // It is an error to attempt to recalculate    
    // if the pointer is already set
    if (surfactConcPtr_)
    {
        FatalErrorIn("freeSurface::makeUs()")
            << "free-surface surfactant concentratio field already exists"
            << abort(FatalError);
    }

    surfactConcPtr_ = new areaScalarField
    (
        IOobject
        (
            "Cs",
            DB().timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh()
    );
}


void freeSurface::makeSurfaceTension() const
{
    if (debug)
    {
        Info<< "freeSurface::makeSurfaceTension() : "
            << "making surface tension field"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (surfaceTensionPtr_)
    {
        FatalErrorIn("freeSurface::makeSurfaceTension()")
            << "surface tension field already exists"
            << abort(FatalError);
    }


    surfaceTensionPtr_ = new areaScalarField
    (
        IOobject
        (
            "surfaceTension",
            DB().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cleanInterfaceSurfTension()
      + surfactant().surfactR()*
        surfactant().surfactT()*
        surfactant().surfactSaturatedConc()*
        log(1.0 - surfactantConcentration()/
        surfactant().surfactSaturatedConc())
    );
}

void freeSurface::makecontactAngle() const
{
    if (debug)
    {
        Info<< "freeSurface::makecontactAngle() : "
            << "making contact angle field"
            << endl;
    }

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
    }
}

void freeSurface::makeSurfactant() const
{
    if (debug)
    {
        Info<< "freeSurface::makeSurfactant() : "
            << "making surfactant properties"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (surfactantPtr_)
    {
        FatalErrorIn("freeSurface::makeSurfactant()")
            << "surfactant properties already exists"
            << abort(FatalError);
    }


    const dictionary& surfactProp = 
        this->subDict("surfactantProperties");

    surfactantPtr_ = new surfactantProperties(surfactProp);
}


void freeSurface::makeFluidIndicator()
{
    if (debug)
    {
        Info<< "freeSurface::makeFluidIndicator() : "
            << "making fluid indicator"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (fluidIndicatorPtr_)
    {
        FatalErrorIn("freeSurface::makeFluidIndicator()")
            << "fluid indicator already exists"
            << abort(FatalError);
    }

    fluidIndicatorPtr_ = new volScalarField
    (
        IOobject
        (
            "fluidIndicator",
            DB().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("1", dimless, 1.0),
        zeroGradientFvPatchScalarField::typeName
    );

    volScalarField& fluidIndicator = *fluidIndicatorPtr_;

    if (twoFluids())
    {
        // find start cell
        label pointOnShadowPatch =
            mesh().boundaryMesh()[bPatchID()][0][0];

        label startCell = mesh().pointCells()[pointOnShadowPatch][0];


        // get cell-cells addressing
        const labelListList& cellCells = mesh().cellCells();

        SLList<label> slList(startCell);

        while (slList.size())
        {
            label curCell = slList.removeHead();

            if (fluidIndicator[curCell] == 1)
            {
                fluidIndicator[curCell] = 0.0;

                for (int i = 0; i < cellCells[curCell].size(); i++)
                {
                    slList.append(cellCells[curCell][i]);
                }
            }
        }
    }

    fluidIndicator.correctBoundaryConditions();
}

void freeSurface::makeNGradUn() const
{
    if (debug)
    {
        Info<< "freeSurface::makeNGradUn() : "
            << "making free-surface normal derivative of normal velocity"
            << endl;
    }


    // It is an error to attempt to recalculate    
    // if the pointer is already set
    if (nGradUnPtr_)
    {
        FatalErrorIn("freeSurface::makeNGradUn()")
            << "free-surface normal derivative of normal velocity "
                << "field already exists"
                << abort(FatalError);
    }

    nGradUnPtr_ = new scalarField
    (
        aMesh().nFaces(),
        0
    );
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const IOpatchToPatchInterpolation& freeSurface::interpolatorAB()
{
    if (!interpolatorABPtr_)
    {
        makeInterpolators();
    }
    
    return *interpolatorABPtr_;
}


const IOpatchToPatchInterpolation& freeSurface::interpolatorBA()
{
    if (!interpolatorBAPtr_)
    {
        makeInterpolators();
    }
    
    return *interpolatorBAPtr_;
}


vectorField& freeSurface::controlPoints()
{
    if (!controlPointsPtr_)
    {
        makeControlPoints();
    }

    return *controlPointsPtr_;
}

areaScalarField& freeSurface::Ti()
{
    if (!TiPtr_)
    {
        makeTi();
    }

    return *TiPtr_;
}

areaScalarField& freeSurface::TiDiff()
{
    if (!TiDiffPtr_)
    {
        makeTiDiff();
    }

    return *TiDiffPtr_;
}

areaScalarField& freeSurface::dTi()
{
    if (!dTiPtr_)
    {
        makedTi();
    }

    return *dTiPtr_;
}

areaVectorField& freeSurface::gradTi()
{
    if (!gradTiPtr_)
    {
        makegradTi();
    }

    return *gradTiPtr_;
}

areaScalarField& freeSurface::massFlux()
{
    if (!massFluxPtr_)
    {
        makemassFlux();
    }

    return *massFluxPtr_;
}

areaScalarField& freeSurface::heatFluxl()
{
    if (!heatFluxlPtr_)
    {
        makeheatFluxl();
    }

    return *heatFluxlPtr_;
}

areaScalarField& freeSurface::heatFluxv()
{
    if (!heatFluxvPtr_)
    {
        makeheatFluxv();
    }

    return *heatFluxvPtr_;
}
areaVectorField& freeSurface::Ui()
{
    if (!UiPtr_)
    {
        makeUi();
    }

    return *UiPtr_;

}
areaScalarField& freeSurface::conci()
{
    if (!conciPtr_)
    {
        makeconci();
    }

    return *conciPtr_;
}

labelList& freeSurface::motionPointsMask()
{
    if (!motionPointsMaskPtr_)
    {
        makeMotionPointsMask();
    }

    return *motionPointsMaskPtr_;
}


vectorField& freeSurface::pointsDisplacementDir()
{
    if (!pointsDisplacementDirPtr_)
    {
        makeDirections();
    }

    return *pointsDisplacementDirPtr_;
}


vectorField& freeSurface::facesDisplacementDir()
{
    if (!facesDisplacementDirPtr_)
    {
        makeDirections();
    }

    return *facesDisplacementDirPtr_;
}


vectorField& freeSurface::totalDisplacement()
{
    if (!totalDisplacementPtr_)
    {
        makeTotalDisplacement();
    }

    return *totalDisplacementPtr_;
}


faMesh& freeSurface::aMesh()
{
    if (!aMeshPtr_)
    {
        makeFaMesh();
    }
    
    return *aMeshPtr_;
}

const faMesh& freeSurface::aMesh() const
{
    if (!aMeshPtr_)
    {
        makeFaMesh();
    }
    
    return *aMeshPtr_;
}

// tetDecompositionMotionSolver& freeSurface::meshMotionSolver()
// {
//     if (!mSolverPtr_)
//     {
//         makeMeshMotionSolver();
//     }

//     return *mSolverPtr_;
// }


areaVectorField& freeSurface::Us()
{
    if (!UsPtr_)
    {
        makeUs();
    }
    
    return *UsPtr_;
}

const areaVectorField& freeSurface::Us() const
{
    if (!UsPtr_)
    {
        makeUs();
    }
    
    return *UsPtr_;
}

edgeScalarField& freeSurface::Phis()
{
    if (!phisPtr_)
    {
        makePhis();
    }
    
    return *phisPtr_;
}

areaScalarField& freeSurface::phi()
{
    if (!phiPtr_)
    {
        makePhi();
    }

    phiPtr_->internalField() = phi_.boundaryField()[aPatchID()];
    phiPtr_->correctBoundaryConditions();

    return *phiPtr_;
}

areaScalarField& freeSurface::ddtPhi()
{
    if (!ddtPhiPtr_)
    {
        makeDdtPhi();
    }

    return *ddtPhiPtr_;
}

areaScalarField& freeSurface::surfactantConcentration()
{
    if (!surfactConcPtr_)
    {
        makeSurfactConc();
    }
    
    return *surfactConcPtr_;
}

const areaScalarField& freeSurface::surfactantConcentration() const
{
    if (!surfactConcPtr_)
    {
        makeSurfactConc();
    }
    
    return *surfactConcPtr_;
}

areaScalarField& freeSurface::surfaceTension()
{
    if (!surfaceTensionPtr_)
    {
        makeSurfaceTension();
    }
    
    return *surfaceTensionPtr_;
}

const areaScalarField& freeSurface::surfaceTension() const
{
    if (!surfaceTensionPtr_)
    {
        makeSurfaceTension();
    }
    
    return *surfaceTensionPtr_;
}

edgeScalarField& freeSurface::contactAngle()
{
    if (!contactAnglePtr_)
    {
        makecontactAngle();
    }
    
    return *contactAnglePtr_;
}
const edgeScalarField& freeSurface::contactAngle() const
{
    if (!contactAnglePtr_)
    {
        makecontactAngle();
    }
    
    return *contactAnglePtr_;
}

const surfactantProperties& freeSurface::surfactant() const
{
    if (!surfactantPtr_)
    {
        makeSurfactant();
    }
    
    return *surfactantPtr_;
}


const volScalarField& freeSurface::fluidIndicator()
{    
    if (!fluidIndicatorPtr_)
    {
        makeFluidIndicator();
    }

    return *fluidIndicatorPtr_;
}


tmp<areaVectorField> freeSurface::surfaceTensionGrad()
{
    tmp<areaVectorField> tgrad
    (
        new areaVectorField
        (
            IOobject
            (
                "surfaceTensionGrad",
                DB().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            (-fac::grad(surfactantConcentration())*
            surfactant().surfactR()*surfactant().surfactT()/
            (1.0 - surfactantConcentration()/
            surfactant().surfactSaturatedConc()))()
        )
    );
    
    return tgrad;
}

scalarField& freeSurface::nGradUn()
{
    if (!nGradUnPtr_)
    {
        makeNGradUn();
    }

    return *nGradUnPtr_;
}


const scalarField& freeSurface::nGradUn() const
{    
    if (!nGradUnPtr_)
    {
        makeNGradUn();
    }

    return *nGradUnPtr_;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
