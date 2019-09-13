/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "liquidFilmParallelFvMotionSolver.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(liquidFilmParallelFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        fvMotionSolver,
        liquidFilmParallelFvMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liquidFilmParallelFvMotionSolver::liquidFilmParallelFvMotionSolver
(
    const polyMesh& mesh,
    Istream&
)
:
    fvMotionSolver(mesh),
    pointMotionU_
    (
        IOobject
        (
            "pointMotionU",
            fvMesh_.time().timeName(),
            fvMesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pointMesh::New(fvMesh_)
    ),
    pointPoint_(fvMesh_.nPoints(), -1),
    max_(0),
    min_(0),
    motionDirection_(lookup("motionDirection")),
    freeSurfacePatchID_(mesh.boundaryMesh().findPatchID("freeSurface")),
    freeSurfaceZoneID_(mesh.faceZones().findZoneID("freeSurfaceZone")),
    procToGlobalFZmap_()
{
    motionDirection_ /= mag(motionDirection_) + SMALL;

    boundBox box(mesh.points());

    max_ = (motionDirection_&box.max());
    min_ = (motionDirection_&box.min());

    // Check free surface patch 
    if(freeSurfacePatchID_ == -1)
    {
        FatalErrorIn
        (
            "liquidFilmParallelFvMotionSolver(...)"
        )
            << "Free surface patch not defined." << abort(FatalError);
    }

    // Check free surface face zone
    if(freeSurfaceZoneID_ == -1)
    {
        FatalErrorIn
        (
            "liquidFilmParallelFvMotionSolver(...)"
        )
            << "Free surface face zone not defined." << abort(FatalError);
    }

    const vectorField& points = fvMesh_.points();

    vectorField fzLocalPoints = 
        mesh.faceZones()[freeSurfaceZoneID_]().localPoints();

    forAll(pointPoint_, pointI)
    {
        scalar minDist = GREAT;

        forAll(fzLocalPoints, lpI)
        {
            vector R = fzLocalPoints[lpI] - points[pointI];
            R -= motionDirection_*(motionDirection_&R);
            scalar curDist = mag(R);

            if (curDist < minDist)
            {
                minDist = curDist;
                pointPoint_[pointI] = lpI;
            }
        }
    }

    // face zona map

    procToGlobalFZmap_.setSize(fzLocalPoints.size(), -1);

    if(Pstream::parRun())
    {
        vectorField fzGlobalPoints = 
            mesh.faceZones()[freeSurfaceZoneID_]().localPoints();

        //- set all slave points to zero because only the master order is used
        if(!Pstream::master())
        {
            fzGlobalPoints *= 0.0;
        }

        //- pass points to all procs
        reduce(fzGlobalPoints, sumOp<vectorField>());

        //- now every proc has the master's list of FZ points
       	//- every proc must now find the mapping from their local FZ points to
        //- the global FZ points

        vectorField fzLocalPoints = 
            mesh.faceZones()[freeSurfaceZoneID_]().localPoints();

        const edgeList& fzLocalEdges = 
            mesh.faceZones()[freeSurfaceZoneID_]().edges();

        const labelListList& fzPointEdges = 
            mesh.faceZones()[freeSurfaceZoneID_]().pointEdges();

        scalarField minEdgeLength(fzLocalPoints.size(), GREAT);

        forAll(minEdgeLength, pI)
        {
            const labelList& curPointEdges = fzPointEdges[pI];

            forAll(curPointEdges, eI)
            {
                scalar Le = fzLocalEdges[curPointEdges[eI]].mag(fzLocalPoints);
                if (Le < minEdgeLength[pI])
                {
                    minEdgeLength[pI] = Le;
                }
            }
        }

        forAll(fzGlobalPoints, globalPointI)
        {
            forAll(fzLocalPoints, procPointI)
            {
                if
                (
                    mag
                    (
                        fzLocalPoints[procPointI] 
                      - fzGlobalPoints[globalPointI]
                    )
                  < 1e-4*minEdgeLength[procPointI]
                )
                {
                    procToGlobalFZmap_[globalPointI] = procPointI;
                    break;
                }
            }
        }

        forAll(procToGlobalFZmap_, globalPointI)
        {
            if (procToGlobalFZmap_[globalPointI] == -1)
            {
                FatalErrorIn
                (
                    "liquidFilmParallelFvMotionSolver(...)"
                )
                    << "Global FZ map is not correct" 
                        << abort(FatalError);
            }
        }
    }
    else
    {
        forAll(procToGlobalFZmap_, globalPointI)
        {
            procToGlobalFZmap_[globalPointI] = globalPointI;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::liquidFilmParallelFvMotionSolver::~liquidFilmParallelFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::liquidFilmParallelFvMotionSolver::curPoints() const
{
    scalar deltaT = fvMesh_.time().deltaT().value();
    const vectorField& oldPoints = fvMesh_.points();
    vectorField& motionPointUI = pointMotionU_.internalField();

    const labelList& fzMeshPoints =
        fvMesh_.faceZones()[freeSurfaceZoneID_]().meshPoints();

    // Displacement for face zone points
    vectorField fzMotionPointU(fzMeshPoints.size(), vector::zero);

    // Face zone points
    vectorField fzPoints(fzMeshPoints.size(), vector::zero);

    //- Inter-proc points are shared by multiple procs
    //  pointNumProc is the number of procs which a point lies on
    scalarField pointNumProcs(fzMeshPoints.size(), 0);

    forAll(fzMotionPointU, globalPointI)
    {
        label localPoint = procToGlobalFZmap_[globalPointI];
	    
        if(fzMeshPoints[localPoint] < fvMesh_.nPoints())
        {
            label procPoint = fzMeshPoints[localPoint];

            fzMotionPointU[globalPointI] = motionPointUI[procPoint];

            fzPoints[globalPointI] = oldPoints[procPoint];

            pointNumProcs[globalPointI] = 1;
        }
    }

    if (Pstream::parRun())
    {
        reduce(fzMotionPointU, sumOp<vectorField>());
        reduce(fzPoints, sumOp<vectorField>());
        reduce(pointNumProcs, sumOp<scalarField>());

        //- now average the motionPointU between all procs
        fzMotionPointU /= pointNumProcs;
        fzPoints /= pointNumProcs;
    }

    //- The globalFZnewPoints now contains the correct FZ displacement in
    //  a global order, now convert them back into the local proc order

    vectorField fzLocalMotionPointU(fzMeshPoints.size(), vector::zero);
    vectorField fzLocalPoints(fzMeshPoints.size(), vector::zero);

    forAll(fzMotionPointU, globalPointI)
    {
        label localPoint = procToGlobalFZmap_[globalPointI];
	    
        fzLocalMotionPointU[localPoint] = fzMotionPointU[globalPointI];
        fzLocalPoints[localPoint] = fzPoints[globalPointI];
    }

    forAll(motionPointUI, pointI)
    {
        if
        (
            (motionDirection_ & oldPoints[pointI])
         <= (motionDirection_ & fzLocalPoints[pointPoint_[pointI]])
        )
        {
            motionPointUI[pointI] =
                fzLocalMotionPointU[pointPoint_[pointI]]
               *(
                    (
                        (motionDirection_ & oldPoints[pointI])
                      - min_
                    )
                   /(
                        (
                            motionDirection_ 
                          & fzLocalPoints[pointPoint_[pointI]]
                        ) 
                      - min_
                    )
                );
        }
        else
        {
            motionPointUI[pointI] =
                fzLocalMotionPointU[pointPoint_[pointI]]
               *(
                    (
                        max_ 
                      - (motionDirection_ & oldPoints[pointI])
                    )
                   /(
                        max_ 
                      - (
                            motionDirection_ 
                          & fzLocalPoints[pointPoint_[pointI]]
                        )
                    )
                );
        }
    }

    vectorField newPoints = fvMesh_.allPoints();

    for(label i=0; i<fvMesh_.nPoints(); i++)
    {
        newPoints[i] += deltaT*motionPointUI[i];
    }

//     Pout << "end-1" << endl;
//     ::sleep(5);

    forAll(fzLocalMotionPointU, pointI)
    {
        if (fzMeshPoints[pointI] >= fvMesh_.nPoints())
        {
            newPoints[fzMeshPoints[pointI]] += 
                deltaT*fzLocalMotionPointU[pointI];
        }
    }
 
//     Pout << "end-1a" << endl;
//     ::sleep(5);

    twoDCorrectPoints(newPoints);

//     Pout << "end-2" << endl;
//     ::sleep(5);

    tmp<pointField> tcurPoints(new pointField(newPoints));

//     tmp<pointField> tcurPoints
//     (
//         oldPoints + deltaT*pointMotionU_.internalField()
//     );


//     Pout << "end-2a" << endl;
//     ::sleep(5);

//     twoDCorrectPoints(tcurPoints());

//     Pout << "end-3" << endl;
//     ::sleep(5);

    return tcurPoints;
}

// ************************************************************************* //
