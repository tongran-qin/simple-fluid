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
#include "emptyFaPatch.H"
#include "wedgeFaPatch.H"
#include "wallFvPatch.H"
#include "PstreamCombineReduceOps.H"
#include "coordinateSystem.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
//z0312
tmp<vectorField> freeSurface::pointDisplacement(const scalarField& deltaH)
{
    const pointField& points = aMesh().patch().localPoints();
    const labelListList& pointFaces = aMesh().patch().pointFaces();

    controlPoints() += facesDisplacementDir()*deltaH;
/*
    if (correctCurvature_)
    {
        // Correct controPoints next to fixed patches
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
            const vectorField& fCentres =
                aMesh().areaCentres().internalField();
            
            forAll(eFaces, edgeI)
            {
                const label& curFace = eFaces[edgeI];
                const labelList& curFaceFaces = fFaces[curFace];

                scalar H = 0.0;
                label counter = 0;

                forAll(curFaceFaces, faceI)
                {
                    label index = findIndex(eFaces, curFaceFaces[faceI]);
                    
                    if (index == -1)
                    {
                        H +=
                            facesDisplacementDir()[curFaceFaces[faceI]]
                          & (
                                controlPoints()[curFaceFaces[faceI]]
                              - fCentres[curFaceFaces[faceI]]
                            );

                        counter++;
                    }
                }
                
                H /= counter;

                controlPoints()[curFace] =
                    fCentres[curFace]
                  + facesDisplacementDir()[curFace]*H;
            }
        }
    }
*/

    tmp<vectorField> tdisplacement
    (
        new vectorField
        (
            points.size(),
            vector::zero
        )
    );
    
    vectorField& displacement = tdisplacement();


    // Calculate displacement of internal points
    const vectorField& pointNormals = aMesh().pointAreaNormals();
    const edgeList& edges = aMesh().patch().edges();
    labelList internalPoints = aMesh().internalPoints();

    forAll (internalPoints, pointI)
    {
        label curPoint = internalPoints[pointI];

        const labelList& curPointFaces = pointFaces[curPoint];

        vectorField lsPoints(curPointFaces.size(), vector::zero);

        for (label i=0; i<curPointFaces.size(); i++)
        {
            label curFace = curPointFaces[i];
            
            lsPoints[i] = controlPoints()[curFace];
        }
       
        vectorField pointAndNormal = 
            lsPlanePointAndNormal
            (
                lsPoints, 
                points[curPoint], 
                pointNormals[curPoint]
            );

        vector& P = pointAndNormal[0];
        vector& N = pointAndNormal[1];

        displacement[curPoint] = 
            pointsDisplacementDir()[curPoint]
           *((P - points[curPoint])&N)
           /(pointsDisplacementDir()[curPoint]&N);
    }


    // Mirror control points
    FieldField<Field, vector> patchMirrorPoints(aMesh().boundary().size());

    // Old faMesh points
    vectorField oldPoints(aMesh().nPoints(), vector::zero);
    const labelList& meshPoints = aMesh().patch().meshPoints();
    forAll(oldPoints, pI)
    {
        oldPoints[pI] = 
            mesh().oldPoints()[meshPoints[pI]];
    }

    forAll(patchMirrorPoints, patchI)
    {
        patchMirrorPoints.set
        (
            patchI,
            new vectorField
            (
                aMesh().boundary()[patchI].faPatch::size(), 
                vector::zero
            )
        );

        vectorField N = 
            aMesh().boundary()[patchI].ngbPolyPatchFaceNormals();

        const labelList& eFaces =
            aMesh().boundary()[patchI].edgeFaces();

        // Correct N according to specified contact angle
        if (contactAnglePtr_)
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
                    scalar rotAngle = 
                        average
                        (
                            90
                          - contactAnglePtr_->boundaryField()[patchI]
                        );
                    rotAngle *= M_PI/180.0;

                    vectorField rotationAxis(N.size(), vector::zero);

                    const vectorField& pEdgN =
                        aMesh().edgeAreaNormals().boundaryField()[patchI];

                    rotationAxis = (N^pEdgN);

                    const edgeList::subList patchEdges =
                        aMesh().boundary()[patchI].patchSlice(aMesh().edges());

                    forAll(rotationAxis, edgeI)
                    {
                        vector e = patchEdges[edgeI].vec(oldPoints);
//                         vector e = patchEdges[edgeI].vec(aMesh().points());

                        // Adjust direction
                        rotationAxis[edgeI] = 
                            e*(e&rotationAxis[edgeI])
                           /mag((e&rotationAxis[edgeI]));
                    }
                    rotationAxis /= mag(rotationAxis) + SMALL;

                    vectorField rotationAxis2 = rotationAxis;
                    forAll(rotationAxis2, edgeI)
                    {
                        rotationAxis2[edgeI] = 
                            (N[edgeI]^facesDisplacementDir()[eFaces[edgeI]]);
                        
                        // Adjust direction
                        rotationAxis2[edgeI] = 
                            rotationAxis2[edgeI]
                           *(rotationAxis2[edgeI]&rotationAxis[edgeI])
                           /mag((rotationAxis2[edgeI]&rotationAxis[edgeI]));
                    }
                    rotationAxis2 /= mag(rotationAxis2) + SMALL;

                    // Rodrigues' rotation formula
                    N = N*cos(rotAngle)
                      + rotationAxis*(rotationAxis & N)*(1 - cos(rotAngle))
                      + (rotationAxis^N)*sin(rotAngle);

                    N /= mag(N);

                    N = (rotationAxis^N);
                    
                    N = (N^rotationAxis2);

                    N /= mag(N);
                }
            }
        }

        const labelList peFaces = 
            labelList::subList
            (
                aMesh().edgeOwner(),
                aMesh().boundary()[patchI].faPatch::size(),
                aMesh().boundary()[patchI].start()
            );

        const labelList& pEdges = aMesh().boundary()[patchI];

        vectorField peCentres(pEdges.size(), vector::zero);
        forAll(peCentres, edgeI)
        {
            peCentres[edgeI] = 
                edges[pEdges[edgeI]].centre(points);
        }

        vectorField delta =
            vectorField(controlPoints(), peFaces)
          - peCentres;

        patchMirrorPoints[patchI] = 
            peCentres + ((I - 2*N*N)&delta);
    }


    // Calculate displacement of boundary points 
    labelList boundaryPoints = aMesh().boundaryPoints();

    const labelListList& edgeFaces = aMesh().patch().edgeFaces();
    const labelListList& pointEdges = aMesh().patch().pointEdges();

    forAll (boundaryPoints, pointI)
    {
        label curPoint = boundaryPoints[pointI];
        
        if (motionPointsMask()[curPoint] == 1)
        {
            // Calculating mirror points
            const labelList& curPointEdges = pointEdges[curPoint];

            vectorField mirrorPoints(2, vector::zero);

            label counter = -1;
        
            forAll (curPointEdges, edgeI)
            {
                label curEdge = curPointEdges[edgeI];

                if(edgeFaces[curEdge].size() == 1)
                {
                    label patchID = -1;
                    label edgeID = -1;
                    forAll(aMesh().boundary(), patchI)
                    {
                        const labelList& pEdges =
                            aMesh().boundary()[patchI];
                        label index = findIndex(pEdges, curEdge);
                        if (index != -1)
                        {
                            patchID = patchI;
                            edgeID = index;
                            break;
                        }
                    }

                    mirrorPoints[++counter] = 
                        patchMirrorPoints[patchID][edgeID];
                }
            }

            // Calculating LS plane fit
            const labelList& curPointFaces = pointFaces[curPoint];
        
            vectorField lsPoints
            (
                curPointFaces.size() + mirrorPoints.size(), 
                vector::zero
            );
                
            counter = -1;

            for (label i=0; i<curPointFaces.size(); i++)
            {
                label curFace = curPointFaces[i];

                lsPoints[++counter] = controlPoints()[curFace];
            }

            for (label i=0; i<mirrorPoints.size(); i++)
            {
                lsPoints[++counter] = mirrorPoints[i];
            }

            vectorField pointAndNormal = 
                lsPlanePointAndNormal
                (
                    lsPoints, 
                    points[curPoint], 
                    pointNormals[curPoint]
                );

            vector& P = pointAndNormal[0];
            vector& N = pointAndNormal[1];

            displacement[curPoint] = 
                pointsDisplacementDir()[curPoint]
               *((P - points[curPoint])&N)
               /(pointsDisplacementDir()[curPoint]&N);
        }
    }


    // Calculate displacement of axis point
    forAll (aMesh().boundary(), patchI)
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
                label axisPoint = wedgePatch.axisPoint();
                
                displacement[axisPoint] =
                    pointsDisplacementDir()[axisPoint]
                   *(
                        pointsDisplacementDir()[axisPoint]
                       &(
                            controlPoints()[pointFaces[axisPoint][0]]
                          - points[axisPoint]
                        )
                    );
            }
        }
    }


    // Calculate displacement of processor patch points
    forAll (aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type() 
         == processorFaPatch::typeName
        )
        {
            const processorFaPatch& procPatch =
                refCast<const processorFaPatch>(aMesh().boundary()[patchI]);

            const labelList& patchPointLabels =
                procPatch.pointLabels();

            FieldField<Field, vector> lsPoints(patchPointLabels.size());
            forAll(lsPoints, pointI)
            {
                lsPoints.set(pointI, new vectorField(0, vector::zero));
            }

            const labelList& nonGlobalPatchPoints = 
                procPatch.nonGlobalPatchPoints();

            forAll(nonGlobalPatchPoints, pointI)
            {
                label curPatchPoint =
                    nonGlobalPatchPoints[pointI];

                label curPoint = 
                    patchPointLabels[curPatchPoint];

                const labelList& curPointFaces = pointFaces[curPoint];

                lsPoints[curPatchPoint].setSize(curPointFaces.size());

                forAll(curPointFaces, faceI)
                {
                    label curFace = curPointFaces[faceI];

                    lsPoints[curPatchPoint][faceI] = controlPoints()[curFace];
                }

#               include "boundaryProcessorFaPatchPoints.H"
            }

            scalar lsPointsSize = 0;
            forAll(lsPoints, pointI)
            {
                lsPointsSize +=
                    2*lsPoints[pointI].size()*sizeof(vector);
            }
            
            // Parallel data exchange
            {
                OPstream toNeighbProc
                (
                    Pstream::blocking,
                    procPatch.neighbProcNo(), 
                    lsPointsSize
                );

                toNeighbProc << lsPoints;
            }

            FieldField<Field, vector> ngbLsPoints(patchPointLabels.size());

            {
                IPstream fromNeighbProc
                (
                    Pstream::blocking,
                    procPatch.neighbProcNo(),
                    lsPointsSize
                );

                fromNeighbProc >> ngbLsPoints;
            }

            forAll(nonGlobalPatchPoints, pointI)
            {
                label curPatchPoint = 
                    nonGlobalPatchPoints[pointI];

                label curPoint = 
                    patchPointLabels[curPatchPoint];

                label curNgbPoint = procPatch.neighbPoints()[curPatchPoint];

                vectorField allLsPoints
                (
                    lsPoints[curPatchPoint].size()
                  + ngbLsPoints[curNgbPoint].size(),
                    vector::zero
                );
                
                label counter = -1;
                forAll(lsPoints[curPatchPoint], pointI)
                {
                    allLsPoints[++counter] = lsPoints[curPatchPoint][pointI];
                }
                forAll(ngbLsPoints[curNgbPoint], pointI)
                {
                    allLsPoints[++counter] = ngbLsPoints[curNgbPoint][pointI];
                }

                vectorField pointAndNormal = 
                    lsPlanePointAndNormal
                    (
                        allLsPoints, 
                        points[curPoint], 
                        pointNormals[curPoint]
                    );

                vector& P = pointAndNormal[0];
                vector& N = pointAndNormal[1];
                    
                if (motionPointsMask()[curPoint] != 0)
                {
                    displacement[curPoint] = 
                        pointsDisplacementDir()[curPoint]
                       *((P - points[curPoint])&N)
                       /(pointsDisplacementDir()[curPoint]&N);
                }
            }
        }
    }


    // Calculate displacement of global processor patch points
    if (aMesh().globalData().nGlobalPoints() > 0)
    {
        const labelList& spLabels =
            aMesh().globalData().sharedPointLabels();

        const labelList& addr = aMesh().globalData().sharedPointAddr();

        for (label k=0; k<aMesh().globalData().nGlobalPoints(); k++)
        {
            List<List<vector> > procLsPoints(Pstream::nProcs());

            label curSharedPointIndex = findIndex(addr, k);

            if (curSharedPointIndex != -1)
            {
                label curPoint = spLabels[curSharedPointIndex];

                const labelList& curPointFaces = pointFaces[curPoint];

                procLsPoints[Pstream::myProcNo()] = 
                    List<vector>(curPointFaces.size());

                forAll (curPointFaces, faceI)
                {
                    label curFace = curPointFaces[faceI];

                    procLsPoints[Pstream::myProcNo()][faceI] = 
                        controlPoints()[curFace];
                }
            }

            Pstream::gatherList(procLsPoints);
            Pstream::scatterList(procLsPoints);
                
            if (curSharedPointIndex != -1)
            {
                label curPoint = spLabels[curSharedPointIndex];
                
                label nAllPoints = 0;
                forAll(procLsPoints, procI)
                {
                    nAllPoints += procLsPoints[procI].size();
                }

                vectorField allPoints(nAllPoints, vector::zero);

                label counter = 0;
                forAll(procLsPoints, procI)
                {
                    forAll(procLsPoints[procI], pointI)
                    {
                        allPoints[counter++] =
                            procLsPoints[procI][pointI];
                    }
                }

                vectorField pointAndNormal = 
                    lsPlanePointAndNormal
                    (
                        allPoints, 
                        points[curPoint], 
                        pointNormals[curPoint]
                    );

                const vector& P = pointAndNormal[0];
                const vector& N = pointAndNormal[1];

                displacement[curPoint] = 
                    pointsDisplacementDir()[curPoint]
                   *((P - points[curPoint])&N)
                   /(pointsDisplacementDir()[curPoint]&N);
            }
        }
    }

    return tdisplacement;
}


/*
tmp<vectorField> freeSurface::pointDisplacement(const scalarField& deltaH)
{
    const pointField& points = aMesh().patch().localPoints();
    const labelListList& pointFaces = aMesh().patch().pointFaces();

    controlPoints() += facesDisplacementDir()*deltaH;

    tmp<vectorField> tdisplacement
    (
        new vectorField
        (
            points.size(),
            vector::zero
        )
    );
    
    vectorField& displacement = tdisplacement();


    // Calculate displacement of internal points
    const vectorField& pointNormals = aMesh().pointAreaNormals();
    const edgeList& edges = aMesh().patch().edges();
    labelList internalPoints = aMesh().internalPoints();

    forAll (internalPoints, pointI)
    {
        label curPoint = internalPoints[pointI];

        const labelList& curPointFaces = pointFaces[curPoint];

        vectorField lsPoints(curPointFaces.size(), vector::zero);

        for (label i=0; i<curPointFaces.size(); i++)
        {
            label curFace = curPointFaces[i];
            
            lsPoints[i] = controlPoints()[curFace];
        }
       
        vectorField pointAndNormal = 
            lsPlanePointAndNormal
            (
                lsPoints, 
                points[curPoint], 
                pointNormals[curPoint]
            );

        vector& P = pointAndNormal[0];
        vector& N = pointAndNormal[1];

        displacement[curPoint] = 
            pointsDisplacementDir()[curPoint]
           *((P - points[curPoint])&N)
           /(pointsDisplacementDir()[curPoint]&N);
    }


    // Mirror control points
    FieldField<Field, vector> patchMirrorPoints(aMesh().boundary().size());

    forAll(patchMirrorPoints, patchI)
    {
        patchMirrorPoints.set
        (
            patchI,
            new vectorField
            (
                aMesh().boundary()[patchI].faPatch::size(), 
                vector::zero
            )
        );

        vectorField N = aMesh().boundary()[patchI].ngbPolyPatchFaceNormals();

        // Correct N according to specified contact angle
        if (contactAnglePtr_)
        {
            label ngbPolyPatchID = 
                aMesh().boundary()[patchI].ngbPolyPatchIndex();

            vectorField Nf = aMesh().faceAreaNormals().internalField();

            if (ngbPolyPatchID != -1)
            {
                if
                (
                    mesh().boundary()[ngbPolyPatchID].type()
                 == wallFvPatch::typeName
                )
                {
                    const labelList peFaces = 
                       labelList::subList
                       (
                          aMesh().edgeOwner(),
                          aMesh().boundary()[patchI].faPatch::size(),
                          aMesh().boundary()[patchI].start()
                       );

//Info<< "ngbPolyPatchID: " << ngbPolyPatchID << endl;
//Info<< "patchI: " << patchI << endl;
//Info<< "peFaces "<<endl <<peFaces <<endl ;

                    vectorField nf = vectorField(Nf, peFaces);

//Info<< "nf: " << nf << endl;

                    scalar rotAngle = 
                        average
                        (
                            mag(contactAnglePtr_->boundaryField()[patchI]) 
                          - 90
                        )*M_PI/180.0;
                    
		    scalar sA = sin(rotAngle);
		    scalar cA = cos(rotAngle);
                    
                    vectorField rotationAxis;

//Info<< "N(-) " << N << endl;

                    rotationAxis = - N^nf;
                    rotationAxis /= mag(rotationAxis) + SMALL;

//Info<< "rotationAxis: " << rotationAxis << endl;

                        // Rodrigues' rotation formula
                    N = N*cA
                      + rotationAxis*(rotationAxis&N)*(1-cA)
                      + (rotationAxis^N)*sA;
                    N /= mag(N);
//Info<< "N(+) " << N << endl;
                }
            }
        }

        const labelList peFaces = 
            labelList::subList
            (
                aMesh().edgeOwner(),
                aMesh().boundary()[patchI].faPatch::size(),
                aMesh().boundary()[patchI].start()
            );

        const labelList& pEdges = aMesh().boundary()[patchI];

        vectorField peCentres(pEdges.size(), vector::zero);
        forAll(peCentres, edgeI)
        {
            peCentres[edgeI] = 
                edges[pEdges[edgeI]].centre(points);
        }

        vectorField delta =
            vectorField(controlPoints(), peFaces)
          - peCentres;

        patchMirrorPoints[patchI] = 
            peCentres + ((I - 2*N*N)&delta);
    }


    // Calculate displacement of boundary points 
    labelList boundaryPoints = aMesh().boundaryPoints();

    const labelListList& edgeFaces = aMesh().patch().edgeFaces();
    const labelListList& pointEdges = aMesh().patch().pointEdges();

    forAll (boundaryPoints, pointI)
    {
        label curPoint = boundaryPoints[pointI];
        
        if (motionPointsMask()[curPoint] == 1)
        {
            // Calculating mirror points
            const labelList& curPointEdges = pointEdges[curPoint];

            vectorField mirrorPoints(2, vector::zero);

            label counter = -1;
        
            forAll (curPointEdges, edgeI)
            {
                label curEdge = curPointEdges[edgeI];

                if(edgeFaces[curEdge].size() == 1)
                {
                    label patchID = -1;
                    label edgeID = -1;
                    forAll(aMesh().boundary(), patchI)
                    {
                        const labelList& pEdges =
                            aMesh().boundary()[patchI];
                        label index = findIndex(pEdges, curEdge);
                        if (index != -1)
                        {
                            patchID = patchI;
                            edgeID = index;
                            break;
                        }
                    }

                    mirrorPoints[++counter] = 
                        patchMirrorPoints[patchID][edgeID];
                }
            }

            // Calculating LS plane fit
            const labelList& curPointFaces = pointFaces[curPoint];
        
            vectorField lsPoints
            (
                curPointFaces.size() + mirrorPoints.size(), 
                vector::zero
            );
                
            counter = -1;

            for (label i=0; i<curPointFaces.size(); i++)
            {
                label curFace = curPointFaces[i];

                lsPoints[++counter] = controlPoints()[curFace];
            }

            for (label i=0; i<mirrorPoints.size(); i++)
            {
                lsPoints[++counter] = mirrorPoints[i];
            }

            vectorField pointAndNormal = 
                lsPlanePointAndNormal
                (
                    lsPoints, 
                    points[curPoint], 
                    pointNormals[curPoint]
                );

            vector& P = pointAndNormal[0];
            vector& N = pointAndNormal[1];

            displacement[curPoint] = 
                pointsDisplacementDir()[curPoint]
               *((P - points[curPoint])&N)
               /(pointsDisplacementDir()[curPoint]&N);
        }
    }


    // Calculate displacement of axis point
    forAll (aMesh().boundary(), patchI)
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
                label axisPoint = wedgePatch.axisPoint();
                
                displacement[axisPoint] =
                    pointsDisplacementDir()[axisPoint]
                   *(
                        pointsDisplacementDir()[axisPoint]
                       &(
                            controlPoints()[pointFaces[axisPoint][0]]
                       //         const labelList peFaces = 
//             labelList::subList
//             (
//                 aMesh().edgeOwner(),
//                 aMesh().boundary()[patchI].faPatch::size(),
//                 aMesh().boundary()[patchI].start()
//             );
//         //Info<< "peFaces "<<endl <<peFaces <<endl ;
   - points[axisPoint]
                        )
                    );
            }
        }
    }


    // Calculate displacement of processor patch points
    forAll (aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type() 
         == processorFaPatch::typeName
        )
        {
            const processorFaPatch& procPatch =
                refCast<const processorFaPatch>(aMesh().boundary()[patchI]);

            const labelList& patchPointLabels =
                procPatch.pointLabels();

            FieldField<Field, vector> lsPoints(patchPointLabels.size());
            forAll(lsPoints, pointI)
            {
                lsPoints.set(pointI, new vectorField(0, vector::zero));
            }

            const labelList& nonGlobalPatchPoints = 
                procPatch.nonGlobalPatchPoints();

            forAll(nonGlobalPatchPoints, pointI)
            {
                label curPatchPoint =
                    nonGlobalPatchPoints[pointI];

                label curPoint = 
                    patchPointLabels[curPatchPoint];

                const labelList& curPointFaces = pointFaces[curPoint];

                lsPoints[curPatchPoint].setSize(curPointFaces.size());

                forAll(curPointFaces, faceI)
                {
                    label curFace = curPointFaces[faceI];

                    lsPoints[curPatchPoint][faceI] = controlPoints()[curFace];
                }

#               include "boundaryProcessorFaPatchPoints.H"
            }

            scalar lsPointsSize = 0;
            forAll(lsPoints, pointI)
            {
                lsPointsSize +=
                    2*lsPoints[pointI].size()*sizeof(vector);
            }
            
            // Parallel data exchange
            {
                OPstream toNeighbProc
                (
                    Pstream::blocking,
                    procPatch.neighbProcNo(), 
                    lsPointsSize
                );

                toNeighbProc << lsPoints;
            }

            FieldField<Field, vector> ngbLsPoints(patchPointLabels.size());

            {
                IPstream fromNeighbProc
                (
                    Pstream::blocking,
                    procPatch.neighbProcNo(),
                    lsPointsSize
                );

                fromNeighbProc >> ngbLsPoints;
            }

            forAll(nonGlobalPatchPoints, pointI)
            {
                label curPatchPoint = 
                    nonGlobalPatchPoints[pointI];

                label curPoint = 
                    patchPointLabels[curPatchPoint];

                label curNgbPoint = procPatch.neighbPoints()[curPatchPoint];

                vectorField allLsPoints
                (
                    lsPoints[curPatchPoint].size()
                  + ngbLsPoints[curNgbPoint].size(),
                    vector::zero
                );
                
                label counter = -1;
                forAll(lsPoints[curPatchPoint], pointI)
                {
                    allLsPoints[++counter] = lsPoints[curPatchPoint][pointI];
                }
                forAll(ngbLsPoints[curNgbPoint], pointI)
                {
                    allLsPoints[++counter] = ngbLsPoints[curNgbPoint][pointI];
                }

                vectorField pointAndNormal = 
                    lsPlanePointAndNormal
                    (
                        allLsPoints, 
                        points[curPoint], 
                        pointNormals[curPoint]
                    );

                vector& P = pointAndNormal[0];
                vector& N = pointAndNormal[1];
                    
                if (motionPointsMask()[curPoint] != 0)
                {
                    displacement[curPoint] = 
                        pointsDisplacementDir()[curPoint]
                       *((P - points[curPoint])&N)
                       /(pointsDisplacementDir()[curPoint]&N);
                }
            }
        }
    }


    // Calculate displacement of global processor patch points
    if (aMesh().globalData().nGlobalPoints() > 0)
    {
        const labelList& spLabels =
            aMesh().globalData().sharedPointLabels();

        const labelList& addr = aMesh().globalData().sharedPointAddr();

        for (label k=0; k<aMesh().globalData().nGlobalPoints(); k++)
        {
            List<List<vector> > procLsPoints(Pstream::nProcs());

            label curSharedPointIndex = findIndex(addr, k);

            if (curSharedPointIndex != -1)
            {
                label curPoint = spLabels[curSharedPointIndex];

                const labelList& curPointFaces = pointFaces[curPoint];

                procLsPoints[Pstream::myProcNo()] = 
                    List<vector>(curPointFaces.size());

                forAll (curPointFaces, faceI)
                {
                    label curFace = curPointFaces[faceI];

                    procLsPoints[Pstream::myProcNo()][faceI] = 
                        controlPoints()[curFace];
                }
            }

            Pstream::gatherList(procLsPoints);
            Pstream::scatterList(procLsPoints);
                
            if (curSharedPointIndex != -1)
            {
                label curPoint = spLabels[curSharedPointIndex];
                
                label nAllPoints = 0;
                forAll(procLsPoints, procI)
                {
                    nAllPoints += procLsPoints[procI].size();
                }

                vectorField allPoints(nAllPoints, vector::zero);

                label counter = 0;
                forAll(procLsPoints, procI)
                {
                    forAll(procLsPoints[procI], pointI)
                    {
                        allPoints[counter++] =
                            procLsPoints[procI][pointI];
                    }
                }

                vectorField pointAndNormal = 
                    lsPlanePointAndNormal
                    (
                        allPoints, 
                        points[curPoint], 
                        pointNormals[curPoint]
                    );

                const vector& P = pointAndNormal[0];
                const vector& N = pointAndNormal[1];

                displacement[curPoint] = 
                    pointsDisplacementDir()[curPoint]
                   *((P - points[curPoint])&N)
                   /(pointsDisplacementDir()[curPoint]&N);
            }
        }
    }

    return tdisplacement;
}
*/

// tmp<vectorField> freeSurface::pointDisplacement(const scalarField& deltaH)
// {
//     const pointField& points = aMesh().patch().localPoints();
//     const labelListList& pointFaces = aMesh().patch().pointFaces();
//     const faceList& facePoints = aMesh().patch().localFaces();

//     controlPoints() += facesDisplacementDir()*deltaH;

//     tmp<vectorField> tdisplacement
//     (
//         new vectorField
//         (
//             points.size(),
//             vector::zero
//         )
//     );
    
//     vectorField& displacement = tdisplacement();


//     // Calculate displacement of internal points
//     const vectorField& pointNormals = aMesh().pointAreaNormals();
//     const edgeList& edges = aMesh().patch().edges();
//     labelList internalPoints = aMesh().internalPoints();

//     forAll (internalPoints, pointI)
//     {
//         label curPoint = internalPoints[pointI];

//         const labelList& curPointFaces = pointFaces[curPoint];

//         vectorField lsPoints(curPointFaces.size(), vector::zero);

//         for (label i=0; i<curPointFaces.size(); i++)
//         {
//             label curFace = curPointFaces[i];
            
//             lsPoints[i] = controlPoints()[curFace];
//         }
       
//         vectorField pointAndNormal = 
//             lsPlanePointAndNormal
//             (
//                 lsPoints, 
//                 points[curPoint], 
//                 pointNormals[curPoint]
//             );

//         vector& P = pointAndNormal[0];
//         vector& N = pointAndNormal[1];

//         displacement[curPoint] = 
//             pointsDisplacementDir()[curPoint]
//            *((P - points[curPoint])&N)
//            /(pointsDisplacementDir()[curPoint]&N);
//     }


//     ///RG-TQ Nf
//     vectorField Nf = aMesh().faceAreaNormals().internalField();
//     dTheta_ = 0.0;

//     // Mirror control points
//     FieldField<Field, vector> patchMirrorPoints(aMesh().boundary().size());

//     Info<< "freeSurfacePoints "<<endl<<points<<endl ;
//     Info<< "controlPoints "<<endl<<controlPoints()<<endl ;
//     Info<< "facesDisplacementDir "<<endl<<facesDisplacementDir()<<endl ;

//     forAll(patchMirrorPoints, patchI)
//     {
//         patchMirrorPoints.set
//         (
//             patchI,
//             new vectorField
//             (
//                 aMesh().boundary()[patchI].faPatch::size(), 
//                 vector::zero
//             )
//         );

//         ///RG-TQ
//         //unallocLabelList patchEdgeFaces = 
//             //aMesh().boundary()[patchI].edgeFaces();

//         vectorField N = 
//             aMesh().boundary()[patchI].ngbPolyPatchFaceNormals();

//         Info<< "ngbPolyPatchFaceNormals "<<endl<<N<<endl ;

//         ///RG-TQ
//         //vector j = vector (0, 1, 0);
            
// //        scalar theta = contactangle().value()*(3.1415926/180.0);
//         scalar theta = 88*(3.1415926/180.0);
        
//         const labelList peFaces = 
//             labelList::subList
//             (
//                 aMesh().edgeOwner(),
//                 aMesh().boundary()[patchI].faPatch::size(),
//                 aMesh().boundary()[patchI].start()
//             );
//         //Info<< "peFaces "<<endl <<peFaces <<endl ;

//         const labelList& pEdges = aMesh().boundary()[patchI];
//        // Info<< "pEdges "<<endl <<pEdges <<endl ;

//         vectorField peCentres(pEdges.size(), vector::zero);

//         forAll(peCentres, edgeI)
//         {
//             peCentres[edgeI] = 
//                 edges[pEdges[edgeI]].centre(points);
//         }
//         //Info<< "peCentres "<<endl <<peCentres <<endl ;

//         vectorField delta =
// 	            vectorField(controlPoints(), peFaces)
// 	          - peCentres;

//         if(aMesh().boundary()[patchI].type() == emptyFaPatch::typeName)
//         {
// 	        patchMirrorPoints[patchI] = 
// 	            peCentres + ((I - 2*N*N)&delta);
//         }
//         else
//         {
// 	        vectorField nw(pEdges.size(), vector::zero);
// 	        vectorField nf = vectorField(Nf, peFaces);
// 	        vectorField ni(pEdges.size(), vector::zero);
// 	        vectorField ns(pEdges.size(), vector::zero);
// 	        vectorField t(pEdges.size(), vector::zero);
// 	        vectorField delta1 =
// 	            vectorField(aMesh().areaCentres().internalField(), peFaces)
// 	          - peCentres;

// 		scalarField theta1 = asin(mag(N&delta1)/mag(delta1));

// 	        Info<< "theta0 "<<theta <<endl ;
// 	        Info<< "theta1 "<<theta1 <<endl ;

// // Info<< "unit face normals "<<endl << nf <<endl ;

// 	        forAll(nw, edgeI)
// 	            nw[edgeI] = -N[edgeI];
// // Info<< "unit wall normals "<<endl << nw <<endl ;
	
// 	        forAll(ni, edgeI)
// 		    ni[edgeI] = nw[edgeI] ^ nf[edgeI];
// 	        ni /= mag(ni);
// // Info<< "unit vector along the contact line "<<endl << ni <<endl ;

// 	        forAll(ns, edgeI)
//             {
// 	            ns[edgeI] = (ni[edgeI] ^ nw[edgeI]);
// 	            ns[edgeI] = sign(ns[edgeI] & nf[edgeI])*ns[edgeI];
//             }
// // Info<< "unit vector normal to the contact line "<<endl << ns <<endl ;

// // RG-TQ Normal to the virtual symmetry plane:
// 	        forAll(t, edgeI)
// 	            t[edgeI] = sin(theta)*nw[edgeI] - cos(theta)*ns[edgeI];
// // Info<< "virtual symmetry plane normals "<<endl << t <<endl ;

// //	        Info<< "deltaH"<< endl <<deltaH <<endl ;

// 	        Info<< "Delta "<<endl<<delta<<endl ;
// 	        Info<< "Delta1 "<<endl<<delta1<<endl ;

//                 dTheta_ = max(dTheta_,max(mag(theta1-theta)));

// 	        patchMirrorPoints[patchI] = 
// 	            peCentres + ((I - 2*t*t)&delta);
//         }

//  Info<< "patchMirrorPoints["<<patchI<<"] "<<patchMirrorPoints[patchI]<<endl ;
//     }


//     // Calculate displacement of boundary points 
//     labelList boundaryPoints = aMesh().boundaryPoints();

//     const labelListList& edgeFaces = aMesh().patch().edgeFaces();
//     const labelListList& pointEdges = aMesh().patch().pointEdges();


//     forAll (boundaryPoints, pointI)
//     {
//         label curPoint = boundaryPoints[pointI];

//         if (motionPointsMask()[curPoint] == 1)
//         {
//             // Calculating mirror points
//             const labelList& curPointEdges = pointEdges[curPoint];

//             vectorField mirrorPoints(2, vector::zero);

//             label counter = -1;
        
//             forAll (curPointEdges, edgeI)
//             {
//                 label curEdge = curPointEdges[edgeI];

//                 if(edgeFaces[curEdge].size() == 1)
//                 {
//                     label patchID = -1;
//                     label edgeID = -1;
//                     forAll(aMesh().boundary(), patchI)
//                     {
//                         const labelList& pEdges =
//                             aMesh().boundary()[patchI];
//                         label index = findIndex(pEdges, curEdge);
//                         if (index != -1)
//                         {
//                             patchID = patchI;
//                             edgeID = index;
//                             break;
//                         }
//                     }

//                     mirrorPoints[++counter] = 
//                         patchMirrorPoints[patchID][edgeID];
//                 }
//             }

//             // Calculating LS plane fit
//             const labelList& curPointFaces = pointFaces[curPoint];
        
//             vectorField lsPoints
//             (
//                 curPointFaces.size() + mirrorPoints.size(), 
//                 vector::zero
//             );
                
//             counter = -1;

//             for (label i=0; i<curPointFaces.size(); i++)
//             {
//                 label curFace = curPointFaces[i];

//     Info<< "facePoints["<<curFace<<"] "<<facePoints[curFace] <<endl;

//                 lsPoints[++counter] = controlPoints()[curFace];
//             }

//             for (label i=0; i<mirrorPoints.size(); i++)
//             {
//                 lsPoints[++counter] = mirrorPoints[i];
//             }

//             vectorField pointAndNormal = 
//                 lsPlanePointAndNormal
//                 (
//                     lsPoints, 
//                     points[curPoint], 
//                     pointNormals[curPoint]
//                 );

//             vector& P = pointAndNormal[0];
//             vector& N = pointAndNormal[1];

//             displacement[curPoint] = 
//                 pointsDisplacementDir()[curPoint]
//                *((P - points[curPoint])&N)
//                /(pointsDisplacementDir()[curPoint]&N);

//     Info<< "pointAndNormal "<<endl<<pointAndNormal <<endl;

//     Info<< "lsPoints "<<endl<<lsPoints <<endl;
//     Info<< "points["<<curPoint<<"] "<<points[curPoint] <<endl;
//     Info<< "pointsNormals["<<curPoint<<"] "<<pointNormals[curPoint] <<endl;
//     Info<< "curPointFaces"<<curPointFaces <<endl;
//     Info<< "P "<<endl<<P <<endl;
//     Info<< "N "<<endl<<N <<endl;
//     Info<< "displacement["<<curPoint<<"] "<<displacement[curPoint] <<endl;

//         }

//     }
//     Info<< "displacement "<<endl<<displacement <<endl;

//     // Calculate displacement of axis point
//     forAll (aMesh().boundary(), patchI)
//     {
//         if
//         (
//             aMesh().boundary()[patchI].type() 
//          == wedgeFaPatch::typeName
//         )
//         {
//             const wedgeFaPatch& wedgePatch =
//                 refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

//             if(wedgePatch.axisPoint() > -1)
//             {
//                 label axisPoint = wedgePatch.axisPoint();
                
//                 displacement[axisPoint] =
//                     pointsDisplacementDir()[axisPoint]
//                    *(
//                         pointsDisplacementDir()[axisPoint]
//                        &(
//                             controlPoints()[pointFaces[axisPoint][0]]
//                           - points[axisPoint]
//                         )
//                     );
//             }
//         }
//     }


//     // Calculate displacement of processor patch points
//     forAll (aMesh().boundary(), patchI)
//     {
//         if
//         (
//             aMesh().boundary()[patchI].type() 
//          == processorFaPatch::typeName
//         )
//         {
//             const processorFaPatch& procPatch =
//                 refCast<const processorFaPatch>(aMesh().boundary()[patchI]);

//             const labelList& patchPointLabels =
//                 procPatch.pointLabels();

//             FieldField<Field, vector> lsPoints(patchPointLabels.size());
//             forAll(lsPoints, pointI)
//             {
//                 lsPoints.set(pointI, new vectorField(0, vector::zero));
//             }

//             const labelList& nonGlobalPatchPoints = 
//                 procPatch.nonGlobalPatchPoints();

//             forAll(nonGlobalPatchPoints, pointI)
//             {
//                 label curPatchPoint =
//                     nonGlobalPatchPoints[pointI];

//                 label curPoint = 
//                     patchPointLabels[curPatchPoint];

//                 const labelList& curPointFaces = pointFaces[curPoint];

//                 lsPoints[curPatchPoint].setSize(curPointFaces.size());

//                 forAll(curPointFaces, faceI)
//                 {
//                     label curFace = curPointFaces[faceI];

//                     lsPoints[curPatchPoint][faceI] = controlPoints()[curFace];
//                 }

// #               include "boundaryProcessorFaPatchPoints.H"
//             }

//             scalar lsPointsSize = 0;
//             forAll(lsPoints, pointI)
//             {
//                 lsPointsSize +=
//                     2*lsPoints[pointI].size()*sizeof(vector);
//             }
            
//             // Parallel data exchange
//             {
//                 OPstream toNeighbProc
//                 (
//                     Pstream::blocking,
//                     procPatch.neighbProcNo(), 
//                     lsPointsSize
//                 );

//                 toNeighbProc << lsPoints;
//             }

//             FieldField<Field, vector> ngbLsPoints(patchPointLabels.size());

//             {
//                 IPstream fromNeighbProc
//                 (
//                     Pstream::blocking,
//                     procPatch.neighbProcNo(),
//                     lsPointsSize
//                 );

//                 fromNeighbProc >> ngbLsPoints;
//             }

//             forAll(nonGlobalPatchPoints, pointI)
//             {
//                 label curPatchPoint = 
//                     nonGlobalPatchPoints[pointI];

//                 label curPoint = 
//                     patchPointLabels[curPatchPoint];

//                 label curNgbPoint = procPatch.neighbPoints()[curPatchPoint];

//                 vectorField allLsPoints
//                 (
//                     lsPoints[curPatchPoint].size()
//                   + ngbLsPoints[curNgbPoint].size(),
//                     vector::zero
//                 );
                
//                 label counter = -1;
//                 forAll(lsPoints[curPatchPoint], pointI)
//                 {
//                     allLsPoints[++counter] = lsPoints[curPatchPoint][pointI];
//                 }
//                 forAll(ngbLsPoints[curNgbPoint], pointI)
//                 {
//                     allLsPoints[++counter] = ngbLsPoints[curNgbPoint][pointI];
//                 }

//                 vectorField pointAndNormal = 
//                     lsPlanePointAndNormal
//                     (
//                         allLsPoints, 
//                         points[curPoint], 
//                         pointNormals[curPoint]
//                     );

//                 vector& P = pointAndNormal[0];
//                 vector& N = pointAndNormal[1];
                    
//                 if (motionPointsMask()[curPoint] != 0)
//                 {
//                     displacement[curPoint] = 
//                         pointsDisplacementDir()[curPoint]
//                        *((P - points[curPoint])&N)
//                        /(pointsDisplacementDir()[curPoint]&N);
//                 }
//             }
//         }
//     }


//     // Calculate displacement of global processor patch points
//     if (aMesh().globalData().nGlobalPoints() > 0)
//     {
//         const labelList& spLabels =
//             aMesh().globalData().sharedPointLabels();

//         const labelList& addr = aMesh().globalData().sharedPointAddr();

//         for (label k=0; k<aMesh().globalData().nGlobalPoints(); k++)
//         {
//             List<List<vector> > procLsPoints(Pstream::nProcs());

//             label curSharedPointIndex = findIndex(addr, k);

//             if (curSharedPointIndex != -1)
//             {
//                 label curPoint = spLabels[curSharedPointIndex];

//                 const labelList& curPointFaces = pointFaces[curPoint];

//                 procLsPoints[Pstream::myProcNo()] = 
//                     List<vector>(curPointFaces.size());

//                 forAll (curPointFaces, faceI)
//                 {
//                     label curFace = curPointFaces[faceI];

//                     procLsPoints[Pstream::myProcNo()][faceI] = 
//                         controlPoints()[curFace];
//                 }
//             }

//             Pstream::gatherList(procLsPoints);
//             Pstream::scatterList(procLsPoints);
                
//             if (curSharedPointIndex != -1)
//             {
//                 label curPoint = spLabels[curSharedPointIndex];
                
//                 label nAllPoints = 0;
//                 forAll(procLsPoints, procI)
//                 {
//                     nAllPoints += procLsPoints[procI].size();
//                 }

//                 vectorField allPoints(nAllPoints, vector::zero);

//                 label counter = 0;
//                 forAll(procLsPoints, procI)
//                 {
//                     forAll(procLsPoints[procI], pointI)
//                     {
//                         allPoints[counter++] =
//                             procLsPoints[procI][pointI];
//                     }
//                 }

//                 vectorField pointAndNormal = 
//                     lsPlanePointAndNormal
//                     (
//                         allPoints, 
//                         points[curPoint], 
//                         pointNormals[curPoint]
//                     );

//                 const vector& P = pointAndNormal[0];
//                 const vector& N = pointAndNormal[1];

//                 displacement[curPoint] = 
//                     pointsDisplacementDir()[curPoint]
//                    *((P - points[curPoint])&N)
//                    /(pointsDisplacementDir()[curPoint]&N);
//             }
//         }
//     }

//     return tdisplacement;
// }


tmp<vectorField> freeSurface::lsPlanePointAndNormal
(
    const vectorField& points,
    const vector& origin,
    const vector& axis
) const
{
    // LS in local CS
    vector dir = (points[0] - origin);
    dir -= axis*(axis&dir);
    dir /= mag(dir);
    coordinateSystem cs("cs", origin, axis, dir);

    vectorField localPoints = cs.localPosition(points);
    scalarField W = 1.0/(mag(points - origin) + SMALL);

    scalarRectangularMatrix M
    (
        points.size(),
        3,
        0.0
    );

    for (label i=0; i<localPoints.size(); i++)
    {
        M[i][0] = localPoints[i].x();
        M[i][1] = localPoints[i].y();
        M[i][2] = 1.0;
    }

    scalarSquareMatrix MtM(3, 0.0);
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

    scalarField MtR(3, 0);
    for (label i = 0; i < MtR.size(); i++)
    {
        for (label j = 0; j < M.n(); j++)
        {
            MtR[i] += M[j][i]*localPoints[j].z()*W[j];
        }
    }

    scalarSquareMatrix::LUsolve(MtM, MtR);

    vector n0 = vector(-MtR[0], -MtR[1], 1);
    n0 = cs.globalVector(n0);
    n0 /= mag(n0);

    vector p0 = vector(0, 0, MtR[2]);
    p0 = cs.globalPosition(p0);

    tmp<vectorField> pointAndNormal
    (
        new vectorField(2, vector::zero)
    );

    pointAndNormal()[0] = p0;
    pointAndNormal()[1] = n0;

    return pointAndNormal;
}

//Zeljko 0716 versions for parallelal gebraic mesh motion solver
void freeSurface::correctContactLinePointNormals()
{
    // Correct normals for contact line points 
    // according to specified contact angle

    vectorField& N =
        const_cast<vectorField&>
        (
            aMesh().pointAreaNormals()
        );

    if (contactAnglePtr_)
    {
//        Pout << "Correcting contact line normals" << endl;

        vectorField oldPoints(aMesh().nPoints(), vector::zero);

        const labelList& meshPoints = aMesh().patch().meshPoints();

        forAll(oldPoints, ptI)
        {
            oldPoints[ptI] = 
                mesh().oldPoints()[meshPoints[ptI]];
        }

#       include "createTangentField.H"

        forAll(aMesh().boundary(), patchI)
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
                    scalar rotAngle =
                        average
                        (
                            90
                          - contactAnglePtr_->boundaryField()[patchI]
                        );

                    rotAngle *= M_PI/180.0;

                    vectorField ngbN =
                        aMesh().boundary()[patchI].ngbPolyPatchPointNormals();

                    const labelList& patchPoints = 
                        aMesh().boundary()[patchI].pointLabels();

                    vectorField pN = vectorField(N, patchPoints);

                    vectorField rotationAxis = (ngbN^pN);
                    rotationAxis /= mag(rotationAxis) + SMALL;


                    // Calc rotation axis using edge vectors

                    const edgeList& edges = aMesh().edges();

                    const labelListList& pointEdges =
                        aMesh().boundary()[patchI].pointEdges();

                    forAll (pointEdges, pointI)
                    {
                        vector rotAx = vector::zero;

                        forAll(pointEdges[pointI], eI)
                        {
                            label curEdge = 
                                aMesh().boundary()[patchI].start()
                              + pointEdges[pointI][eI];

                            vector e = edges[curEdge].vec(oldPoints);

                            e *= (e&rotationAxis[pointI])
                               /mag(e&rotationAxis[pointI]);

                            e /= mag(e) + SMALL;

                            rotAx += e;
                        }

                        if (pointEdges[pointI].size() == 1)
                        {
#                           include "addNgbProcessorEdgeTangent.H"
                        }

                        rotationAxis[pointI] = rotAx/(mag(rotAx) + SMALL);
                    }

                    // Rodrigues' rotation formula
                    ngbN = ngbN*cos(rotAngle)
                      + rotationAxis*(rotationAxis & ngbN)*(1 - cos(rotAngle))
                      + (rotationAxis^ngbN)*sin(rotAngle);
                
                    forAll (patchPoints, pointI)
                    {
                        N[patchPoints[pointI]] -= 
                            ngbN[pointI]*(ngbN[pointI]&N[patchPoints[pointI]]);

                        N[patchPoints[pointI]] /= mag(N[patchPoints[pointI]]);
                    }
                }
            }            
        }
    }
}


/*
//z0312
void freeSurface::correctContactLinePointNormals()
{
    // Correct normals for contact line points 
    // according to specified contact angle

    vectorField& N =
        const_cast<vectorField&>
        (
            aMesh().pointAreaNormals()
        );

    if (contactAnglePtr_)
    {
        Info << "Correcting contact line normals" << endl;

        vectorField oldPoints(aMesh().nPoints(), vector::zero);

        const labelList& meshPoints = aMesh().patch().meshPoints();

        forAll(oldPoints, pI)
        {
            oldPoints[pI] = 
                mesh().oldPoints()[meshPoints[pI]];
        }

        forAll(aMesh().boundary(), patchI)
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
                    scalar rotAngle =
                        average
                        (
                            90
                          - contactAnglePtr_->boundaryField()[patchI]
                        );

                    rotAngle *= M_PI/180.0;
                    
                    vectorField ngbN =
                        aMesh().boundary()[patchI].ngbPolyPatchPointNormals();

                    const labelList& patchPoints = 
                        aMesh().boundary()[patchI].pointLabels();

                    vectorField pN = vectorField(N, patchPoints);

                    vectorField rotationAxis = (ngbN^pN);
                    rotationAxis /= mag(rotationAxis) + SMALL;


                    // Calc rotation axis using edge vectors

                    const edgeList& edges = aMesh().edges();

                    const labelListList& pointEdges =
                        aMesh().boundary()[patchI].pointEdges();

                    forAll (pointEdges, pointI)
                    {
                        vector rotAx = vector::zero;

                        forAll(pointEdges[pointI], eI)
                        {
                            label curEdge = 
                                aMesh().boundary()[patchI].start()
                              + pointEdges[pointI][eI];

//                             vector t = edges[curEdge].vec(aMesh().points());

                            vector e = edges[curEdge].vec(oldPoints);

                            e *= (e&rotationAxis[pointI])
                               /mag(e&rotationAxis[pointI]);

                            e /= mag(e) + SMALL;

                            rotAx += e;
                        }

                        rotationAxis[pointI] = rotAx/(mag(rotAx) + SMALL);
                    }

                    // Rodrigues' rotation formula
                    ngbN = ngbN*cos(rotAngle)
                      + rotationAxis*(rotationAxis & ngbN)*(1 - cos(rotAngle))
                      + (rotationAxis^ngbN)*sin(rotAngle);
                
                    forAll (patchPoints, pointI)
                    {
                        N[patchPoints[pointI]] -= 
                            ngbN[pointI]*(ngbN[pointI]&N[patchPoints[pointI]]);

                        N[patchPoints[pointI]] /= mag(N[patchPoints[pointI]]);
                    }
                }
            }            
        }
    }
}
*/

/*
void freeSurface::correctContactLinePointNormals()
{
    // Correct normals for contact line points 
    // according to specified contact angle

    vectorField& N =
        const_cast<vectorField&>
        (
            aMesh().pointAreaNormals()
        ); 

    if (contactAnglePtr_)
    {
        forAll(aMesh().boundary(), patchI)
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
                    scalar rotAngle = 
                        average
                        (
                            mag(contactAnglePtr_->boundaryField()[patchI]) 
                          - 90
                        )*M_PI/180.0;

		    scalar sA = sin(rotAngle);
		    scalar cA = cos(rotAngle);
                    
                    vector corrN;
                    vector rotationAxis;

                    vectorField ngbN =
                        aMesh().boundary()[patchI].ngbPolyPatchPointNormals();
                
                    labelList patchPoints = 
                        aMesh().boundary()[patchI].pointLabels();

                    forAll (patchPoints, pointI)
                    {
                        corrN = ngbN[pointI];

                        rotationAxis = - corrN^N[patchPoints[pointI]];
                        rotationAxis /= mag(rotationAxis) + SMALL;

                        // Rodrigues' rotation formula
                        corrN = corrN*cA
                              + rotationAxis*(rotationAxis&corrN)*(1-cA)
                              + (rotationAxis^corrN)*sA;

                        N[patchPoints[pointI]] -= 
                            corrN*(corrN&N[patchPoints[pointI]]);
                    }
               }
            }            
        }
    }

    N /= mag(N);
}
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

