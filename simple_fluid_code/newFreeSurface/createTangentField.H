areaVectorField tangent
(
    IOobject
    (
        "tangent",
        DB().timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    aMesh(),
    dimensioned<vector>("0", dimless, vector::zero)
);

if (Pstream::parRun())
{
    const labelListList& edgeFaces = aMesh().patch().edgeFaces();
    const labelListList& pointEdges = aMesh().patch().pointEdges();
    const labelListList& pointFaces = aMesh().patch().pointFaces();
    const edgeList& edges = aMesh().edges();

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

            forAll(patchPointLabels, pointI)
            {
                label curPoint = patchPointLabels[pointI];

                // Check if processor point is boundary point

                label patchID = -1;
                label edgeID = -1;

                const labelList& curPointEdges = pointEdges[curPoint];

                forAll (curPointEdges, edgeI)
                {
                    label curEdge = curPointEdges[edgeI];

                    if (edgeFaces[curEdge].size() == 1)
                    {
                        forAll(aMesh().boundary(), pI)
                        {
                            const labelList& curEdges = aMesh().boundary()[pI];

                            label index = findIndex(curEdges, curEdge);

                            if (index != -1)
                            {
                                if
                                (
                                    aMesh().boundary()[pI].type()
                                 != processorFaPatch::typeName
                                )
                                {
                                    patchID = pI;
                                    edgeID = index;
                                    break;
                                }
                            }
                        }
                    }
                }

                if (patchID != -1)
                {
                    label curEdge =
                        aMesh().boundary()[patchID].start() + edgeID;
                    
                    vector t = edges[curEdge].vec(oldPoints);
                    t /= mag(t) + SMALL;
                    
//                     Pout << aMesh().boundary()[patchID].name() << ", " 
//                         << t << endl;

                    const labelList& curPointFaces = pointFaces[curPoint];

                    forAll(curPointFaces, fI)
                    {
                        tangent.internalField()[curPointFaces[fI]] = t;
                            
                    }
                }
            }
        }
    }

    tangent.correctBoundaryConditions();
}
