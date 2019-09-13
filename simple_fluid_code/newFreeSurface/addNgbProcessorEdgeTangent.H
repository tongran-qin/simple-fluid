{
//     const labelList& patchPointLabels =
//         aMesh().boundary()[patchI].pointLabels();

    label curPoint = patchPoints[pointI];

    const labelListList& ptEdges = aMesh().patch().pointEdges();
    const labelList& curPointEdges = ptEdges[curPoint];

    label procPatchID = -1;
    label edgeID = -1;

    const labelListList& edgeFaces = aMesh().patch().edgeFaces();

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
                     == processorFaPatch::typeName
                    )
                    {
                        procPatchID = pI;
                        edgeID = index;
                        break;
                    }
                }
            }
        }
    }

    if (procPatchID != -1)
    {
        vector t = 
            tangent.boundaryField()[procPatchID]
           .patchNeighbourField()()[edgeID];

        t *= (t&rotationAxis[pointI])
            /(mag(t&rotationAxis[pointI]) + SMALL);

        t /= mag(t) + SMALL;
        
//         Pout << "--------->" << t << endl;

        rotAx += t;
    }
}
