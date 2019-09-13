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

Class
    Foam::liquidFilmParallelFvMotionSolver

Description
    Mesh motion solver for an fvMesh.

SourceFiles
    liquidFilmParallelFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef liquidFilmParallelFvMotionSolver_H
#define liquidFilmParallelFvMotionSolver_H

#include "fvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class liquidFilmParallelFvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class liquidFilmParallelFvMotionSolver
:
    public fvMotionSolver
{
    // Private data

        //- Point motion field
        mutable pointVectorField pointMotionU_;

        //- Point to point addressing
        labelList pointPoint_;

        //- Top plane y-coordinate
        scalar max_;

        //- Bottom plane y-coordinate
        scalar min_;

        //- Motion direction
        vector motionDirection_;

        //- Free-surface patch ID
        label freeSurfacePatchID_;

        //- Free-surface face zone ID
        label freeSurfaceZoneID_;

        //- map from the proc faceZone point order to the master proc 
        //  point order
        labelList procToGlobalFZmap_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        liquidFilmParallelFvMotionSolver
        (
            const liquidFilmParallelFvMotionSolver&
        );

        //- Disallow default bitwise assignment
        void operator=(const liquidFilmParallelFvMotionSolver&);


public:

    //- Runtime type information
    TypeName("liquidFilmParallel");


    // Constructors

        //- Construct from polyMesh and data stream
        liquidFilmParallelFvMotionSolver
        (
            const polyMesh&,
            Istream& msDataUnused
        );


    // Destructor

        virtual ~liquidFilmParallelFvMotionSolver();


    // Member Functions

        //- Return reference to the point motion velocity field
        pointVectorField& pointMotionU()
        {
            return pointMotionU_;
        }

        //- Return const reference to the point motion velocity field
        const pointVectorField& pointMotionU() const
        {
            return pointMotionU_;
        }

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
