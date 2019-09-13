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
    Foam::liquidFilmFvMotionSolver

Description
    Mesh motion solver for an fvMesh.

SourceFiles
    liquidFilmFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef liquidFilmFvMotionSolver_H
#define liquidFilmFvMotionSolver_H

#include "fvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class liquidFilmFvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class liquidFilmFvMotionSolver
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

    // Private Member Functions

        //- Disallow default bitwise copy construct
        liquidFilmFvMotionSolver
        (
            const liquidFilmFvMotionSolver&
        );

        //- Disallow default bitwise assignment
        void operator=(const liquidFilmFvMotionSolver&);


public:

    //- Runtime type information
    TypeName("liquidFilm");


    // Constructors

        //- Construct from polyMesh and data stream
        liquidFilmFvMotionSolver
        (
            const polyMesh&,
            Istream& msDataUnused
        );


    // Destructor

        virtual ~liquidFilmFvMotionSolver();


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
