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

Description
    Simple central-difference snGrad scheme with non-orthogonal correction.

\*---------------------------------------------------------------------------*/

#include "skewCorrectedSnGrad.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "linear.H"
#include "reverseLinear.H"
#include "fvcGrad.H"
#include "gaussGrad.H"

#include "skewCorrectionVectors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
skewCorrectedSnGrad<Type>::~skewCorrectedSnGrad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
skewCorrectedSnGrad<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    // construct GeometricField<Type, fvsPatchField, surfaceMesh>
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tssf
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "snGradCorr("+vf.name()+')',
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            vf.dimensions()*mesh.deltaCoeffs().dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& ssf = tssf();

    typedef typename
        outerProduct<vector, typename pTraits<Type>::cmptType>::type 
        CmptGradType;

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        GeometricField<CmptGradType, fvPatchField, volMesh> cmptGrad =
            gradScheme<typename pTraits<Type>::cmptType>::New
            (
                mesh,
                mesh.gradScheme(ssf.name())
            )()
            .grad(vf.component(cmpt));

        // Non-orthogonal correction
        ssf.replace
        (
            cmpt,
            mesh.correctionVectors()
          & reverseLinear
            <
                typename 
                outerProduct<vector, typename pTraits<Type>::cmptType>::type
            >(mesh).interpolate(cmptGrad)
        );

        // Skewness correction
        if (skewCorrectionVectors::New(mesh).skew())
        {
            const skewCorrectionVectors& scv = 
                skewCorrectionVectors::New(mesh);
            
            Field<CmptGradType> DCmptGrad
            (
                ssf.internalField().size(), 
                pTraits<CmptGradType>::zero
            );

            const unallocLabelList& owner = mesh.owner();
            const unallocLabelList& neighbour = mesh.neighbour();

            forAll(DCmptGrad, faceI)
            {
                DCmptGrad[faceI] = 
                    cmptGrad[neighbour[faceI]]
                  - cmptGrad[owner[faceI]];
            }

            ssf.internalField().replace
            (
                cmpt,
                ssf.internalField().component(cmpt)
              + (scv().internalField()&DCmptGrad)
            );
        }
    }

    return tssf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
