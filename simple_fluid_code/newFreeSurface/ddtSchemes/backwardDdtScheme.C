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

#include "backwardDdtScheme.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvMatrices.H"
#include "slipFvPatchFields.H"
#include "symmetryFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
scalar backwardDdtScheme<vector>::deltaT_() const
{
    return mesh().time().deltaT().value();
}


template<>
scalar backwardDdtScheme<vector>::deltaT0_() const
{
    return mesh().time().deltaT0().value();
}


template<>
template<>
scalar backwardDdtScheme<vector>::deltaT0_
(
    const GeometricField<vector, fvPatchField, volMesh>& vf
) const
{
    if (vf.oldTime().timeIndex() == vf.oldTime().oldTime().timeIndex())
    {
        return GREAT;
    }
    else
    {
        return deltaT0_();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<GeometricField<vector, fvPatchField, volMesh> >
backwardDdtScheme<vector>::fvcDdt
(
    const dimensioned<vector>& dt
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+dt.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_();

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        tmp<GeometricField<vector, fvPatchField, volMesh> > tdtdt
        (
            new GeometricField<vector, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                dimensioned<vector>
                (
                    "0",
                    dt.dimensions()/dimTime,
                    pTraits<vector>::zero
                )
            )
        );

        tdtdt().internalField() = rDeltaT.value()*dt.value()*
        (
            coefft - (coefft0*mesh().V0() - coefft00*mesh().V00())/mesh().V()
        );

        return tdtdt;
    }
    else
    {
        return tmp<GeometricField<vector, fvPatchField, volMesh> >
        (
            new GeometricField<vector, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                dimensioned<vector>
                (
                    "0",
                    dt.dimensions()/dimTime,
                    pTraits<vector>::zero
                ),
                calculatedFvPatchField<vector>::typeName
            )
        );
    }
}


template<>
tmp<GeometricField<vector, fvPatchField, volMesh> >
backwardDdtScheme<vector>::fvcDdt
(
    const GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        return tmp<GeometricField<vector, fvPatchField, volMesh> >
        (
            new GeometricField<vector, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*vf.dimensions(),
                rDeltaT.value()*
                (
                    coefft*vf.internalField() -
                    (
                        coefft0*vf.oldTime().internalField()*mesh().V0()
                      - coefft00*vf.oldTime().oldTime().internalField()
                       *mesh().V00()
                    )/mesh().V()
                ),
                rDeltaT.value()*
                (
                    coefft*vf.boundaryField() -
                    (
                        coefft0*vf.oldTime().boundaryField()
                      - coefft00*vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<vector, fvPatchField, volMesh> >
        (
            new GeometricField<vector, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*
                (
                    coefft*vf
                  - coefft0*vf.oldTime()
                  + coefft00*vf.oldTime().oldTime()
                )
            )
        );
    }
}


template<>
tmp<GeometricField<vector, fvPatchField, volMesh> >
backwardDdtScheme<vector>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        return tmp<GeometricField<vector, fvPatchField, volMesh> >
        (
            new GeometricField<vector, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.value()*rho.value()*
                (
                    coefft*vf.internalField() -
                    (
                        coefft0*vf.oldTime().internalField()*mesh().V0()
                      - coefft00*vf.oldTime().oldTime().internalField()
                       *mesh().V00()
                    )/mesh().V()
                ),
                rDeltaT.value()*rho.value()*
                (
                    coefft*vf.boundaryField() -
                    (
                        coefft0*vf.oldTime().boundaryField()
                      - coefft00*vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<vector, fvPatchField, volMesh> >
        (
            new GeometricField<vector, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*rho*
                (
                    coefft*vf
                  - coefft0*vf.oldTime()
                 + coefft00*vf.oldTime().oldTime()
                )
            )
        );
    }
}

template<>
tmp<GeometricField<vector, fvPatchField, volMesh> >
backwardDdtScheme<vector>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        return tmp<GeometricField<vector, fvPatchField, volMesh> >
        (
            new GeometricField<vector, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.value()*
                (
                    coefft*rho.internalField()*vf.internalField() -
                    (
                        coefft0*rho.oldTime().internalField()
                       *vf.oldTime().internalField()*mesh().V0()
                      - coefft00*rho.oldTime().oldTime().internalField()
                       *vf.oldTime().oldTime().internalField()*mesh().V00()
                    )/mesh().V()
                ),
                rDeltaT.value()*
                (
                    coefft*rho.boundaryField()*vf.boundaryField() -
                    (
                        coefft0*rho.oldTime().boundaryField()
                       *vf.oldTime().boundaryField()
                      - coefft00*rho.oldTime().oldTime().boundaryField()
                       *vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<vector, fvPatchField, volMesh> >
        (
            new GeometricField<vector, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*
                (
                    coefft*rho*vf
                  - coefft0*rho.oldTime()*vf.oldTime()
                  + coefft00*rho.oldTime().oldTime()*vf.oldTime().oldTime()
                )
            )
        );
    }
}


template<>
tmp<fvMatrix<vector> >
backwardDdtScheme<vector>::fvmDdt
(
    GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<vector> > tfvm
    (
        new fvMatrix<vector>
        (
            vf,
            vf.dimensions()*dimVol/dimTime
        )
    );

    fvMatrix<vector>& fvm = tfvm();

    scalar rDeltaT = 1.0/deltaT_();

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    fvm.diag() = (coefft*rDeltaT)*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT*
        (
            coefft0*vf.oldTime().internalField()*mesh().V0()
          - coefft00*vf.oldTime().oldTime().internalField()
           *mesh().V00()
        );
    }
    else
    {
        fvm.source() = rDeltaT*mesh().V()*
        (
            coefft0*vf.oldTime().internalField()
          - coefft00*vf.oldTime().oldTime().internalField()
        );
    }

    return tfvm;
}


template<>
tmp<fvMatrix<vector> >
backwardDdtScheme<vector>::fvmDdt
(
    const dimensionedScalar& rho,
    GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<vector> > tfvm
    (
        new fvMatrix<vector>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<vector>& fvm = tfvm();

    scalar rDeltaT = 1.0/deltaT_();

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    fvm.diag() = (coefft*rDeltaT*rho.value())*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT*rho.value()*
        (
            coefft0*vf.oldTime().internalField()*mesh().V0()
          - coefft00*vf.oldTime().oldTime().internalField()
           *mesh().V00()
        );
    }
    else
    {
        fvm.source() = rDeltaT*mesh().V()*rho.value()*
        (
            coefft0*vf.oldTime().internalField()
          - coefft00*vf.oldTime().oldTime().internalField()
        );
    }

    return tfvm;
}


template<>
tmp<fvMatrix<vector> >
backwardDdtScheme<vector>::fvmDdt
(
    const volScalarField& rho,
    GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<vector> > tfvm
    (
        new fvMatrix<vector>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<vector>& fvm = tfvm();

    scalar rDeltaT = 1.0/deltaT_();

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    fvm.diag() = (coefft*rDeltaT)*rho.internalField()*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT*
        (
            coefft0*rho.oldTime().internalField()
           *vf.oldTime().internalField()*mesh().V0()
          - coefft00*rho.oldTime().oldTime().internalField()
           *vf.oldTime().oldTime().internalField()*mesh().V00()
        );
    }
    else
    {
        fvm.source() = rDeltaT*mesh().V()*
        (
            coefft0*rho.oldTime().internalField()
           *vf.oldTime().internalField()
          - coefft00*rho.oldTime().oldTime().internalField()
           *vf.oldTime().oldTime().internalField()
        );
    }

    return tfvm;
}


template<>
tmp<surfaceScalarField>
backwardDdtScheme<vector>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const GeometricField<vector, fvPatchField, volMesh>& U,
    const surfaceScalarField& phi
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddtPhiCorr(" + rA.name() + ',' + U.name() + ',' + phi.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(U);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            ddtIOobject,
            rDeltaT*fvcDdtPhiCoeff(U.oldTime(), phi.oldTime())
           *(
                fvc::interpolate(rA)
               *(
                   coefft0*phi.oldTime()
                 - coefft00*phi.oldTime().oldTime()
                )
              - (
                    fvc::interpolate
                    (
                        rA*
                        (
                            coefft0*U.oldTime()
                          - coefft00*U.oldTime().oldTime()
                        )
                    ) & mesh().Sf()
                )
            )
        )
    );
}


template<>
tmp<surfaceScalarField>
backwardDdtScheme<vector>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const volScalarField& rho,
    const GeometricField<vector, fvPatchField, volMesh>& U,
    const surfaceScalarField& phi
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddtPhiCorr("
      + rA.name() + ','
      + rho.name() + ','
      + U.name() + ','
      + phi.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(U);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;


    if (mesh().moving())
    {
        volScalarField V0oV
        (
            IOobject
            (
                "V0oV",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimless,
            zeroGradientFvPatchScalarField::typeName
        );

        V0oV.internalField() = mesh().V0()/mesh().V();
        V0oV.correctBoundaryConditions();

        volScalarField V00oV
        (
            IOobject
            (
                "V00oV",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimless,
            zeroGradientFvPatchScalarField::typeName
        );

        V00oV.internalField() = mesh().V00()/mesh().V();
        V00oV.correctBoundaryConditions();

        surfaceScalarField ddtPhiCoeff
        (
            IOobject
            (
                "ddtPhiCoeff",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensioned<scalar>("1", dimless, 1.0)
        );


//         forAll (U.oldTime().boundaryField(), patchI)
//         {
//             if (!U.boundaryField()[patchI].coupled())
//             {
//                 ddtPhiCoeff.boundaryField()[patchI] = 0.0;
//             }
//         }

        forAll (U.boundaryField(), patchI)
        {
            if 
            (
                U.boundaryField()[patchI].fixesValue()
            )
            {
                ddtPhiCoeff.boundaryField()[patchI] = 0.0;
            }
            else if 
            (
                U.boundaryField()[patchI].type()
             == slipFvPatchVectorField::typeName
            )
            {
                ddtPhiCoeff.boundaryField()[patchI] = 0.0;
            }
            else if 
            (
                U.boundaryField()[patchI].type()
             == symmetryFvPatchVectorField::typeName
            )
            {
                ddtPhiCoeff.boundaryField()[patchI] = 0.0;
            }
        }

        if
        (
            U.dimensions() == dimVelocity
         && phi.dimensions() == dimVelocity*dimArea
        )
        {
            surfaceVectorField dU0 = fvc::interpolate(U.oldTime());
            forAll(dU0.boundaryField(), patchI)
            {
                if (!U.boundaryField()[patchI].coupled())
                {
                    dU0.boundaryField()[patchI] = 
                        U.oldTime().boundaryField()[patchI]
                       .patchInternalField();
                }
            }

            surfaceVectorField dU00 = fvc::interpolate(U.oldTime().oldTime());
            forAll(dU00.boundaryField(), patchI)
            {
                if (!U.boundaryField()[patchI].coupled())
                {
                    dU00.boundaryField()[patchI] = 
                        U.oldTime().oldTime().boundaryField()[patchI]
                       .patchInternalField();
                }
            }

            if(mesh().objectRegistry::foundObject<surfaceVectorField>("Sf"))
            {
                Info << "ZT, backwardDdtPhiCorr" << endl;

                const surfaceVectorField& Sf = 
                    mesh().objectRegistry::lookupObject<surfaceVectorField>
                    (
                        "Sf"
                    );

                dU0 = (phi.oldTime() - (Sf.oldTime()&dU0))
                    *Sf.oldTime()/magSqr(Sf.oldTime());

                dU00 = Sf.oldTime().oldTime()
                   *(phi.oldTime().oldTime() - (Sf.oldTime().oldTime()&dU00))
                   /magSqr(Sf.oldTime().oldTime());

//                 surfaceVectorField U0 = fvc::interpolate(U.oldTime());
//                 U0 -= (Sf.oldTime()&U0)*Sf.oldTime()/magSqr(Sf.oldTime());
//                 U0 += phi.oldTime()*Sf.oldTime()/magSqr(Sf.oldTime());

//                 surfaceVectorField U00 = 
//                     fvc::interpolate(U.oldTime().oldTime());
//                 U00 -= (Sf.oldTime().oldTime()&U00)*Sf.oldTime().oldTime()
//                     /magSqr(Sf.oldTime().oldTime());
//                 U00 += phi.oldTime().oldTime()*Sf.oldTime().oldTime()
//                     /magSqr(Sf.oldTime().oldTime());


//                 return tmp<surfaceScalarField>
//                 (
//                     new surfaceScalarField
//                     (
//                         ddtIOobject,
//                         rDeltaT*ddtPhiCoeff
//                        *(
//                             coefft0*fvc::interpolate(rho.oldTime()*V0oV)
//                            *(mesh().Sf()&U0)
//                           - coefft00
//                            *fvc::interpolate(rho.oldTime().oldTime()*V00oV)
//                            *(mesh().Sf()&U00)
//                           - (
//                                 fvc::interpolate
//                                 (
//                                     coefft0*rho.oldTime()*U.oldTime()*V0oV
//                                   - coefft00*rho.oldTime().oldTime()
//                                    *U.oldTime().oldTime()*V00oV
//                                 ) 
//                               & mesh().Sf()
//                             )
//                         )
//                        /fvc::interpolate(1.0/rA) 
//                     )
//                 );
            }
            else
            {
                Info << "ZT, backwardDdtPhiCorr,2" << endl;

                dU0 = (phi.oldTime() - (mesh().Sf()&dU0))
                    *mesh().Sf()/magSqr(mesh().Sf());

                dU00 = (phi.oldTime().oldTime() - (mesh().Sf()&dU00))
                   *mesh().Sf()/magSqr(mesh().Sf());
            }

            return tmp<surfaceScalarField>
            (
                new surfaceScalarField
                (
                    ddtIOobject,
                    rDeltaT*ddtPhiCoeff
                   *(
                        coefft0*fvc::interpolate(rho.oldTime()*V0oV)
                       *(mesh().Sf()&dU0)
                      - coefft00
                       *fvc::interpolate(rho.oldTime().oldTime()*V00oV)
                       *(mesh().Sf()&dU00)
                    )
                   /fvc::interpolate(1.0/rA) 
                )
            );

//             return tmp<surfaceScalarField>
//             (
//                 new surfaceScalarField
//                 (
//                     ddtIOobject,
//                     rDeltaT*ddtPhiCoeff
//                    *(
//                         coefft0*fvc::interpolate(rA*rho.oldTime()*V0oV)
//                        *phi.oldTime()
//                       - coefft00
//                        *fvc::interpolate(rA*rho.oldTime().oldTime()*V00oV)
//                        *phi.oldTime().oldTime()
//                       - (
//                             fvc::interpolate
//                             (
//                                 coefft0*rA*rho.oldTime()*U.oldTime()*V0oV
//                               - coefft00*rA*rho.oldTime().oldTime()
//                                *U.oldTime().oldTime()*V00oV
//                             ) & mesh().Sf()
//                         )
//                     )
//                 )
//             );
        }
        else if 
        (
            U.dimensions() == dimVelocity
         && phi.dimensions() == rho.dimensions()*dimVelocity*dimArea
        )
        {
            if(mesh().objectRegistry::foundObject<surfaceVectorField>("Sf"))
            {
                const surfaceVectorField& Sf = 
                    mesh().objectRegistry::lookupObject<surfaceVectorField>
                    (
                        "Sf"
                    );

                surfaceVectorField U0 = fvc::interpolate(U.oldTime());
                U0 -= (Sf.oldTime()&U0)*Sf.oldTime()/magSqr(Sf.oldTime());
                U0 += phi.oldTime()*Sf.oldTime()/magSqr(Sf.oldTime());

                surfaceVectorField U00 = 
                    fvc::interpolate(U.oldTime().oldTime());
                U00 -= (Sf.oldTime().oldTime()&U00)*Sf.oldTime().oldTime()
                    /magSqr(Sf.oldTime().oldTime());
                U00 += phi.oldTime().oldTime()*Sf.oldTime().oldTime()
                    /magSqr(Sf.oldTime().oldTime());


                return tmp<surfaceScalarField>
                (
                    new surfaceScalarField
                    (
                        ddtIOobject,
                        rDeltaT*ddtPhiCoeff
                       *(
                            coefft0*fvc::interpolate(rA*rho*V0oV)
                           *(mesh().Sf()&U0)
                          - coefft00*fvc::interpolate(rA*rho*V00oV)
                           *(mesh().Sf()&U00)
                          - (
                                fvc::interpolate
                                (
                                    rho*rA*
                                    (
                                        coefft0*rho.oldTime()*U.oldTime()*V0oV
                                      - coefft00*rho.oldTime().oldTime()
                                       *U.oldTime().oldTime()*V00oV
                                    )
                                ) & mesh().Sf()
                            )
                        )
                    )
                );
            }

            return tmp<surfaceScalarField>
            (
                new surfaceScalarField
                (
                    ddtIOobject,
                    rDeltaT*ddtPhiCoeff
                   *(
                        coefft0*fvc::interpolate(rho*rA*V0oV)*phi.oldTime()
                      - coefft00*fvc::interpolate(rho*rA*V00oV)
                       *phi.oldTime().oldTime()
                      - (
                            fvc::interpolate
                            (
                                rho*rA
                               *(
                                    coefft0*rho.oldTime()*U.oldTime()*V0oV
                                  - coefft00*rho.oldTime().oldTime()
                                   *U.oldTime().oldTime()*V00oV
                                )
                            ) & mesh().Sf()
                        )
                    )
                )
            );
        }
        else if 
        (
            U.dimensions() == rho.dimensions()*dimVelocity
         && phi.dimensions() == rho.dimensions()*dimVelocity*dimArea
        )
        {
            if(mesh().objectRegistry::foundObject<surfaceVectorField>("Sf"))
            {
                const surfaceVectorField& Sf = 
                    mesh().objectRegistry::lookupObject<surfaceVectorField>
                    (
                        "Sf"
                    );

                surfaceVectorField U0 = fvc::interpolate(U.oldTime());
                U0 -= (Sf.oldTime()&U0)*Sf.oldTime()/magSqr(Sf.oldTime());
                U0 += phi.oldTime()*Sf.oldTime()/magSqr(Sf.oldTime());

                surfaceVectorField U00 = 
                    fvc::interpolate(U.oldTime().oldTime());
                U00 -= (Sf.oldTime().oldTime()&U00)*Sf.oldTime().oldTime()
                    /magSqr(Sf.oldTime().oldTime());
                U00 += phi.oldTime().oldTime()*Sf.oldTime().oldTime()
                    /magSqr(Sf.oldTime().oldTime());

                
                return tmp<surfaceScalarField>
                (
                    new surfaceScalarField
                    (
                        ddtIOobject,
                        rDeltaT*ddtPhiCoeff
                       *(
                            coefft0*fvc::interpolate(rA*V0oV)*(mesh().Sf()&U0)
                          - coefft00*fvc::interpolate(rA*V00oV)
                           *(mesh().Sf()&U00)
                          - (
                                fvc::interpolate
                                (
                                    coefft0*rA*U.oldTime()*V0oV
                                  - coefft00*rA*U.oldTime().oldTime()*V00oV
                                ) & mesh().Sf()
                            )
                        )
                    )
                );
            }

            return tmp<surfaceScalarField>
            (
                new surfaceScalarField
                (
                    ddtIOobject,
                    rDeltaT*ddtPhiCoeff
                   *(
                        coefft0*fvc::interpolate(rA*V0oV)*phi.oldTime()
                      - coefft00*fvc::interpolate(rA*V00oV)
                       *phi.oldTime().oldTime()
                      - (
                            fvc::interpolate
                            (
                                coefft0*rA*U.oldTime()*V0oV
                              - coefft00*rA*U.oldTime().oldTime()*V00oV
                            ) & mesh().Sf()
                        )
                    )
                )
            );
        }
        else
        {
            FatalErrorIn
            (
                "backwardDdtScheme<vector>::fvcDdtPhiCorr"
            )   << "dimensions of phi are not correct"
                << abort(FatalError);

            return surfaceScalarField::null();
        }
    }
    else
    {
        if
        (
            U.dimensions() == dimVelocity
         && phi.dimensions() == dimVelocity*dimArea
        )
        {
            return tmp<surfaceScalarField>
            (
                new surfaceScalarField
                (
                    ddtIOobject,
                    rDeltaT*fvcDdtPhiCoeff(U.oldTime(), phi.oldTime())
                   *(
                        coefft0*fvc::interpolate(rA*rho.oldTime())
                       *phi.oldTime()
                      - coefft00*fvc::interpolate(rA*rho.oldTime().oldTime())
                       *phi.oldTime().oldTime()
                      - (
                            fvc::interpolate
                            (
                                rA*
                                (
                                    coefft0*rho.oldTime()*U.oldTime()
                                  - coefft00*rho.oldTime().oldTime()
                                   *U.oldTime().oldTime()
                                )
                            ) & mesh().Sf()
                        )
                    )
                )
            );
        }
        else if 
        (
            U.dimensions() == dimVelocity
         && phi.dimensions() == rho.dimensions()*dimVelocity*dimArea
        )
        {
            return tmp<surfaceScalarField>
            (
                new surfaceScalarField
                (
                    ddtIOobject,
                    rDeltaT
                   *fvcDdtPhiCoeff
                    (
                        rho.oldTime(),
                        rho.oldTime()*U.oldTime(),
                        phi.oldTime()
                    )
                   *(
                        fvc::interpolate(rA*rho)
                       *(
                           coefft0*phi.oldTime()
                         - coefft00*phi.oldTime().oldTime()
                        )
                      - (
                            fvc::interpolate
                            (
                                rA*rho*
                                (
                                    coefft0*rho.oldTime()*U.oldTime()
                                  - coefft00*rho.oldTime().oldTime()
                                   *U.oldTime().oldTime()
                                )
                            ) & mesh().Sf()
                        )
                    )
                )
            );
        }
        else if 
        (
            U.dimensions() == rho.dimensions()*dimVelocity
         && phi.dimensions() == rho.dimensions()*dimVelocity*dimArea
        )
        {
            return tmp<surfaceScalarField>
            (
                new surfaceScalarField
                (
                    ddtIOobject,
                    rDeltaT
                   *fvcDdtPhiCoeff(rho.oldTime(), U.oldTime(), phi.oldTime())
                   *(
                        fvc::interpolate(rA)
                       *(
                           coefft0*phi.oldTime()
                         - coefft00*phi.oldTime().oldTime()
                        )
                      - (
                            fvc::interpolate
                            (
                                rA*
                                (
                                    coefft0*U.oldTime()
                                  - coefft00*U.oldTime().oldTime()
                                )
                            ) & mesh().Sf()
                        )
                    )
                )
            );
        }
        else
        {
            FatalErrorIn
            (
                "backwardDdtScheme<vector>::fvcDdtPhiCorr"
            )   << "dimensions of phi are not correct"
                << abort(FatalError);

            return surfaceScalarField::null();
        }
    }

//     if
//     (
//         U.dimensions() == dimVelocity
//      && phiAbs.dimensions() == dimVelocity*dimArea
//     )
//     {
//         return tmp<surfaceScalarField>
//         (
//             new surfaceScalarField
//             (
//                 ddtIOobject,
//                 rDeltaT*fvcDdtPhiCoeff(U.oldTime(), phiAbs.oldTime())
//                *(
//                     coefft0*fvc::interpolate(rA*rho.oldTime())
//                    *phiAbs.oldTime()
//                   - coefft00*fvc::interpolate(rA*rho.oldTime().oldTime())
//                    *phiAbs.oldTime().oldTime()
//                   - (
//                         fvc::interpolate
//                         (
//                             rA*
//                             (
//                                 coefft0*rho.oldTime()*U.oldTime()
//                               - coefft00*rho.oldTime().oldTime()
//                                *U.oldTime().oldTime()
//                             )
//                         ) & mesh().Sf()
//                     )
//                 )
//             )
//         );
//     }
//     else if
//     (
//         U.dimensions() == dimVelocity
//      && phiAbs.dimensions() == rho.dimensions()*dimVelocity*dimArea
//     )
//     {
//         return tmp<surfaceScalarField>
//         (
//             new surfaceScalarField
//             (
//                 ddtIOobject,
//                 rDeltaT
//                *fvcDdtPhiCoeff
//                 (
//                     U.oldTime(),
//                     phiAbs.oldTime()/fvc::interpolate(rho.oldTime())
//                 )
//                *(
//                     fvc::interpolate(rA*rho.oldTime())
//                    *(
//                        coefft0*phiAbs.oldTime()
//                       /fvc::interpolate(rho.oldTime())
//                      - coefft00*phiAbs.oldTime().oldTime()
//                       /fvc::interpolate(rho.oldTime().oldTime())
//                     )
//                   - (
//                         fvc::interpolate
//                         (
//                             rA*rho.oldTime()*
//                             (
//                                 coefft0*U.oldTime()
//                               - coefft00*U.oldTime().oldTime()
//                             )
//                         ) & mesh().Sf()
//                     )
//                 )
//             )
//         );
//     }
//     else if
//     (
//         U.dimensions() == rho.dimensions()*dimVelocity
//      && phiAbs.dimensions() == rho.dimensions()*dimVelocity*dimArea
//     )
//     {
//         return tmp<surfaceScalarField>
//         (
//             new surfaceScalarField
//             (
//                 ddtIOobject,
//                 rDeltaT
//                *fvcDdtPhiCoeff(rho.oldTime(), U.oldTime(), phiAbs.oldTime())
//                *(
//                     fvc::interpolate(rA)
//                    *(
//                        coefft0*phiAbs.oldTime()
//                      - coefft00*phiAbs.oldTime().oldTime()
//                     )
//                   - (
//                         fvc::interpolate
//                         (
//                             rA*
//                             (
//                                 coefft0*U.oldTime()
//                               - coefft00*U.oldTime().oldTime()
//                             )
//                         ) & mesh().Sf()
//                     )
//                 )
//             )
//         );
//     }
//     else
//     {
//         FatalErrorIn
//         (
//             "backwardDdtScheme<vector>::fvcDdtPhiCorr"
//         )   << "dimensions of phiAbs are not correct"
//             << abort(FatalError);

//         return surfaceScalarField::null();
//     }
}


template<>
tmp<surfaceScalarField> backwardDdtScheme<vector>::meshPhi
(
    const GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    // Coefficient for t-3/2 (between times 0 and 00)
    scalar coefft0_00 = deltaT/(deltaT + deltaT0);

    // Coefficient for t-1/2 (between times n and 0)
    scalar coefftn_0 = 1 + coefft0_00;

    return coefftn_0*mesh().phi() - coefft0_00*mesh().phi().oldTime();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
