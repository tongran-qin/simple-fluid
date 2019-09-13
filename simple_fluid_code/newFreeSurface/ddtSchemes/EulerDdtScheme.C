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

#include "EulerDdtScheme.H"
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
tmp<GeometricField<vector, fvPatchField, volMesh> >
EulerDdtScheme<vector>::fvcDdt
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

        tdtdt().internalField() =
            rDeltaT.value()*dt.value()*(1.0 - mesh().V0()/mesh().V());

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
EulerDdtScheme<vector>::fvcDdt
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
                    vf.internalField()
                  - vf.oldTime().internalField()*mesh().V0()/mesh().V()
                ),
                rDeltaT.value()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
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
                rDeltaT*(vf - vf.oldTime())
            )
        );
    }
}


template<>
tmp<GeometricField<vector, fvPatchField, volMesh> >
EulerDdtScheme<vector>::fvcDdt
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
                    vf.internalField()
                  - vf.oldTime().internalField()*mesh().V0()/mesh().V()
                ),
                rDeltaT.value()*rho.value()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
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
                rDeltaT*rho*(vf - vf.oldTime())
            )
        );
    }
}


template<>
tmp<GeometricField<vector, fvPatchField, volMesh> >
EulerDdtScheme<vector>::fvcDdt
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
                    rho.internalField()*vf.internalField()
                  - rho.oldTime().internalField()
                   *vf.oldTime().internalField()*mesh().V0()/mesh().V()
                ),
                rDeltaT.value()*
                (
                    rho.boundaryField()*vf.boundaryField()
                  - rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
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
                rDeltaT*(rho*vf - rho.oldTime()*vf.oldTime())
            )
        );
    }
}


template<>
tmp<fvMatrix<vector> >
EulerDdtScheme<vector>::fvmDdt
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

    scalar rDeltaT = 1.0/mesh().time().deltaT().value();

    fvm.diag() = rDeltaT*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT*vf.oldTime().internalField()*mesh().V0();
    }
    else
    {
        fvm.source() = rDeltaT*vf.oldTime().internalField()*mesh().V();
    }

    return tfvm;
}


template<>
tmp<fvMatrix<vector> >
EulerDdtScheme<vector>::fvmDdt
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

    scalar rDeltaT = 1.0/mesh().time().deltaT().value();

    fvm.diag() = rDeltaT*rho.value()*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT
            *rho.value()*vf.oldTime().internalField()*mesh().V0();
    }
    else
    {
        fvm.source() = rDeltaT
            *rho.value()*vf.oldTime().internalField()*mesh().V();
    }

    return tfvm;
}


template<>
tmp<fvMatrix<vector> >
EulerDdtScheme<vector>::fvmDdt
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

    scalar rDeltaT = 1.0/mesh().time().deltaT().value();

    fvm.diag() = rDeltaT*rho.internalField()*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT
            *rho.oldTime().internalField()
            *vf.oldTime().internalField()*mesh().V0();
    }
    else
    {
        fvm.source() = rDeltaT
            *rho.oldTime().internalField()
            *vf.oldTime().internalField()*mesh().V();
    }

    return tfvm;
}


template<>
tmp<surfaceScalarField> EulerDdtScheme<vector>::fvcDdtPhiCorr
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

        surfaceScalarField ddtPhiCoeff
        (
            IOobject
            (
                "ddtCouplingCoeff",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensioned<scalar>("1", dimless, 1.0)
        );

        forAll (U.oldTime().boundaryField(), patchI)
        {
            if (U.boundaryField()[patchI].fixesValue())
            {
                ddtPhiCoeff.boundaryField()[patchI] = 0.0;
            }
//             else if 
//             (
//                 U.boundaryField()[patchI].type()
//              == slipFvPatchVectorField::typeName
//             )
//             {
//                 ddtPhiCoeff.boundaryField()[patchI] = 0.0;
//             }
        }

        if(mesh().objectRegistry::foundObject<surfaceVectorField>("Sf"))
        {
            const surfaceVectorField& Sf = 
                mesh().objectRegistry::lookupObject<surfaceVectorField>("Sf");

            surfaceVectorField U0 = fvc::interpolate(U.oldTime());
            U0 -= (Sf.oldTime()&U0)*Sf.oldTime()/magSqr(Sf.oldTime());
            U0 += phi.oldTime()*Sf.oldTime()/magSqr(Sf.oldTime());

            return tmp<surfaceScalarField>
            (
                new surfaceScalarField
                (
                    ddtIOobject,
                    rDeltaT*ddtPhiCoeff
                   *(
                        fvc::interpolate(rA*V0oV)*(mesh().Sf()&U0)
                      - (
                            fvc::interpolate(rA*U.oldTime()*V0oV) 
                          & mesh().Sf()
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
                    phi.oldTime()*fvc::interpolate(rA*V0oV)
                  - (
                        fvc::interpolate(rA*U.oldTime()*V0oV) 
                      & mesh().Sf()
                    )
                )
            )
        );
    }
    else
    {
        Info << "Zt, ddtPhiCorr" << endl;

        surfaceScalarField ddtPhiCoeff
        (
            IOobject
            (
                "ddtCouplingCoeff",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensioned<scalar>("1", dimless, 1.0)
        );

        forAll (U.oldTime().boundaryField(), patchI)
        {
            if (U.boundaryField()[patchI].fixesValue())
            {
                ddtPhiCoeff.boundaryField()[patchI] = 0.0;
            }
        }

        tmp<surfaceScalarField> phiCorr =
            phi.oldTime() - (fvc::interpolate(U.oldTime()) & mesh().Sf());

        return tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                ddtIOobject,
                ddtPhiCoeff
               *fvc::interpolate(rDeltaT*rA)*phiCorr
            )
        );
    }
}


template<>
tmp<surfaceScalarField> EulerDdtScheme<vector>::fvcDdtPhiCorr
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
      + rA.name() + ',' + rho.name() + ',' + U.name() + ',' + phi.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
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

            if(mesh().objectRegistry::foundObject<surfaceVectorField>("Sf"))
            {
                Info << "ZT, EulerDdtPhiCorr" << endl;

                const surfaceVectorField& Sf = 
                    mesh().objectRegistry::lookupObject<surfaceVectorField>
                    (
                        "Sf"
                    );

                dU0 = (phi.oldTime() - (Sf.oldTime()&dU0))
                    *Sf.oldTime()/magSqr(Sf.oldTime());
            }
            else
            {
                dU0 = (phi.oldTime() - (mesh().Sf()&dU0))
                    *mesh().Sf()/sqr(mesh().magSf());
            }

            return tmp<surfaceScalarField>
            (
                new surfaceScalarField
                (
                    ddtIOobject,
                    rDeltaT*ddtPhiCoeff
                   *fvc::interpolate(rho.oldTime()*V0oV)
                   *(mesh().Sf()&dU0)
                   /fvc::interpolate(1.0/rA) 
                )
            );
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

                return tmp<surfaceScalarField>
                (
                    new surfaceScalarField
                    (
                        ddtIOobject,
                        rDeltaT*ddtPhiCoeff
                       *(
                            fvc::interpolate(rho*rA*V0oV)*(mesh().Sf()&U0)
                          - (
                                fvc::interpolate
                                (
                                    rA*rho*rho.oldTime()*U.oldTime()*V0oV
                                )
                              & mesh().Sf()
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
                       fvc::interpolate(rA*rho*V0oV)*phi.oldTime()
                     - (
                           fvc::interpolate
                           (
                               rA*rho*rho.oldTime()*U.oldTime()*V0oV
                           )
                         & mesh().Sf()
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

                return tmp<surfaceScalarField>
                (
                    new surfaceScalarField
                    (
                        ddtIOobject,
                        rDeltaT*ddtPhiCoeff
                       *(
                            fvc::interpolate(rA*V0oV)*(mesh().Sf()&U0)
                          - (
                                fvc::interpolate(rA*U.oldTime()*V0oV) 
                              & mesh().Sf()
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
                       fvc::interpolate(rA*V0oV)*phi.oldTime()
                     - (
                           fvc::interpolate(rA*U.oldTime()*V0oV) 
                         & mesh().Sf()
                       )
                    )
                )
            );
        }
        else
        {
            FatalErrorIn
            (
                "EulerDdtScheme<vector>::fvcDdtPhiCorr"
            )   << "dimensions of phi are not correct"
                << abort(FatalError);

            return surfaceScalarField::null();
        }
    }
    else
    {
        Info << "ddtPhiCorr fixed mesh" << endl;
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
                        fvc::interpolate(rA*rho.oldTime())*phi.oldTime()
                      - (fvc::interpolate(rA*rho.oldTime()*U.oldTime())
                      & mesh().Sf())
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
                        fvc::interpolate(rA*rho)*phi.oldTime()
                      - (fvc::interpolate(rA*rho*rho.oldTime()*U.oldTime())
                      & mesh().Sf())
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
                        fvc::interpolate(rA)*phi.oldTime()
                      - (fvc::interpolate(rA*U.oldTime()) & mesh().Sf())
                    )
                )
            );
        }
        else
        {
            FatalErrorIn
            (
                "EulerDdtScheme<vector>::fvcDdtPhiCorr"
            )   << "dimensions of phi are not correct"
                << abort(FatalError);

            return surfaceScalarField::null();
        }
    }
}


template<>
tmp<surfaceScalarField> EulerDdtScheme<vector>::meshPhi
(
    const GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    return mesh().phi();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
