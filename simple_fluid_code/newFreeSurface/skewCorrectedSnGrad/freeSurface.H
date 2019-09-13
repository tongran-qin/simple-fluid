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

Class
    freeSurface

Description
    Implementation of Free surface (interface) tracking method 
    based on moving mesh.

SourceFiles
    freeSurface.C
    freeSurfacePointDisplacement.C
    makeFreeSurfaceData.C

\*---------------------------------------------------------------------------*/

#ifndef FreeSurface_H
#define FreeSurface_H


#include "fvCFD.H"
#include "IOdictionary.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "faCFD.H"
#include "IOpatchToPatchInterpolation.H"
#include "dynamicFvMesh.H"

#include "surfactantProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class freeSurface Declaration
\*---------------------------------------------------------------------------*/

class freeSurface
:
    public IOdictionary
{
    // Private data
    
        //- Reference to fvMesh 
        dynamicFvMesh& mesh_;

        //- Reference to density field
        const volScalarField& rho_;
        volScalarField& rho1_;
        volScalarField& nt_;
        //- Reference to velocity field
        volVectorField& U_;

        //- Reference to pressure field
        volScalarField& p_;

        ///RG-TQ Reference to temperature field
        volScalarField& T_;

        //- Reference to pressure field
//        volScalarField& conLiquid_;

        ///RG-TQ Reference to temperature field
        volScalarField& con_;

        //- Reference to fluid flux field
        const surfaceScalarField& phi_;

        //- Current time index for free-surface
        label curTimeIndex_;

        //- Are there exist two fluids in the free-surface model
        Switch twoFluids_;

        //- Is it free-surface points displacement direction
        //  parallel with free-surface point normals
        Switch normalMotionDir_;

        //- Is it free-surface clean of surfactant
        Switch cleanInterface_;

        //- ID of the fvMesh patch which represents 
        //  A side of the free-surface
        label aPatchID_; 

        //- ID of the fvMesh patch which represents 
        //  B side of the free-surface
        label bPatchID_; 

        //- Viscosity of fluid A
        dimensionedScalar muFluidA_;
    
        //- Viscosity of fluid B
        dimensionedScalar muFluidB_;
        dimensionedScalar muFluidB1_;
        //- Density of fluid A
        dimensionedScalar rhoFluidA_;
    
        //- Density of fluid B
        dimensionedScalar rhoFluidB_;

        //- Density of fluid B
        dimensionedScalar rhoFluidB1_;

        ///RG-TQ Heat conductivity of fluid A
        dimensionedScalar kFluidA_;
    
        ///- Heat conductivity of fluid B
        dimensionedScalar kFluidB_;
        dimensionedScalar kFluidB1_;
        ///RG-TQ Heat conductivity of Wall
        dimensionedScalar kWall_;

        ///RG-TQ Thickness of Wall
        dimensionedScalar thWall_;

        dimensionedScalar TLeft_;
        dimensionedScalar TRight_;
        dimensionedScalar acmCoef_;

        ///- Heat capacity of fluid A
        dimensionedScalar CpFluidA_;
    
        ///- Heat capacity of fluid B
        dimensionedScalar CpFluidB_;
        dimensionedScalar CpFluidB1_;
        ///- Volumetric Thermal expansion of fluid A
        dimensionedScalar betaFluidA_;
    
        ///- Volumetric Thermal expansion of fluid B
        //dimensionedScalar betaFluidB_;

        ///- Diffusion coeff of fluid A
        dimensionedScalar DfFluidA_;
    
        ///- Diffusion coeff of fluid B in air
        dimensionedScalar DfFluidB_;
        dimensionedScalar DfFluidB0_;

        dimensionedScalar antoineA_;
        dimensionedScalar antoineB_;
        dimensionedScalar antoineC_;

        dimensionedScalar pTotal_;

        dimensionedScalar c0_;
        ///- Absolute pressure offset
        dimensionedScalar pCorr_;

        dimensionedScalar MB1_;
        dimensionedScalar MBD_;
        dimensionedScalar Mair0_;

        ///- Initial Volume of Fluid
        dimensionedScalar VvInitial_;
        dimensionedScalar MvInitial_;
        ///- Instant Volume of Fluid
        dimensionedScalar Vv_;
        dimensionedScalar Mv_;
        dimensionedScalar molarMass1_;
        dimensionedScalar molarMassD_;

        ///- Latent heat of evaporation
        dimensionedScalar latentHeat_;

        ///- Universal gas constant
        dimensionedScalar gasConstant1_;

        ///- Universal gas constant
        dimensionedScalar gasConstantD_;

        ///- Reference Temperature and reference vapor density
        dimensionedScalar TRef_;
//        dimensionedScalar rhoVaporRef_;

        //- Reference to gravitational acceleration
        dimensionedVector g_;

        ///- Contact Angle
        dimensionedScalar ctangle_;

        //- Surface tension for the clean free-surface
        dimensionedScalar cleanInterfaceSurfTension_;

        ///RG-TQ Temperature coefficient of surface tension
        dimensionedScalar tempCoeffSurfTension_;

        ///- Phase change mass flux field
        scalarField J;
        scalarField Vn;
        ///- Interfacial temperature
        scalarField TFs;

        ///- Vapor density corresponds to equilibrium state at TFs
        scalarField rhoB1;
        scalarField ncB1;
        scalarField ncB1b;
        scalarField ncBD;
        scalarField ncBDb;
        scalarField gradn1B;
        ///- Initialization of correction pointfield
        pointField displacementInitial;

        //- Free surface patches which do not move
        wordList fixedFreeSurfacePatches_;

        //- Free surface patches for witch point normals 
        //  must be corrected
        wordList pointNormalsCorrectionPatches_;
//z0312s
        int nFreeSurfCorr_;
//z0312e
        ///RG-TQ Deviation of contact angle
//        scalar dTheta_;

        //- Interface smoothing at the begining of time step
        Switch smoothing_;
//z0312s
        Switch     correctPointNormals_;
        Switch     correctDisplacement_;
        Switch     correctCurvature_;
        Switch     curvExtrapOrder_;
//z0312e
        //- Calc normal derivative of normal velocity using fvc
        Switch fvcNGradUn_;//0917-2012

    // Demand-driven data

        //- Patch to patch interpolation object which deals with
        //  interpolation of flow variables between two free-surface 
        //  (interface) patches
        mutable IOpatchToPatchInterpolation* interpolatorABPtr_;
    
        //- Patch to patch interpolation object which deals with
        //  interpolation of flow variables between two free-surface 
        //  (interface) patches
        mutable IOpatchToPatchInterpolation* interpolatorBAPtr_;

        //- Points which are attached to the free-surface A side faces

        //  and which defines the free-surface shape
        mutable vectorIOField* controlPointsPtr_;
        mutable areaScalarField* TiPtr_;
        mutable areaScalarField* TiDiffPtr_;
        mutable areaScalarField* dTiPtr_;


        mutable areaVectorField* gradTiPtr_;
        mutable areaScalarField* massFluxPtr_;
        mutable areaScalarField* heatFluxlPtr_;
        mutable areaScalarField* heatFluxvPtr_;
        mutable areaVectorField* UiPtr_;
        mutable areaScalarField* conciPtr_;

        //- Field which additionally determines 
        //  the motion of free-surface points
        mutable labelList* motionPointsMaskPtr_;

        //- Displacement direction of free-surface points
        mutable vectorField* pointsDisplacementDirPtr_;

        //- Displacement direction of free-surface control points
        mutable vectorField* facesDisplacementDirPtr_;

        //- Total displacement of free-surface points 
        //  in one time step
        mutable vectorIOField* totalDisplacementPtr_;

        //- Finite area mesh which consists
        //  free-surface side A faces
        mutable faMesh* aMeshPtr_;

        //- Free-surface velocity field
        mutable areaVectorField* UsPtr_;

        //- Free-surface fluid flux
        mutable edgeScalarField* phisPtr_;

        //- Free-surface flux
        mutable areaScalarField* phiPtr_;

        //- Free-surface flux
        mutable areaScalarField* ddtPhiPtr_;

        //- Free-surface surfactant concetration
        mutable areaScalarField* surfactConcPtr_;

        //- Surface tension field
        mutable areaScalarField* surfaceTensionPtr_;

        //- Surfactant properties
        mutable surfactantProperties* surfactantPtr_;

        //- Fluid indicator
        mutable volScalarField* fluidIndicatorPtr_;

        //- Contact angle
//        mutable edgeVectorField* contactAnglePtr_;
        mutable edgeScalarField* contactAnglePtr_;

        mutable scalarField* nGradUnPtr_;//0917-2012

    // Private Member Functions

        // Make demand-driven data

            void makeInterpolators();
            void makeControlPoints();

            void makeTi();
            void makeTiDiff();
            void makedTi();            
            void makegradTi();  
            void makemassFlux(); 
            void makeheatFluxl(); 
            void makeheatFluxv(); 
            void makeUi(); 
            void makeconci();               
            
            void makeMotionPointsMask();
            void makeDirections();
            void makeTotalDisplacement();
            void readTotalDisplacement();
            void makeFaMesh() const;
            void makeUs() const;
            void makePhis();
            void makePhi();
            void makeDdtPhi();
            void makeSurfactConc() const;
            void makeSurfaceTension() const;
            void makecontactAngle() const;
            void makeSurfactant() const;
            void makeFluidIndicator();
            void makeNGradUn() const;//0917-2012

        //- Clear all demand-driven data
        void clearOut();

        //- Disallow default bitwise copy construct
        freeSurface(const freeSurface&);

        //- Disallow default bitwise assignment
        void operator=(const freeSurface&);

        //- Move control ponts for deltaH and calculate
        //  free surface points displacement
        //  for given new control points position
        tmp<vectorField> pointDisplacement(const scalarField& deltaH);

        // Calc least sqare plane point and normal
        tmp<vectorField> lsPlanePointAndNormal
        (
            const vectorField& points,
            const vector& origin,
            const vector& axis
        ) const;

        // Set surface normals at solid walls according to contact angle
        void correctContactLinePointNormals();

        // Calc surface gradient of surface scalar field
        tmp<areaVectorField> grad(areaScalarField& Ts);

public:

    // Declare name of the class and it's debug switch
    ClassName("freeSurface");


    // Constructors

    freeSurface
    (
        dynamicFvMesh& m,
        const volScalarField& rho,
        volScalarField& rho1b,
        volScalarField& ntb,
        volVectorField& Ub, 
        volScalarField& Pb,
        ///RG-TQ
        volScalarField& Tb,
        volScalarField& conb,
        const surfaceScalarField& phi
    );


    // Destructor

    ~freeSurface();  


    // Member Functions


        //- Return reference to fvMesh
        dynamicFvMesh& mesh()
        {
            return mesh_;
        }

        //- Return reference to fvMesh
        const dynamicFvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return reference to DB
        const Time& DB() const
        {
            return U_.mesh().time();
        }

        //- Return reference to velocity field
        const volScalarField& rho() const
        {
            return rho_;
        }

        volScalarField& rho1()
        {
            return rho1_;
        }
        const volScalarField& rho1() const
        {
            return rho1_;
        }

        //- Return reference to velocity field

        volScalarField& nt()
        {
            return nt_;
        }
        const volScalarField& nt() const
        {
            return nt_;
        }

        //- Return reference to velocity field
        volVectorField& U()
        {
            return U_;
        }

        //- Return reference to velocity field
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return reference to pressure field
        volScalarField& p()
        {
            return p_;
        }

        //- Return reference to pressure field
        const volScalarField& p() const
        {
            return p_;
        }

        ///RG-TQ- Return reference to temperature field
        volScalarField& T()
        {
            return T_;
        }

        ///RG-TQ- Return reference to temperature field
        const volScalarField& T() const
        {
            return T_;
         }

        ///RG-TQ- Return reference to temperature field
        volScalarField& con()
        {
            return con_;
        }

        ///RG-TQ- Return reference to temperature field
        const volScalarField& con() const
        {
            return con_;
        }
         //- Return reference to flux field
         const surfaceScalarField& phi() const
         {
             return phi_;
         }

        //- Are there exist two fluid
        //  in the free-surface model
        bool twoFluids() const
        {
            return twoFluids_;
        }

        //- Motion direction swithc
        bool normalMotionDir() const
        {
            return normalMotionDir_;
        }

        //- Is it interface clean
        bool cleanInterface() const
        {
            return cleanInterface_;
        }

        //- ID of patch which represent A side
        //  of free-surface
        const label& aPatchID() const
        {
            return aPatchID_;
        }

        //- ID of patch which represent B side
        //  of free-surface
        const label& bPatchID() const
        {
            return bPatchID_;
        }

        //- Viscosity of fluid A
        const dimensionedScalar& muFluidA() const
        {
            return muFluidA_;
        }

        //- Viscosity of fluid B
        const dimensionedScalar& muFluidB() const
        {
            return muFluidB_;
        }
        const dimensionedScalar& muFluidB1() const
        {
            return muFluidB1_;
        }
        //- Density of fluid A
        const dimensionedScalar& rhoFluidA() const
        {
            return rhoFluidA_;
        }

        //- Densiti of fluid B
        const dimensionedScalar& rhoFluidB() const
        {
            return rhoFluidB_;
        }

        //- Densiti of fluid B
        const dimensionedScalar& rhoFluidB1() const
        {
            return rhoFluidB1_;
        }


        ///- Heat conductivity of fluid A
        const dimensionedScalar& kFluidA() const
        {
            return kFluidA_;
        }

        ///- Heat conductivity of fluid B
        const dimensionedScalar& kFluidB() const
        {
            return kFluidB_;
        }
        const dimensionedScalar& kFluidB1() const
        {
            return kFluidB1_;
        }
        const dimensionedScalar& kWall() const
        {
            return kWall_;
        }
        const dimensionedScalar& thWall() const
        {
            return thWall_;
        }
        const dimensionedScalar& TLeft() const
        {
            return TLeft_;
        }
        const dimensionedScalar& TRight() const
        {
            return TRight_;
	}
        const dimensionedScalar& acmCoef() const
        {
            return acmCoef_;
        }
        ///- Heat capacity of fluid A
        const dimensionedScalar& CpFluidA() const
        {
            return CpFluidA_;
        }

        ///- Heat capacity of fluid B
        const dimensionedScalar& CpFluidB() const
        {
            return CpFluidB_;
        }
        const dimensionedScalar& CpFluidB1() const
        {
            return CpFluidB1_;
        }
        ///- Volumetric thermal expansion of fluid A
        const dimensionedScalar& betaFluidA() const
        {
            return betaFluidA_;
        }

        ///- Volumetric thermal expansion of fluid B
        /*const dimensionedScalar& betaFluidB() const
        {
            return betaFluidB_;
        }*/

        ///- Heat capacity of fluid A

        const dimensionedScalar& DfFluidA() const
        {
            return DfFluidA_;
        }

        ///- Heat capacity of fluid B
        const dimensionedScalar& DfFluidB() const
        {
            return DfFluidB_;
        }
        const dimensionedScalar& DfFluidB0() const
        {
            return DfFluidB0_;
        }

        const dimensionedScalar& antoineA() const
        {
            return antoineA_;
        }
        const dimensionedScalar& antoineB() const
        {
            return antoineB_;
        }
        const dimensionedScalar& antoineC() const
        {
            return antoineC_;
        }

        const dimensionedScalar& pTotal() const
        {
            return pTotal_;
        }
        const dimensionedScalar& c0() const
        {
            return c0_;
        }

        ///- Absolute pressure offset
        const dimensionedScalar& pCorr() const
        {
            return pCorr_;
        }

        const dimensionedScalar& MB1() const
        {
            return MB1_;
        }
        const dimensionedScalar& MBD() const
        {
            return MBD_;
        }
        const dimensionedScalar& Mair0() const
        {
            return Mair0_;
        }

        ///- Initial water volume
        const dimensionedScalar& VvInitial() const
        {
            return VvInitial_;
        }
        ///- Instant water volume
        const dimensionedScalar& Vv() const
        {
            return Vv_;
        }
        ///- Initial water volume
        const dimensionedScalar& MvInitial() const
        {
            return MvInitial_;
        }
        ///- Instant water volume
        const dimensionedScalar& Mv() const
        {
            return Mv_;
        }

        const dimensionedScalar& molarMass1() const
        {
            return molarMass1_;
        }
        const dimensionedScalar& molarMassD() const
        {
            return molarMassD_;
        }

        ///- Latent heat
        const dimensionedScalar& latentHeat() const
        {
            return latentHeat_;
        }

        ///- Universal gas constant
        const dimensionedScalar& gasConstant1() const
        {
            return gasConstant1_;
        }

        ///- Universal gas constant
        const dimensionedScalar& gasConstantD() const
        {
            return gasConstantD_;
        }

        ///- Reference Temperature
        const dimensionedScalar& TRef() const
        {
            return TRef_;
        }
//        const dimensionedScalar& rhoVaporRef() const
//        {
//            return rhoVaporRef_;
//        }

        ///- Contact Angle
        const dimensionedScalar& ctangle() const
        {
            return ctangle_;
        }

        //- Gravity
        const dimensionedVector& g() const
        {
            return g_;
        }

        //- Surface tension for clean interface
        const dimensionedScalar& cleanInterfaceSurfTension() const
        {
            return cleanInterfaceSurfTension_;
        }

        ///RG-TQ- Temperature coefficient of surface tension
        const dimensionedScalar& tempCoeffSurfTension() const
        {
            return tempCoeffSurfTension_;
        }

        ///RG-TQ Deviation of contact angle
//        const scalar& dTheta() const
//        {
//            return dTheta_;
//        }

        //- Return reference to interpolator from free surface side A to
        //  free surface side B
        const IOpatchToPatchInterpolation& interpolatorAB();

        //- Return reference to interpolator from free surface side B to
        //  free surface side A
        const IOpatchToPatchInterpolation& interpolatorBA();

        //- Return control points
        vectorField& controlPoints();

        areaScalarField& Ti();
        areaScalarField& TiDiff();
        areaScalarField& dTi();
	areaVectorField& gradTi();
        areaScalarField& massFlux();
        areaScalarField& heatFluxl();
        areaScalarField& heatFluxv();
	areaVectorField& Ui();
        areaScalarField& conci();
        
        //- Return reference to motion points mask field
        labelList& motionPointsMask();

        //- Return reference to point displacement direction field
        vectorField& pointsDisplacementDir();

        //- Return reference to control points displacement direction field
        vectorField& facesDisplacementDir();

        //- Return total points displacement
        vectorField& totalDisplacement();

        //- Return reference to finite area mesh
        faMesh& aMesh();

        //- Return reference to finite area mesh
        const faMesh& aMesh() const;

        //- Return free-surface velocity field
        areaVectorField& Us();

        //- Return free-surface velocity field
        const areaVectorField& Us() const;

        //- Return free-surface fluid flux field
        edgeScalarField& Phis();

        //- Return free-surface fluid flux field
        areaScalarField& phi();

        //- Return free-surface fluid flux field
        areaScalarField& ddtPhi();

        //- Return free-surface surfactant concentration field
        areaScalarField& surfactantConcentration();

        //- Return free-surface surfactant concentration field
        const areaScalarField& surfactantConcentration() const;

        //- Return surface tension field
        areaScalarField& surfaceTension();

        //- Return surface tension field
        const areaScalarField& surfaceTension() const;

        //- Return contact angle field
        edgeScalarField& contactAngle();
        //- Return contact angle field
        const edgeScalarField& contactAngle() const;

        //- Return surface tension gradient
        tmp<areaVectorField> surfaceTensionGrad();

        //- Return surfactant properties
        const surfactantProperties& surfactant() const;

        //- Return referenct to fluid indicator
        const volScalarField& fluidIndicator();

        //- Return normal derivative of normal velicty at the interface
        scalarField& nGradUn();

        //- Return normal derivative of normal velicty at the interface
        const scalarField& nGradUn() const;

        //- Update control points end displacement directions
        void updateDisplacementDirections();

        //- Predict interface points position
        bool predictPoints();

        //- Correct interface points position
        bool correctPoints();

        //- Move only free-surface points
        bool movePoints(const scalarField& interfacePHi);

        ///RG-TQ- Move the whole plane to make sure volume is conserved
        bool moveWholeFreeSurface(dimensionedScalar,dimensionedScalar);

        //- Move mesh points for old time free-surface points displacement
        bool moveMeshPointsForOldFreeSurfDisplacement();

        //- Move mesh points
        bool moveMeshPoints();

        //- Update free-surface boundary conditions
        void updateBoundaryConditions();

        ///RG-TQ- Correct free-surface temperature
        void updateTemperature();

        //- Correct normal derivative of normal velocity
        void updateNGradUn();

        //- Correct free-surface pressure
        void updatePressure();

        //- Correct free-surface velocity
        void updateVelocity();

        ///RG-TQ- Update free-surface mass flux due to phase change
        void updateMassFlux();

        ///RG-TQ- Update free-surface concentration due to phase change
        void updateCon();

        //- Update free-surface flux
        void updateSurfaceFlux();

        //- Update free-surface surfactant concentration 
        void updateSurfactantConcentration();

        //- Correct surface velocity boundary conditions
        void correctUsBoundaryConditions();

        ///- Calculate free-surface undulation indicator
        tmp<scalarField> divPhiSUndulationIndicator();
        tmp<scalarField> gradGradDivPhiS();

        ///- Calculate geometric mean position of the interface
        dimensionedVector geometricMeanPosition() const;

        //- Calculate total pressure force
        vector totalPressureForce() const;

        //- Calculate total viscous force
        vector totalViscousForce() const;

        //- Calculate total surface tension force
        vector totalSurfaceTensionForce() const;

        ///- Calculate free-surface undulation indicator
        tmp<scalarField> undulationIndicator();

        ///- Smooth free-surface
        void smooth();

        ///- Chech free-surface faces flatness
        tmp<scalarField> checkFaceFlatness();

        //- Initialize control points position
        void initializeControlPointsPosition();
//z0312s
        void smoothing();
//z0312e
        //- Maximal surface tension based Courant number
        scalar maxCourantNumber();

        ///RG-TQ- Update absolute pressure offset and vapor density
//        void updatePressureDensity(dimensionedScalar, dimensionedScalar,dimensionedScalar);

//        void updatePressureDensity(dimensionedScalar,dimensionedScalar);
        void updatePressureDensity(dimensionedScalar);
        dimensionedScalar updatePressureDensity2();
        ///RG-TQ- Update water volume and initial water volume
        void updateVolume(dimensionedScalar,dimensionedScalar);

	void initializeConcentration();
	void initializeMass();

        //- Update properties
        void updateProperties();

        //~- Write VTK freeSurface mesh
        void writeVTK() const;

        //~- Write VTK freeSurface control points
        void writeVTKControlPoints();

//z0312s
        void correctCurvature();

	void correctPointNormals();

	void correctPointDisplacement
	(
    	const scalarField& , vectorField& 
	);

	void smoothCurvature();

//z0312e


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

