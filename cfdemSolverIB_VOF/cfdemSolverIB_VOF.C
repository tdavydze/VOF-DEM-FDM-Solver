/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright (C) 1991-2009 OpenCFD Ltd.
                                Copyright (C) 2009-2012 JKU, Linz
                                Copyright (C) 2012-     DCS Computing GmbH,Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling.  If not, see <http://www.gnu.org/licenses/>.

Application
    cfdemSolverIB

Description
    Transient solver for incompressible flow.
    The code is an evolution of the solver pisoFoam in OpenFOAM(R) 1.6, 
    where additional functionality for CFD-DEM coupling using immersed body
    (fictitious domain) method is added.
Contributions
    Alice Hager
\*---------------------------------------------------------------------------*/


#include "fvCFD.H"
//#include "singlePhaseTransportModel.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "subCycle.H"
#include "CrankNicolsonDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "pimpleControl.H"
#include "/home/tsimur/OpenFOAM/OpenFOAM-5.x/src/transportModels/incompressible/incompressibleTwoPhaseMixture/incompressibleTwoPhaseMixture.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "interfaceProperties.H"
#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "SLTSDdtScheme.H"
#include <assert.h>


#include "OFversion.H"

#if defined(version30)
    #include "turbulentTransportModel.H"
    //#include "pisoControl.H"
   // #include "createTimeControls.H"
#else
    #include "turbulenceModel.H"
#endif
#include "cfdemCloudIB.H"
//#if defined(superquadrics_flag)
//#include "cfdemCloudIBSuperquadric.H"
//#endif
#include "implicitCouple.H"

#include "averagingModel.H"
#include "voidFractionModel.H"

#include "dynamicFvMesh.H"

#include "cellSet.H"

#if defined(version22)
    #include "meshToMeshNew.H"
    #include "fvIOoptionList.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"

    #include "createDynamicFvMesh.H"

    #if defined(version30)
        //pisoControl piso(mesh);
        #include "createPimpleControl.H"
        #include "createMRF.H"
        #include "createTimeControls.H"
    #endif

    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "createAlphaFluxes.H"
    #include "createFvOptions.H"
//#include "correctPhi.H:"

    #if defined(version22)
        #include "createFvOptions.H"
    #endif

    // create cfdemCloud
    //#include "readGravitationalAcceleration.H"
    #if defined(superquadrics_flag)
        cfdemCloudIBSuperquadric particleCloud(mesh);
    #else
        cfdemCloudIB particleCloud(mesh);
    #endif

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        //=== dyM =================
        interFace = mag(mesh.lookupObject<volScalarField>("voidfractionNext"));
        particleCloud.setMeshHasUpdatedFlag(mesh.update()); //dyM

        #if defined(version30)
            #include "readTimeControls.H"
            #include "alphaCourantNo.H"
            #include "CourantNo.H"
            #include "setDeltaT.H"
        #else
            //#include "readPISOControls.H"
            #include "pimpleControl.H"
            #include "CourantNo.H"
        #endif

        // do particle stuff
        //Info << "- evolve()" << endl;
        particleCloud.evolve(voidfraction, interFace);

        int nCorr = readInt(pimple.dict().lookup("nCorrectors"));
        int nOuterCorr = readInt(pimple.dict().lookup("nOuterCorrectors"));
        int nNonOrthCorr = readInt(pimple.dict().lookup("nNonOrthogonalCorrectors"));

        // Pressure-velocity PISO corrector
        if(particleCloud.solveFlow())
        {
            while (pimple.loop())
            {
                #include "alphaControls.H"
                #include "alphaEqnSubCycle.H"

                mixture.correct();

                // Momentum predictor
                MRF.correctBoundaryVelocity(U);    // not used in unresolved

                fvVectorMatrix UEqn
                (
                    fvm::ddt(rho, U) //fvm::ddt(voidfraction,U)
                  + fvm::div(rhoPhi, U)
                  + MRF.DDt(rho, U)                // Not used in unresolved
                  + turbulence->divDevRhoReff(rho, U)
                    #if defined(version22)
                    ==
                    fvOptions(rho, U)
                    #endif
                );

                UEqn.relax();

                #if defined(version22)
                fvOptions.constrain(UEqn);
                #endif

                #if defined(version30)
                    if (pimple.momentumPredictor())
                #else
                    if (momentumPredictor)
                #endif
                {
                    solve
                    (
                        UEqn
                     ==
                        fvc::reconstruct
                        (
                            (
                                mixture.surfaceTensionForce()
                              - ghf*fvc::snGrad(rho)
                              - fvc::snGrad(pd)
                            ) * mesh.magSf()
                        )
                    );
                    fvOptions.correct(U);
                }

                // --- PISO loop ---
                // pressure correction
                #if defined(version30)
                    while (pimple.correct())
                #else
                    for (int corr=0; corr<nCorr; corr++)
                #endif
                {
                    volScalarField rUA = 1.0/UEqn.A();
                    surfaceScalarField rUAf(fvc::interpolate(rUA));

                    volVectorField HbyA(constrainHbyA(rUA*UEqn.H(), U, pd));    // NEW !!

                    #ifdef version23
                    surfaceScalarField phiHbyA              //NEW !!
                    (
                        "phiHbyA",
                        fvc::flux(HbyA)
                       +fvc::interpolate(rho*rUA)*fvc::ddtCorr(U,phi)
                    );

                    #else
                    surfaceScalarField phiHbyA
                    (
                        "phiHbyA",
                        fvc::flux(HbyA)
                       +fvc::interpolate(rho*rUA)*fvc::ddtPhiCorr(rUA, U, phi)
                    );

                    #endif
                    MRF.makeRelative(phiHbyA);

                    adjustPhi(phiHbyA, U, pd);

                    surfaceScalarField phig
                    (
                        (
                            mixture.surfaceTensionForce()
                           -ghf*fvc::snGrad(rho)
                        )*rUAf*mesh.magSf()

                    );

                    //phi = phiHbyA + phig;
                    phiHbyA += phig;

                    //update the pressure BC to ensure flux consistency
                    constrainPressure(pd,U, phiHbyA, rUAf, MRF);

                    #if defined(version22)
                    fvOptions.relativeFlux(phi);
                    #endif


                    // Non-orthogonal pressure corrector loop
                    #if defined(version30)
                        while (pimple.correctNonOrthogonal())
                    #else
                        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                    #endif
                    {
                        // Pressure corrector

                        fvScalarMatrix pdEqn
                        (
                            fvm::laplacian(rUAf, pd) == fvc::div(phiHbyA)
                        );

                        pdEqn.setReference(pRefCell, getRefCellValue(pd, pRefValue));

                        #if defined(version30)
                            pdEqn.solve(mesh.solver(pd.select(pimple.finalInnerIter())));
                            if (pimple.finalNonOrthogonalIter())
                            {
                                phi = phiHbyA - pdEqn.flux();
                                pd.relax();
                                U = HbyA + rUA*fvc::reconstruct((phig - pdEqn.flux())/rUAf);
                            }

                        #else
                            if( corr == nCorr-1 && nonOrth == nNonOrthCorr )
                                #if defined(versionExt32)
                                    pdEqn.solve(mesh.solutionDict().solver("pdFinal"));
                                #else
                                    pdEqn.solve(mesh.solver("pdFinal"));
                                #endif

                            else
                                pdEqn.solve();

                            if (nonOrth == nNonOrthCorr)       // just final pd loop
                            {
                                phi = phiHbyA - pdEqn.flux();
                                pd.relax();
                                U = HbyA + rUA*fvc::reconstruct((phig - pdEqn.flux())/rUAf);
                            }
                        #endif
                    } //end correct nonOrthogonal

                } // end nCorrect/pimple.correct

                //U -= rUA*fvc::grad(p);
                U.correctBoundaryConditions();
                fvOptions.correct(U);

            } //end solveFlow

            #include "continuityErrs.H"

            p == pd + rho*gh;

            mixture.correct();
            turbulence->correct();
        }  // end of particleCloud.solveFlow


        Info << "particleCloud.calcVelocityCorrection() " << endl;
        volScalarField voidfractionNext=mesh.lookupObject<volScalarField>("voidfractionNext");

        particleCloud.calcVelocityCorrection(p,U,phiIB,voidfractionNext);
        #if defined(version22)
        fvOptions.correct(U);
        #endif

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
