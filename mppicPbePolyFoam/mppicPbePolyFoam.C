/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    PBEDPMFoam

Description
    Transient solver for compressible, turbulent flow with a reacting,
    multiphase particle cloud, and optional sources/constraints.

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "rhoReactionThermo.H"
#include "SLGThermo.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

#include "basicKinematicMPPICCloud.H"
#define basicKinematicTypeCloud basicKinematicMPPICCloud

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        if (!LTS)
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        mu = thermo.mu();
        eps = turbulence->epsilon();

        volScalarField x_monomer(Y[0]) ;
	volScalarField x_water(Y[1]) ;

        alphaMonomer = x_monomer / 886.0 / (x_monomer / 886.0 + x_water / 971) ;
        alphaWater = 1 - alphaMonomer ;

        Info<< "Evolving " << kinematicCloud.name() << endl;

        kinematicCloud.evolve();

        // Update continuous phase volume fraction field
        alpha = max(1.0 - kinematicCloud.theta(), alphaMin);
        alphaParticle = 1.0 - alpha ;
        alpha.correctBoundaryConditions();
        alphaf = fvc::interpolate(alpha);
        alphaphi = alphaf*phi;
        rhophi = fvc::interpolate(rho)*phi;
        alpharhophi = alphaf*rhophi;

        fvVectorMatrix cloudSU(kinematicCloud.SU(U));
        volVectorField cloudVolSUSu
        (
            IOobject
            (
                "cloudVolSUSu",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedVector
            (
                "0",
                cloudSU.dimensions()/dimVolume,
                Zero
            ),
            zeroGradientFvPatchVectorField::typeName
        );

        cloudVolSUSu.primitiveFieldRef() = -cloudSU.source()/mesh.V();
        cloudVolSUSu.correctBoundaryConditions();
        cloudSU.source() = Zero;

        R_monomer.correctBoundaryConditions();

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {

            #include "UEqn.H"
            #include "EEqn.H"
            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        rho = thermo.rho();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
