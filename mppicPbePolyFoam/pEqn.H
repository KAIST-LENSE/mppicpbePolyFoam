if (!pimple.simpleRho())
{
    rho = thermo.rho();
}
    rhof = fvc::interpolate(rho);
    const volScalarField psip0(psi*p);
    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

    surfaceScalarField phig(-alphaf*rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());
    surfaceScalarField phiForce(-rhorAUf*fvc::flux(cloudVolSUSu));

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
           fvc::interpolate(rho)*fvc::flux(HbyA)
         + MRF.zeroFilter(alphaf*rhorAUf*fvc::ddtCorr(rho, U, alpharhophi))
        )
        + phig
        + phiForce
    );

    MRF.makeRelative(fvc::interpolate(alpha*rho), phiHbyA);
    if (p_rgh.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p_rgh);
        fvc::makeAbsolute(phiHbyA, U);
    }

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, MRF);

    fvScalarMatrix pDDtEqn
    (
        fvc::ddt(alpha, rho)
      + alpha*psi*correction(fvm::ddt(p_rgh))
      + fvc::div(alphaf*phiHbyA)
      ==
      - R_monomer
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn(pDDtEqn - fvm::laplacian(alphaf*rhorAUf, p_rgh));

        p_rghEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi = (phiHbyA + p_rghEqn.flux()/alphaf)/rhof;
            alphaphi = alphaf*phi;
            rhophi = rhof*phi;
            alpharhophi = alphaf*rhof*phi;
            p_rgh.relax();
            U = HbyA + rAU*fvc::reconstruct((p_rghEqn.flux() + phig + phiForce)/rhorAUf);
            U.correctBoundaryConditions();
        }
    }

// Thermodynamic density update
thermo.correctRho(psi*p - psip0);

//#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

K = 0.5*magSqr(U);

p_rgh.correctBoundaryConditions();

p == p_rgh + rho*gh + oP;

if (pressureControl.limit(p))
{
    p_rgh = p - rho*gh - oP;
}

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
}
/*
rho = thermo.rho();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
Info<< "rho min/max = " << min(rho).value() << ", " << max(rho).value() << endl;
*/
