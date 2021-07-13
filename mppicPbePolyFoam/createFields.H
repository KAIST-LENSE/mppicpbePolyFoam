#include "createRDeltaT.H"

Info<< "Creating themophysical properties model\n" << endl;

autoPtr<rhoReactionThermo> Thermo
(
    rhoReactionThermo::New(mesh)
);

rhoReactionThermo& thermo = Thermo();

thermo.validate(args.executable(), "h", "e");

SLGThermo slgThermo(mesh, thermo);

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

const word inertSpecie(thermo.lookup("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species()
        << exit(FatalIOError);
}

const volScalarField& T = thermo.T();
volScalarField Tc("Tc", T);

volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField oP
(
    IOobject
    (
        "oP",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.mu()
);

Info<< "\nReading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(U) & mesh.Sf()
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

volScalarField& p = thermo.p();

p_rgh = p - rho*gh - oP;

pressureControl pressureControl(p, rho, pimple.dict(), false);

mesh.setFluxRequired(p_rgh.name());

dimensionedScalar rhoMax
(
    dimensionedScalar::lookupOrDefault
    (
        "rhoMax",
        pimple.dict(),
        dimDensity,
        GREAT
    )
);

dimensionedScalar rhoMin
(
    dimensionedScalar::lookupOrDefault
    (
        "rhoMin",
        pimple.dict(),
        dimDensity,
        0
    )
);

Info<< "Creating multi-variate interpolation scheme\n" << endl;
multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

Info << "Creating field alphac\n" << endl;
// alphac must be constructed before the cloud
// so that the drag-models can find it
volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimless, 0)
);

volScalarField alphaParticle
(
    IOobject
    (
        "alphaParticle",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimless, 0)
);

volScalarField eps
(
    IOobject
    (
        "eps",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("eps", dimensionSet (0, 2,-3,0,0,0,0), 0)
);

volScalarField alphaMonomer
(
    IOobject
    (
        "alphaMonomer",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimless, 0)
);

volScalarField alphaWater
(
    IOobject
    (
        "alphaWater",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimless, 0)
);

volScalarField R_monomer
(
    IOobject
    (
        "R_monomer",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("R_monomer", dimensionSet (1,-3,-1,0,0,0,0), 0.0)
);

#include "createClouds.H"

// Particle fraction upper limit
scalar alphaMin
(
    1.0
  - readScalar
    (
        kinematicCloud.particleProperties().subDict("constantProperties")
       .lookup("alphaMax")
    )
);

// Update alphac from the particle locations
alpha = max(1.0 - kinematicCloud.theta(), alphaMin);
alpha.correctBoundaryConditions();
//R_monomer = kinematicCloud.R_total() ;
//R_monomer.correctBoundaryConditions();

surfaceScalarField alphaf("alphaf", fvc::interpolate(alpha));
surfaceScalarField rhof("rhof", fvc::interpolate(rho));
surfaceScalarField alphaphi("alphaphi", alphaf*phi);
surfaceScalarField rhophi("rhophi", rhof*phi);
surfaceScalarField alpharhophi("alpharhophi", alphaf*rhof*phi);

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        alpharhophi,
        thermo
    )
);

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

#include "createMRF.H"

