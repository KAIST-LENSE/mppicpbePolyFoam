{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(alpha, rho, he) + fvm::div(alpharhophi, he)
      + fvc::ddt(alpha, rho, K) + fvc::div(alpharhophi, K)
      - fvm::laplacian(alpha*turbulence->alphaEff(), he)
      + (
            he.name() == "e"
          ? fvc::div(fvc::absolute(alphaphi, U), p, "div(phiv,p)")
          : -alpha*dpdt
        )
    ==
        alpha*rho*(U&g)
    );
    EEqn.relax();

    EEqn.solve();

    thermo.correct();

    Info<< "T min/max   " << min(T).value() << ", "
        << max(T).value() << endl;
}

