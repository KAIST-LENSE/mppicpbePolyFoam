fvVectorMatrix UEqn
(
    fvm::ddt(alpha, rho, U) + fvm::div(alpharhophi, U)
  + MRF.DDt(alpha*rho, U)
  + turbulence->divDevRhoReff(U)
 ==
    cloudSU
  - R_monomer*U
);

UEqn.relax();

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
      ==
        fvc::reconstruct
        (
                -fvc::flux(cloudVolSUSu)+
            (
                -alphaf*ghf*fvc::snGrad(rho)
                -fvc::snGrad(p_rgh)
            )*mesh.magSf()
        )
    );

    K = 0.5*magSqr(U);
}
