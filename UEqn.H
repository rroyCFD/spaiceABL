    // Solve the momentum equation
    
    // Update Coriolis forcing
    Coriolis.update();

    // Update buoyancy term
    Boussinesq.updateBuoyancyTerm();

    // Classical CrankNicolson scheme split for advection and dissipation terms
    // - Temporal derivative shall be Euler
    // - Using old-nuEff skipped as noticeable benefit not found +(save storage)
    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(U) 
      // + fvm::div(phi, U)
      // + turbulence->divDevReff(U)
      + 0.5*(fvm::div(phi,U) + fvc::div(phi.oldTime(), U.oldTime())) 
      // - written for eddy viscosity as it's easy to split
      - 0.5*(
            fvm::laplacian(turbulence->nuEff(), U) 
          + fvc::laplacian(turbulence->nuEff(), U.oldTime())
      )
      - fvc::div(turbulence->nuEff() * dev2((fvc::grad(U))().T()))  
      + fvc::div(Rwall)
      - momentumSourceTerm.force()
      - Coriolis.force()
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  - fvc::snGrad(p_rgh)
                  + Boussinesq.buoyancyTerm()
                )*mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
