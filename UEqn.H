    // Solve the momentum equation
    
    // Update Coriolis forcing
    Coriolis.update();

    // Update buoyancy term
    Boussinesq.updateBuoyancyTerm();

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(U) + fvm::div(phi, U)
      + turbulence->divDevReff(U)
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
