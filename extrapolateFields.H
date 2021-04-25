// Extrapolate fields using Adam-Bashforth like extrapolation method
if(adjustTimeStep)
{
	extrapCoeff =  0.5 * (runTime.deltaTValue()/runTime.deltaT0Value());
	// bound between 0 and 1; (assuming dt(n+1) is no larger than twice of dt(n))
	extrapCoeff = max(min(extrapCoeff, 1.0), 0.0);
	Info << "Extrapolation coeffcient: " << extrapCoeff << endl;
}

// extrapolate fields to next-time step to obtain a better initial condition
p_rgh = ((1 + extrapCoeff) * p_rgh.oldTime() - extrapCoeff * p_rgh.oldTime().oldTime());
if(extrapolP_rgh_BC)
{
    p_rgh.correctBoundaryConditions();
}

// Update U and phi with their respective extrapolated value
phi = ((1 + extrapCoeff) * phi.oldTime() - extrapCoeff * phi.oldTime().oldTime());

U = ((1 + extrapCoeff) * U.oldTime() - extrapCoeff * U.oldTime().oldTime());
// U.correctBoundaryConditions(); updated at UEqn matrix formation

// T = ((1 + extrapCoeff) * T.oldTime() - extrapCoeff * T.oldTime().oldTime());
// T.correctBoundaryConditions();
