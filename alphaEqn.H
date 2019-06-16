plicVofSolver.preProcess();

plicVofSolver.orientation();
plicVofSolver.reconstruction();
plicVofSolver.advection();

rhoPhi = plicVofSolver.getRhoPhi(rho1, rho2);

Info<< "plicVofSolving: Mass conservation Error = "
    << plicVofSolver.massConservationError() << endl;

Info<< "Phase-1 volume fraction = "
    << alpha1.weightedAverage(mesh.Vsc()).value()
    << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
    << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
    << endl;