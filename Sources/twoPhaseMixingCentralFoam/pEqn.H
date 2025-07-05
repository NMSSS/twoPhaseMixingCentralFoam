//if it is the first pimple / piso iteration
// {
//     p_own = fvc::interpolate(p, own, "reconstruct(p)");
//     p_nei = fvc::interpolate(p, nei, "reconstruct(p)");
// }

volScalarField rhoT (thermo.rho());
rhoHat = rhoT - psi*p;
volVectorField HbyA ("HbyA", U);
volScalarField rAU ("rAU", 1.0 / UEqn.A());
HbyA = UEqn.H() * rAU;
HbyA.boundaryFieldRef() == U.boundaryField();

{
    #include "updateCentralWeights.H"
    #include "rhoUpdateCentralFields.H"

    //surfaceVectorField rhoUf = linearInterpolate(rho*U);
    surfaceScalarField Dpf("Dpf", linearInterpolate(rho*rAU));
    surfaceScalarField phiRhoHbyA
    (
        "phiRhoHbyA",
        fvc::interpolate(rho)*fvc::flux(HbyA) + Dpf*fvc::ddtCorr(rho, U, phi)
    );
    surfaceScalarField onemkappa (1.0 - kappa);
    phiRhoHbyA *= onemkappa;
    Dpf *= onemkappa;

    // blend KNP coeffs with linear interpolation
    Dp_own *= kappa;
    Dp_own += Dpf;
    Dp_nei *= kappa;

    aphiv_own *= kappa;
    aphiv_nei *= kappa;

    phid_own      = aphiv_own * psi_own;
    phid_nei      = aphiv_nei * psi_nei;

    phiRhoHat_own =  aphiv_own * rhoHat_own;
    phiRhoHat_own += phiRhoHbyA;
    phiRhoHat_nei =  aphiv_nei * rhoHat_nei;

    //Update the pressure BCs to ensure flux consistency; added by me
    phiRhoHbyA  = phiRhoHat_own + phiRhoHat_nei + phid_own*p_own + phid_nei*p_nei;
    Dpf         = Dp_own + Dp_nei;
    constrainPressure(p, rhoT, U, phiRhoHbyA, Dpf);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn_nei
        (
            fvm::div(phid_nei,p) + fvc::div(phiRhoHat_nei) - fvm::laplacian(Dp_nei,p)
        );
        
        fvScalarMatrix pEqn_own
        (
            fvm::div(phid_own,p) + fvc::div(phiRhoHat_own) - fvm::laplacian(Dp_own,p)
        );

        fvScalarMatrix pEqn
        (
            fvm::ddt(psi,p)
            +
            fvc::ddt(rhoHat)
            +
            pEqn_nei
            +
            pEqn_own
            ==
            fvOptions(psi,p,rho.name())
        );

        fvOptions.constrain(pEqn);

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi_own = pEqn_own.flux() + phiRhoHat_own;
            phi_nei = pEqn_nei.flux() + phiRhoHat_nei;
            phi = phi_own + phi_nei;
        }
    }

    // surfaceScalarField pf(linearInterpolate(p));
    p_own = fvc::interpolate(p, own, "reconstruct(p)");
    p_nei = fvc::interpolate(p, nei, "reconstruct(p)");

    if (fluxScheme == "HLL" || fluxScheme == "Kurganov" || fluxScheme == "Tadmor")
    {
        //gradp = fvc::div((kappa*(alpha_own*p_own + alpha_nei*p_nei)+onemkappa*pf)*mesh.Sf());
        gradp = fvc::div((alpha_own*p_own + alpha_nei*p_nei)*mesh.Sf());
    }

    #include "limitPressureGradient.H"
}

U = HbyA - rAU * gradp;
U.correctBoundaryConditions();
fvOptions.correct(U);

Info << "max/min(mag(U)): " << max(mag(U)).value() <<"/" << min(mag(U)).value() << endl;
Info << "max/min p: " << max(p).value() << "/" << min(p).value() << endl;

// #include "massEqn.H"
p = max(p, pMin);
thermo.correctRealDensities();
// rhoT = 1.0 / (YLiq / thermo.thermoLiq().rho() + YGas / thermo.thermoGas().rho());
rhoT = rhoHat + psi*p;

#include "massEqn.H"
#include "centralContinuityErrs.H"

rho = rhoT;

rho_own = fvc::interpolate(rho, own, "reconstruct(rho)");
rho_nei = fvc::interpolate(rho, nei, "reconstruct(rho)");
correctCentralACMIInterpolation(rho_nei);


//
//END-OF-FILE
//