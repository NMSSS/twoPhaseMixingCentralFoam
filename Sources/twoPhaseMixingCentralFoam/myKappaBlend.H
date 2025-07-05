{
    aphiv_own = phi_own / rho_own;
    aphiv_nei = phi_nei / rho_nei;
	
    // surfaceScalarField rhof = agrp_own*rho_own + agrp_nei*rho_nei;
	// surfaceScalarField rhof = alpha_own*rho_own + alpha_nei*rho_nei;
    // surfaceScalarField phiv = phi / rhof;
    
    // dimensionedScalar amaxSmall ("amaxSmall", amaxSf.dimensions(), SMALL * min(mesh.magSf()).value());
    // surfaceScalarField cSf (max(cSf_own, cSf_nei));
    // amaxSf = max(mag(phiv + cSf), mag(phiv - cSf));
    // amaxSf.setOriented(true);
    // YbarLiqf    = fvc::interpolate(YbarLiq, phi, "reconstruct(YbarLiq)");
    if (fluxScheme == "HLL" || fluxScheme == "Kurganov" || fluxScheme == "Tadmor")
    {
        agrp_own = alpha_own;
        agrp_nei = alpha_nei;
    }
    kappaFuncPtr().update();
    kappa.primitiveFieldRef() = kappaFuncPtr().kappa()();
    kappa = kappaFuncPtr().kappa()();

    kappa.setOriented(false);
    //const scalar globalKappa = 1.0; //1.0, 0.9, 0.5, 0.25;
    //kappa.primitiveFieldRef() = globalKappa; // kappaFuncPtr().kappa()();
    //kappa.boundaryFieldRef() =  globalKappa;
    
    //Info << "max/min kappa: " << max(kappa).value() << "/" << min(kappa).value() << endl;
}

