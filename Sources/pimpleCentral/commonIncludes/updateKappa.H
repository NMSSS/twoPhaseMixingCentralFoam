{
    aphiv_own = phi_own / rho_own;
    aphiv_nei = phi_nei / rho_nei;
    surfaceScalarField rhof = agrp_own*rho_own + agrp_nei*rho_nei;
    surfaceScalarField phiv = phi / rhof;
    
    dimensionedScalar amaxSmall ("amaxSmall", amaxSf.dimensions(), SMALL * min(mesh.magSf()).value());
    surfaceScalarField cSf (max(cSf_own, cSf_nei));
    amaxSf = max(mag(phiv + cSf), mag(phiv - cSf));
    amaxSf.setOriented(true);
    
    kappaFuncPtr().update();
    kappa.primitiveFieldRef() = kappaFuncPtr().kappa()();
    Info << "max/min kappa: " << max(kappa).value() << "/" << min(kappa).value() << endl;
    phi_own += (1.0 - kappa) * phi_nei;
    phi_nei = kappa * phi_nei;
    
}

