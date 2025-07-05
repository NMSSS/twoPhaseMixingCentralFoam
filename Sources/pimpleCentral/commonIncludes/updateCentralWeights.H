    {
        if (noFluxWeighting)
        {
            phiv_own    = (fvc::interpolate(U, own, "reconstruct(U)")) & mesh.Sf();
            phiv_nei    = (fvc::interpolate(U, nei, "reconstruct(U)")) & mesh.Sf();
        }
        else
        {
            rho_own     = fvc::interpolate(rho, own, "reconstruct(rho)");
            rho_nei     = fvc::interpolate(rho, nei, "reconstruct(rho)");
            p_own     = fvc::interpolate(p, own, "reconstruct(p)");
            p_nei     = fvc::interpolate(p, nei, "reconstruct(p)");
            pS_own       = p_own*mesh.magSf()*mesh.magSf();
            pS_nei       = p_nei*mesh.magSf()*mesh.magSf();
            correctCentralACMIInterpolation(rho_nei);
            phiv_own    = (fvc::interpolate(rho*U, own, "reconstruct(U)") / rho_own) & mesh.Sf();
            phiv_nei    = (fvc::interpolate(rho*U, nei, "reconstruct(U)") / rho_nei) & mesh.Sf();
        }

        if (mesh.moving())
        {
            phiv_own   -= mesh_phi;
            phiv_nei   -= mesh_phi;
        }

        c         = sqrt(thermo.Cp()/thermo.Cv() / psi);  //this is only updated for the following analysis
        cf_own      = fvc::interpolate(c, own, "reconstruct(psi)");
        cf_nei      = fvc::interpolate(c, nei, "reconstruct(psi)");
        // surfaceScalarField T_own = fvc::interpolate(T, own, "reconstruct(psi)");
        // surfaceScalarField T_nei = fvc::interpolate(T, nei, "reconstruct(psi)");

        // /*
        // const simplePhasePsiThermo& thermoGas = thermo.thermoGas();
        // Info<< "thermoLiq Gas: " << thermoGas.rho() << endl;
        // Info<< "thermoLiq rho: " << thermo.thermoLiq().rho() << endl;
        // */

        // surfaceScalarField rho_liq_own = fvc::interpolate(thermo.thermoLiq().rho(), own, "reconstruct(rho");
        // surfaceScalarField rho_liq_nei = fvc::interpolate(thermo.thermoLiq().rho(), nei, "reconstruct(rho)");
        // surfaceScalarField rho_gas_own = fvc::interpolate(thermo.thermoGas().rho(), own, "reconstruct(rho)");
        // surfaceScalarField rho_gas_nei = fvc::interpolate(thermo.thermoGas().rho(), nei, "reconstruct(rho)");
        // surfaceScalarField YbarLiqf    = fvc::interpolate(YbarLiq, phi, "reconstruct(YbarLiq)");
        // cf_own = thermo.Cfm(p_own, T_own, YbarLiqf, rho_liq_own, rho_gas_own);
        // cf_nei = thermo.Cfm(p_nei, T_nei, YbarLiqf, rho_liq_nei, rho_gas_nei);

        Info << "max/min cf_own:" << max(cf_own).value() << "/" << min(cf_own).value() << endl;
        Info << "max/min cf_nei:" << max(cf_nei).value() << "/" << min(cf_nei).value() << endl;

        correctCentralACMIInterpolation(cf_nei);
        cSf_own     = cf_own * mesh.magSf();
        cSf_own.setOriented(true);
        cSf_nei     = cf_nei * mesh.magSf();
        cSf_nei.setOriented(true);

        ap = max(max(phiv_own + cSf_own, phiv_nei + cSf_nei), v_zero);
        am = min(min(phiv_own - cSf_own, phiv_nei - cSf_nei), v_zero);

        if (fluxScheme == "Kurganov")
        {
            alpha_own   = ap/(ap - am);
            aSf_own     = am*alpha_own;
            aSf_nei		= -aSf_own;
            alpha_nei   = 1.0 - alpha_own;
        }
        else if (fluxScheme == "Tadmor")
        {
            aSf_own     = -0.5*max(mag(ap), mag(am));
            aSf_nei     = -aSf_own;
			aSf_own.setOriented(true);
			aSf_nei.setOriented(true);
        }
        else if (fluxScheme == "HLL")
        {
           alpha_own   = ap/(ap - am);
           aSf_own     = am*alpha_own;
           alpha_nei   = 1.0 - alpha_own;
           aSf_nei     = -aSf_own;
        } 
        else
        {
            FatalErrorIn
            (
                    "updateCentralWeights.H"
            )   << "Unknown central flux scheme " << fluxScheme << nl << nl
                << "Valid flux schemes are:" << nl
                << "1) Kurganov" << nl
                << "2) Tamodr" << nl
			    << "3) HLL"<<endl;
        }
    }
