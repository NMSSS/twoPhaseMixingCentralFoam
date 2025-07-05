{

    // rho = 1.0 / (YLiq / thermo.thermoLiq().rho() + YGas / thermo.thermoGas().rho());
    volScalarField& he = thermo.he();
    volScalarField alphahEff (turbulence->alphaEff());

    volScalarField Cp (thermo.Cp());
    // Info<< "Cp min/max = " << min(Cp).value() << "/"
    //     << max(Cp).value() << endl;
    upwind<scalar> UDs(mesh, phi);

    fvScalarMatrix hConvection
    (
        fvm::ddt(rho,he)
        // +
        // fv::gaussConvectionScheme<scalar>(mesh, phi, UDs).fvmDiv(phi, he)
        +
        fvm::div(phi_own,he)
        +
        fvm::div(phi_nei,he)

    );
    
    /*
     *
     * Explicit flux correction, bounded with MULES
     *
     */
    // surfaceScalarField hFluxCorr
    // (
    //     lambdaCoeffs
    //     *
    //     (
    //         fvc::flux
    //         (
    //             phi,
    //             he,
    //             "div(phi," + he.name() + ")"
    //         )
    //         -
    //         hConvection.flux()
    //     )
    // );

    surfaceScalarField hCmptFlux
    (
        "hCmptFlux",
        fvc::flux
        (
            fvc::snGrad(YLiq) * mesh.magSf(),
            thermo.thermoLiq().he() * alphahEff,
            "div(alphaEff*hi,grad(Yi))"
        )
        +
        fvc::flux
        (
            fvc::snGrad(YGas) * mesh.magSf(),
            thermo.thermoGas().he() * alphahEff,
            "div(alphaEff*hi,grad(Yi))"
        )
    );
    fvScalarMatrix EEqn
    (
        hConvection  + fvc::div(hCmptFlux) //+ fvc::div(hFluxCorr)
        +
        KChange
        ==
        fvm::laplacian(alphahEff, he)
        +
        KSource
        +
        fvOptions(rho,he)
    );
    fvOptions.constrain(EEqn);
    EEqn.solve();
    fvOptions.correct(he);
    thermo.correct();
}
