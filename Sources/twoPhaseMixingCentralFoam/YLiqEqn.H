scalarField allFacesLambda(mesh.nFaces(), 1.0);
slicedSurfaceScalarField lambdaCoeffs
(
    IOobject
    (
        "lambdaCoeffs",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    ),
    mesh,
    dimless,
    allFacesLambda,
    false   // Use slices for the couples
);

lambdaCoeffs.primitiveFieldRef() = 1.0;

{
    word YiOwnScheme ("div(phi_own," + YLiq.name() +")");
    word YiNeiScheme ("div(phi_nei," + YLiq.name() +")");
    word YiScheme ("div(phi," + YLiq.name() +")");
    
    surfaceScalarField rhoPhiYCorr
    (
        "rhoPhiYLiqCorr",
        // fvc::flux
        // (
        //     phi_own,
        //     YLiq,
        //     YiOwnScheme
        // )
        // +
        // fvc::flux
        // (
        //     phi_nei,
        //     YLiq,
        //     YiNeiScheme
        // )
        fvc::flux
        (
            phi,
            YLiq,
            YiScheme
        )
    );
    
    surfaceScalarField mDiffC
    (
        "mDiffC",
        fvc::interpolate(turbulence->muEff())
    );
    
    Info << "Creating matrix for fvOptions" << endl;
    
    #warning "change call to fvOptions to compressible"
    
    fvScalarMatrix SuSp
    (
        fvOptions(rho,YLiq)
    );

    volScalarField  Sp
    (
       - rho / runTime.deltaT() * 0
       +
       fvc::laplacian(mDiffC, YLiq)
    );
    
    mulesWithDiffusionImplicitLimiter
    (
        rho,
        YLiq,
        phi,
        allFacesLambda,
        rhoPhiYCorr,
        liqDiffusionFlux,
        mDiffC,
        SuSp,
        Sp
    );
	
    upwind<scalar> UDsPos(mesh, phi_own);
    upwind<scalar> UDsNeg(mesh, phi_nei);
    upwind<scalar> UDs(mesh, phi);

    fvScalarMatrix YLiqEqn
    (
        fvm::ddt(rho,YLiq)
        +
        fv::gaussConvectionScheme<scalar>(mesh, phi, UDs).fvmDiv(phi, YLiq)
        // fv::gaussConvectionScheme<scalar>(mesh, phi, UDsPos).fvmDiv(phi_own, YLiq)
        // +
        // fv::gaussConvectionScheme<scalar>(mesh, phi, UDsNeg).fvmDiv(phi_nei, YLiq)
        // fvc::div(phi_own,YLiq)
        // +
        // fvc::div(phi_nei,YLiq)
        -
        fvm::laplacian(mDiffC, YLiq) 
        +
        fvc::div(lambdaCoeffs * rhoPhiYCorr)
        ==
        SuSp
    );
    
    YLiqEqn.relax();
    
    YLiqEqn.solve(mesh.solver(YLiq.name()));
    

    Info<< "Mass fraction min/max = " << min(YLiq).value() << "/"
        << max(YLiq).value() << endl;
    
    YLiq = max(YLiq, 0.0);
    YLiq = min(YLiq, 1.0);
    YGas = 1.0 - YLiq;
    
}
