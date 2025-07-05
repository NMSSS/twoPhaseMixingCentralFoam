/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
-------------------------------------------------------------------------------
       hybridCentralSolvers | Copyright (C) 2016-2018 ISP RAS (www.unicfd.ru)
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "standardAndSonicInverseMachKappaFunction.H"
#include "fvMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "coupledFvsPatchFields.H"
#include "correctCentralACMIInterpolation.H"
#include "localEulerDdtScheme.H"

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(standardAndSonicInverseMachKappaFunction, 0);
    addToRunTimeSelectionTable
    (
        kappaFunction,
        standardAndSonicInverseMachKappaFunction,
        dictionary
    );
}
}

namespace Foam
{
namespace fv
{

standardAndSonicInverseMachKappaFunction::standardAndSonicInverseMachKappaFunction
(
    const word& name,
    const word& type,
    const dictionary& parentDict,
    const fvMesh& mesh
)
:
    kappaFunction (name, type, parentDict, mesh)
{
    this->read(this->coeffs_);
}

standardAndSonicInverseMachKappaFunction::~standardAndSonicInverseMachKappaFunction()
{
}

void standardAndSonicInverseMachKappaFunction::update()
{
}

tmp<surfaceScalarField> standardAndSonicInverseMachKappaFunction::kappa()
{
    const surfaceScalarField& amaxSf  = mesh_.thisDb().lookupObject<surfaceScalarField>("amaxSf");
    const surfaceScalarField& phi     = mesh_.thisDb().lookupObject<surfaceScalarField>("phi");
    const surfaceScalarField& rho_own = mesh_.thisDb().lookupObject<surfaceScalarField>("rho_own");
    const surfaceScalarField& rho_nei = mesh_.thisDb().lookupObject<surfaceScalarField>("rho_nei");
    const surfaceScalarField& a_own   = mesh_.thisDb().lookupObject<surfaceScalarField>("agrp_own");
    const surfaceScalarField& a_nei   = mesh_.thisDb().lookupObject<surfaceScalarField>("agrp_nei");
    // const surfaceScalarField& a_own   = mesh_.thisDb().lookupObject<surfaceScalarField>("alpha_own");
    // const surfaceScalarField& a_nei   = mesh_.thisDb().lookupObject<surfaceScalarField>("alpha_nei");
    const surfaceScalarField& cf_own  = mesh_.thisDb().lookupObject<surfaceScalarField>("cf_own");
    const surfaceScalarField& cf_nei  = mesh_.thisDb().lookupObject<surfaceScalarField>("cf_nei");
    const surfaceScalarField& uMagSf  = mesh_.thisDb().lookupObject<surfaceScalarField>("uMagSf");
	const surfaceScalarField& YbarLiqf   = mesh_.thisDb().lookupObject<surfaceScalarField>("YbarLiqf");

    dimensionedScalar v_eps_("v_eps", dimVelocity, 1e-15);
    dimensionedScalar rho_eps_("rho_eps", dimMass/dimVolume, 1e-15);

    surfaceScalarField amaxSfbyDelta
    (
        mesh_.surfaceInterpolation::deltaCoeffs()*amaxSf
    );

    autoPtr<surfaceScalarField> FaceAcCourantPtr;

    if (fv::localEulerDdt::enabled(mesh_))
    {
        dimensionedScalar cDeltaT = runTime_.deltaT();

        cDeltaT.value() = 1.0 / gMax
            (
                mesh_.thisDb().lookupObject<volScalarField>(fv::localEulerDdt::rDeltaTName).internalField()
            );

        FaceAcCourantPtr.reset
        (
            new surfaceScalarField
            (
                "FaceAcCourant",
                amaxSfbyDelta/uMagSf * cDeltaT
            )
        );
    }
    else
    {
        FaceAcCourantPtr.reset
        (
            new surfaceScalarField
            (
                "FaceAcCourant",
                amaxSfbyDelta/uMagSf * runTime_.deltaT()
            )
        );
    }

    const surfaceScalarField& FaceAcCourant = FaceAcCourantPtr();

    Info << "max/min FaceAcCourant: " << max(FaceAcCourant).value() << "/" << min(FaceAcCourant).value() << endl;
    surfaceScalarField rhof = max(rho_own*a_own + rho_nei*a_nei, rho_eps_);
    surfaceScalarField cf = max(cf_own*a_own + cf_nei*a_nei,v_eps_);
    // Info << "max/min a_own: " << max(a_own).value() << "/" << min(a_own).value() << endl;
    // Info << "max/min a_nei: " << max(a_nei).value() << "/" << min(a_nei).value() << endl;
    // Info << "max/min rhof: " << max(rhof).value() << "/" << min(rhof).value() << endl;
    // Info << "max/min cf: " << max(cf).value() << "/" << min(cf).value() << endl;
    // Info << "max/min UnSf: " << max(uMagSf).value() << "/" << min(uMagSf).value() << endl;

    surfaceScalarField Maf
    (
        "Maf",
        max
        (
            (mag(phi) / rhof / uMagSf)
            /cf,
            scalar(0)
            //(mag(phi) / (rho_own*a_own + rho_nei*a_nei) / uMagSf)
            /// (cf_own*a_own + cf_nei*a_nei),
            //scalar(0)
        )
    );

    Info << "max/min Maf: " << max(Maf).value() << "/" << min(Maf).value() << endl;

	surfaceScalarField cfbyDelta
    (
        mesh_.surfaceInterpolation::deltaCoeffs()
        *
        (
            cf_own*a_own + cf_nei*a_nei
        )
    );

    dimensionedScalar cDeltaT = runTime_.deltaT();

    if (fv::localEulerDdt::enabled(mesh_))
    {
        cDeltaT.value() = 1.0 / gMax
        (
            mesh_.thisDb().lookupObject<volScalarField>(fv::localEulerDdt::rDeltaTName)
        );
    }

    surfaceScalarField FaceSonicCourant
    (
        "FaceSonicCourant",
        (cfbyDelta * cDeltaT)
    );

    Info << "max/min FaceSonicCourant: " << max(FaceSonicCourant).value() << "/" << min(FaceSonicCourant).value() << endl;
	
	//Maf.setOriented(true);
    FaceSonicCourant.setOriented(false);
    tmp<surfaceScalarField> tKappa
    (
        min
        (
            scalar(0.2) / FaceSonicCourant * pos(YbarLiqf-0.0001)
            +
            Maf / FaceAcCourant * pos(Maf / FaceAcCourant - scalar(0.2) / FaceSonicCourant),
        //  Maf / FaceAcCourant
		//  + 
		//  pos(0.05 / FaceSonicCourant - scalar(1)) * pos0((1-YbarLiqf) * (YbarLiqf)),
            scalar(1.0)
        )
    );
	
    surfaceScalarField& kappa = tKappa.ref();
	kappa.setOriented(true);
    //resetCoupledBoundaries(kappa);

    writeMaxMinKappa(kappa);

    return tKappa;
}

void standardAndSonicInverseMachKappaFunction::writeData (Ostream& os) const
{
    kappaFunction::writeData(os);
}

bool standardAndSonicInverseMachKappaFunction::read(const dictionary& dict)
{
    if (kappaFunction::read(dict))
    {
        return true;
    }
    else
    {
        return false;
    }

    return true;
}

}; //namespace fv

}; //namespace Foam


//END-OF-FILE
