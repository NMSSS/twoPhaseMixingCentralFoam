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

#include "compressibleTwoPhaseMixtureThermo.H"
#include "fixedEnergyFvPatchScalarField.H"
#include "gradientEnergyFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(compressibleTwoPhaseMixtureThermo, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::compressibleTwoPhaseMixtureThermo::heBoundaryCorrection(volScalarField& h)
{
    volScalarField::Boundary& hbf = h.boundaryFieldRef();

    forAll(hbf, patchi)
    {
        if (isA<gradientEnergyFvPatchScalarField>(hbf[patchi]))
        {
            refCast<gradientEnergyFvPatchScalarField>(hbf[patchi]).gradient()
                = hbf[patchi].fvPatchField::snGrad();
        }
        else if (isA<mixedEnergyFvPatchScalarField>(hbf[patchi]))
        {
            refCast<mixedEnergyFvPatchScalarField>(hbf[patchi]).refGrad()
                = hbf[patchi].fvPatchField::snGrad();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressibleTwoPhaseMixtureThermo::compressibleTwoPhaseMixtureThermo
(
    const fvMesh& mesh
)
:
    rhoThermo(mesh, word::null),
    compressibleTwoPhaseMixture(mesh, *this),
    /*pMin_("pMin",dimensionSet(1,-1,-2,0,0),0.0,  *this),*/
    pMin_("pMin", dimPressure, static_cast<const rhoThermo&>(*this)),
    thermoLiq_(nullptr),
    thermoGas_(nullptr),
    he_
    (
        IOobject
        (
            "h",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass,
        this->heBoundaryTypes(),
        this->heBoundaryBaseTypes()
    ),
    
    saturationModelPtr_
    (
        saturationModel::New
        (
            static_cast<const rhoThermo&>(*this).subDict("saturationModel"),
            mesh
        )
    ),
    Pc_("Pc", dimPressure, static_cast<const rhoThermo&>(*this)),
    Tc_("Tc", dimTemperature, static_cast<const rhoThermo&>(*this)),
    Pt_("Pt", dimPressure, static_cast<const rhoThermo&>(*this)),
    Tt_("Tt", dimTemperature, static_cast<const rhoThermo&>(*this))
    
{
    //limit pressure before proceeding
    p_ = max(p_, pMin_);

    thermoLiq_ .reset
    (
        new simplePhasePsiThermo(mesh, this->subDict(liqPhaseName()))
    );

    thermoGas_.reset
    (
        new simplePhasePsiThermo(mesh, this->subDict(gasPhaseName()))
    );

    he_ = YLiq()*thermoLiq_->he() + YGas()*thermoGas_->he();

    heBoundaryCorrection(he_);

    correct();

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::compressibleTwoPhaseMixtureThermo::~compressibleTwoPhaseMixtureThermo()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

static void findTemperature
(
    Foam::scalarField& T,
    const Foam::scalarField& p,
    const Foam::scalarField& YLiq,
    const Foam::scalarField& YGas,
    const Foam::scalarField& href,
    const Foam::simplePhasePsiThermo& liq,
    const Foam::simplePhasePsiThermo& gas
)
{
    static const Foam::scalar EpsilonH = 1.0e-5;
    Foam::scalar F     = 0.0;
    Foam::scalar dFdT  = 0.0;
    Foam::scalar Cpm   = 0.0;
    Foam::scalar deltaH= 0.0;

    forAll(T, elemi)
    {
        Cpm =  liq.Cp_inline(p[elemi], T[elemi]) * YLiq[elemi] +
            gas.Cp_inline(p[elemi], T[elemi]) * YGas[elemi];

        deltaH = YLiq[elemi]*liq.deltah(p[elemi], T[elemi], href[elemi])
            + YGas[elemi]*gas.deltah(p[elemi], T[elemi], href[elemi]);
        F = deltaH;
        dFdT = Cpm;

        while (Foam::mag(deltaH) >= EpsilonH)
        {
            T[elemi] = T[elemi] - F / dFdT;

            Cpm =  liq.Cp_inline(p[elemi], T[elemi]) * YLiq[elemi] + gas.Cp_inline(p[elemi], T[elemi]) * YGas[elemi];
            deltaH = YLiq[elemi]*liq.deltah(p[elemi], T[elemi], href[elemi])
            + YGas[elemi]*gas.deltah(p[elemi], T[elemi], href[elemi]);
            F = deltaH;
            dFdT = Cpm;
        }
    }
}

void Foam::compressibleTwoPhaseMixtureThermo::correct()
{
    //update temperature and heat capacities of liquid and gas

    findTemperature
    (
        T_.primitiveFieldRef(),
        p_.primitiveField(),
        YLiq_.primitiveField(),
        YGas_.primitiveField(),
        he_.primitiveField(),
        thermoLiq_(),
        thermoGas_()
    );

    forAll(T_.boundaryField(), patchi)
    {
        fvPatchScalarField&       pT  = T_.boundaryFieldRef()[patchi];
        const fvPatchScalarField& pp  = p_.boundaryField()[patchi];
        fvPatchScalarField&       phe = he_.boundaryFieldRef()[patchi];
        const fvPatchScalarField& pyl = YLiq_.boundaryField()[patchi];
        const fvPatchScalarField& pyg = YGas_.boundaryField()[patchi];

        if (pT.fixesValue())
        {
            phe == YLiq().boundaryField()[patchi] * thermoLiq_->he(pp, pT, patchi) +
                YGas().boundaryField()[patchi] * thermoGas_->he(pp, pT, patchi);
        }
        else
        {
            findTemperature
            (
                pT,
                pp,
                pyl,
                pyg,
                phe,
                thermoLiq_(),
                thermoGas_()
            );
        }
    }

    Info << "max/min T: " << max(T_).value() << "/" << min(T_).value() << endl;

    //correct properties of liquid
    thermoLiq_->correct();

    //correct properties of gas
    thermoGas_->correct();

    //correct mixture properties
    rhoEff() = 1.0 / (YLiq()/thermoLiq_->rho() + YGas()/thermoGas_->rho());

    updateVolFrac(thermoLiq_->rho(), thermoGas_->rho());

    const volScalarField& rhoGas = thermoGas_->rho();
    const volScalarField& rhoLiq = thermoLiq_->rho();
    const volScalarField& psiGas = thermoGas_->psi();
    const volScalarField& psiLiq = thermoLiq_->psi();

    volScalarField YLiqByRhoLiq (YLiq() / rhoLiq);
    volScalarField YGasByRhoGas (YGas() / rhoGas);
    volScalarField OneByRhoSqr(Foam::pow(YLiqByRhoLiq + YGasByRhoGas,-2.0));

    psiLiqEff_ = psiLiq / (rhoLiq * rhoLiq);
    psiGasEff_ = psiGas / (rhoGas * rhoGas);
    psi_ = (YLiq()*psiLiqEff_ + YGas()*psiGasEff_)*OneByRhoSqr;

    /*
    linear
    */
    // psiLiqEff_ = psiLiq / (rhoLiq * rhoLiq);
    // psiGasEff_ = psiGas / (rhoGas * rhoGas);
    // psi_ = (YbarLiq()*psiLiq + YbarGas()*psiGas);


    mu_ = YbarLiq()*thermoLiq_->mu() + YbarGas()*thermoGas_->mu();
    alpha_ = YbarLiq()*thermoLiq_->alpha() + YbarGas()*thermoGas_->alpha();
	
	//lambda_= YbarLiq()*thermoLiq_->lambda() + YbarGas()*thermoGas_->lambda();
	
    //mu_ = YLiq()*thermoLiq_->mu() + YGas()*thermoGas_->mu();
    //alpha_ = YLiq()*thermoLiq_->alpha() + YGas()*thermoGas_->alpha();
}


Foam::word Foam::compressibleTwoPhaseMixtureThermo::thermoName() const
{
    return "compressibleTwoPhaseMixtureThermo";
}

bool Foam::compressibleTwoPhaseMixtureThermo::incompressible() const
{
    return false;
}


bool Foam::compressibleTwoPhaseMixtureThermo::isochoric() const
{
    return false;
}


Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixtureThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return YLiq()*thermoLiq_->he(p, T) + YGas()*thermoGas_->he(p, T);
}


Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    return
        scalarField(YLiq(), cells)*thermoLiq_->he(p, T, cells)
      + scalarField(YGas(), cells)*thermoGas_->he(p, T, cells);
}


Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        YLiq().boundaryField()[patchi]*thermoLiq_->he(p, T, patchi)
      + YGas().boundaryField()[patchi]*thermoGas_->he(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixtureThermo::hc() const
{
    return YLiq()*thermoLiq_->hc() + YGas()*thermoGas_->hc();
}


Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const labelList& cells
) const
{
    notImplemented("compressibleTwoPhaseMixtureThermo::THE(...)");
    return T0;
}


Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const label patchi
) const
{
    notImplemented("compressibleTwoPhaseMixtureThermo::THE(...)");
    return T0;
}


Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixtureThermo::Cp() const
{
    return YLiq()*thermoLiq_->Cp() + YGas()*thermoGas_->Cp();
}


Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::Cp
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        YLiq().boundaryField()[patchi]*thermoLiq_->Cp(p, T, patchi)
      + YGas().boundaryField()[patchi]*thermoGas_->Cp(p, T, patchi);
}

Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::Cp
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    volScalarField CpAll(this->Cp());
    tmp<scalarField> tCp( new scalarField(cells.size()));
    scalarField& Cp = tCp.ref();
    forAll(cells,icell)
    {
        Cp[icell] = CpAll.primitiveField()[cells[icell]];
    }
    
    return tCp;
}


Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixtureThermo::Cv() const
{
    return YLiq()*thermoLiq_->Cv() + YGas()*thermoGas_->Cv();
}


Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::Cv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        YLiq().boundaryField()[patchi]*thermoLiq_->Cv(p, T, patchi)
      + YGas().boundaryField()[patchi]*thermoGas_->Cv(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixtureThermo::gamma() const
{
    return Cp() / Cv();
}


Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::gamma
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return Cp(p, T, patchi) / Cv(p, T, patchi);
//        YLiq().boundaryField()[patchi]*thermoLiq_->gamma(p, T, patchi)
//      + YGas().boundaryField()[patchi]*thermoGas_->gamma(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixtureThermo::Cpv() const
{
    return YLiq()*thermoLiq_->Cpv() + YGas()*thermoGas_->Cpv();
}


Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::Cpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        YLiq().boundaryField()[patchi]*thermoLiq_->Cpv(p, T, patchi)
      + YGas().boundaryField()[patchi]*thermoGas_->Cpv(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixtureThermo::CpByCpv() const
{
    return
        YLiq()*thermoLiq_->CpByCpv()
      + YGas()*thermoGas_->CpByCpv();
}


Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::CpByCpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        YLiq().boundaryField()[patchi]*thermoLiq_->CpByCpv(p, T, patchi)
      + YGas().boundaryField()[patchi]*thermoGas_->CpByCpv(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixtureThermo::kappa() const
{
    return YLiq()*thermoLiq_->kappa() + YGas()*thermoGas_->kappa();
}


Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::kappa
(
    const label patchi
) const
{
    return
        YLiq().boundaryField()[patchi]*thermoLiq_->kappa(patchi)
      + YGas().boundaryField()[patchi]*thermoGas_->kappa(patchi);
}


Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixtureThermo::kappaEff
(
    const volScalarField& alphat
) const
{
    return
        YLiq()*thermoLiq_->kappaEff(alphat)
      + YGas()*thermoGas_->kappaEff(alphat);
}


Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        YLiq().boundaryField()[patchi]*thermoLiq_->kappaEff(alphat, patchi)
      + YGas().boundaryField()[patchi]*thermoGas_->kappaEff(alphat, patchi)
    ;
}


Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixtureThermo::alphaEff
(
    const volScalarField& alphat
) const
{
    return
        YLiq()*thermoLiq_->alphaEff(alphat)
      + YGas()*thermoGas_->alphaEff(alphat);
}


Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        YLiq().boundaryField()[patchi]*thermoLiq_->alphaEff(alphat, patchi)
      + YGas().boundaryField()[patchi]*thermoGas_->alphaEff(alphat, patchi)
    ;
}


Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixtureThermo::mu() const
{
    return tmp<volScalarField>(this->mu_);
}

Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::mu(const label patchi) const
{
    return mu_.boundaryField()[patchi];
}

Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixtureThermo::alphahe() const
{
    return YLiq()*thermoLiq_->alpha() + YGas()*thermoGas_->alpha();
}

Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::alphahe(const label patchi) const
{
    return
        YLiq().boundaryField()[patchi]*thermoLiq_->alpha(patchi)
      + YGas().boundaryField()[patchi]*thermoGas_->alpha(patchi);
}

Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixtureThermo::rho() const
{
    return rhoEff_;
}

Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::rho(const label patchi) const
{
    return rhoEff_.boundaryField()[patchi];
}

Foam::tmp<Foam::scalarField>
Foam::compressibleTwoPhaseMixtureThermo::rhoEoS
(
    const scalarField& p, 
    const scalarField& T,
    const labelList& cells
) const
{
    volScalarField rho (this->rho());
    tmp<scalarField> tRho(new scalarField(cells.size()));
    scalarField& rhoEoS = tRho.ref();
    forAll(cells,icell)
    {
        rhoEoS[icell] = rho.primitiveField()[cells[icell]];
    }
    return tRho;
}

const Foam::dimensionedScalar& Foam::compressibleTwoPhaseMixtureThermo::pMin() const
{
    return pMin_;
}

Foam::volScalarField& Foam::compressibleTwoPhaseMixtureThermo::h()
{
    return he_;
}

const Foam::volScalarField& Foam::compressibleTwoPhaseMixtureThermo::h() const
{
    return he_;
}

Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::h
(
    const scalarField& T,
    const label patchi
) const
{
    return tmp<scalarField>(he_.boundaryField()[patchi]);
}

void Foam::compressibleTwoPhaseMixtureThermo::correctRealDensities()
{
    volScalarField pLimited
    (
        "pLim",
        max(pMin_, p_)
    );
    thermoLiq_->correctDensity(pLimited);
    thermoGas_->correctDensity(pLimited);
}

Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixtureThermo::Tsat
(
    const volScalarField& p
) const
{
    return saturationModelPtr_->Tsat(p);
}


Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixtureThermo::pSat
(
    const volScalarField& T
) const
{
    return saturationModelPtr_->pSat(T);
}

Foam::tmp<Foam::surfaceScalarField> 
Foam::compressibleTwoPhaseMixtureThermo::Cfm
(
    const surfaceScalarField& pf,
    const surfaceScalarField& Tf,
    const surfaceScalarField& YbarLiqf,
    const surfaceScalarField& rhoLiqf,
    const surfaceScalarField& rhoGasf
) const
{

    tmp<surfaceScalarField> tCfm
    (
       YbarLiqf*dimensionedScalar("Cfm0", dimLength / dimTime, 0.0)
    );
    surfaceScalarField &CfmRef = tCfm.ref();
    CfmRef.rename("Cfm");

    
    const simplePhasePsiThermo& Liq = this->thermoLiq_;
    const simplePhasePsiThermo& Gas = this->thermoGas_;
    auto CsMixture = [&Liq, &Gas](const scalar& p, const scalar& T, const scalar& ybarliq, const scalar& rholiq, const scalar& rhogas)
    {
        scalar psiliq = 0.0;
        scalar psigas = 0.0;
        scalar psiliq_eff = 0.0;
        scalar psigas_eff = 0.0;
        scalar rhom = 0.0;
        scalar yliq = 0.0;
        scalar ygas = 0.0;
        scalar psim = 0.0;
        scalar onebyrhosqr = 0.0;
        scalar Cpm = 0.0;
        scalar Cvm = 0.0;
        scalar gammam = 0.0;

        psiliq = Liq.psi(p, T);
        psigas = Gas.psi(p, T);
        rhom = rholiq*ybarliq + rhogas*(1.0 -  ybarliq);

        yliq = rholiq*ybarliq / rhom;
        ygas = 1.0 - yliq;
        psiliq_eff = psiliq / sqr(rholiq);
        psigas_eff = psigas / sqr(rhogas);
        onebyrhosqr = 1.0 / sqr((yliq / rholiq) + (ygas / rhogas));

        Cpm  = yliq*Liq.Cp_inline(p, T)+
               (1.0 - yliq)*Gas.Cp_inline(p, T);
        Cvm  = yliq*Liq.Cv(p, T)+
               (1.0 - yliq)*Gas.Cv(p, T);
        gammam = Cpm / Cvm;
        psim = (yliq*psiliq_eff + ygas*psigas_eff)*onebyrhosqr;
        //psim = ybarliq*psiliq + (1.0 - ybarliq)*psigas; //linear model
        return sqrt(gammam / psim);
    };

    
    forAll(Tf.primitiveField(), iface)
    {
        CfmRef.primitiveFieldRef()[iface] = CsMixture
        (
            pf.primitiveField()[iface],
            Tf.primitiveField()[iface],
            YbarLiqf.primitiveField()[iface],
            rhoLiqf.primitiveField()[iface],
            rhoGasf.primitiveField()[iface]
        );
    }

    forAll(Tf.boundaryField(), ipatch)
    {
        scalarField& Cfp = CfmRef.boundaryFieldRef()[ipatch];
        
        const scalarField& pfp = pf.boundaryField()[ipatch];
        const scalarField& Tfp = Tf.boundaryField()[ipatch];
        const scalarField& YbarLiqfp = YbarLiqf.boundaryField()[ipatch];
        const scalarField& rhoLiqfp = rhoLiqf.boundaryField()[ipatch];
        const scalarField& rhoGasfp = rhoGasf.boundaryField()[ipatch];

        forAll(Cfp, iface)
        {
            Cfp[iface] = CsMixture
            (
                pfp[iface],
                Tfp[iface],
                YbarLiqfp[iface],
                rhoLiqfp[iface],
                rhoGasfp[iface]
            );
        }
    }

    return tCfm;

}


// ************************************************************************* //
