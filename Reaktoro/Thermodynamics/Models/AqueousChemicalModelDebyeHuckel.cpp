// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include "AqueousChemicalModelDebyeHuckel.hpp"

// C++ includes
#include <map>
#include <string>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/ConvertUtils.hpp>
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Math/BilinearInterpolator.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/AqueousMixture.hpp>
#include <Reaktoro/Thermodynamics/Species/AqueousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>

namespace Reaktoro {

auto aqueousChemicalModelDebyeHuckel(const AqueousMixture& mixture, const DebyeHuckelParams& params) -> PhaseChemicalModel
{
    // The natural log of 10
    const double ln10 = std::log(10);

    // The molar mass of water
    const double Mw = waterMolarMass;

    // The number of moles of water per kg
    const double nwo = 1/Mw;

    // The number of species in the mixture
    const Index num_species = mixture.numSpecies();

    // The number of charged and neutral species in the mixture
    const Index num_charged_species = mixture.numChargedSpecies();
    const Index num_neutral_species = mixture.numNeutralSpecies();

    // The indices of the charged and neutral species
    const Indices icharged_species = mixture.indicesChargedSpecies();
    const Indices ineutral_species = mixture.indicesNeutralSpecies();

    // The index of the water species
    const Index iwater = mixture.indexWater();

    // The electrical charges of the charged species only
    const Vector charges = mixture.chargesChargedSpecies();

    // The Debye-Huckel parameters a and b of the charged species
    std::vector<double> aions, bions;

    // The Debye-Huckel parameter b of the neutral species
    std::vector<double> bneutral;

    // Collect the Debye-Huckel parameters a and b of the charged species
    for(Index i : icharged_species)
    {
        const AqueousSpecies& species = mixture.species(i);
        aions.push_back(params.aion(species.name()));
        bions.push_back(params.bion(species.name()));
    }

    // Collect the Debye-Huckel parameter b of the neutral species
    for(Index i : ineutral_species)
    {
        const AqueousSpecies& species = mixture.species(i);
        bneutral.push_back(params.bneutral(species.name()));
    }

    // The state of the aqueous mixture
    AqueousMixtureState state;

    // Auxiliary variables
    ChemicalScalar xw, ln_xw, I2, sqrtI, mSigma, sigma(num_species), sigmacoeff, Lambda;
    ChemicalVector ln_m;
    ThermoScalar A, B, sqrt_rho, T_epsilon, sqrt_T_epsilon;

    // Define the intermediate chemical model function of the aqueous mixture
    PhaseChemicalModel model = [=](PhaseChemicalModelResult& res, Temperature T, Pressure P, VectorConstRef n) mutable
    {
        // Evaluate the state of the aqueous mixture
        state = mixture.state(T, P, n);

        // Auxiliary constant references
        const auto& I = state.Ie;            // ionic strength
        const auto& x = state.x;             // mole fractions of the species
        const auto& m = state.m;             // molalities of the species
        const auto& rho = state.rho/1000;    // density in units of g/cm3
        const auto& epsilon = state.epsilon; // dielectric constant

        // Auxiliary references
        auto& ln_g = res.ln_activity_coefficients;
        auto& ln_a = res.ln_activities;

        // Update auxiliary variables
		ln_m = log(m);
		xw = x[iwater];
		ln_xw = log(xw);
		mSigma = nwo * (1 - xw)/xw;
		I2 = I*I;
		sqrtI = sqrt(I);
		sqrt_rho = sqrt(rho);
		T_epsilon = T * epsilon;
		sqrt_T_epsilon = sqrt(T_epsilon);
		A = 1.824829238e+6 * sqrt_rho/(T_epsilon*sqrt_T_epsilon);
		B = 50.29158649 * sqrt_rho/sqrt_T_epsilon;
		sigmacoeff = (2.0/3.0)*A*I*sqrtI;

        // Loop over all neutral species in the mixture
        for(Index i = 0; i < num_neutral_species; ++i)
        {
            // The index of the current neutral species
            const Index ispecies = ineutral_species[i];

            // Calculate the ln activity coefficient of the current neutral species
            ln_g[ispecies] = ln10 * bneutral[i] * I;

            // Calculate the ln activity coefficient of the current neutral species
            ln_a[ispecies] = ln_g[ispecies] + ln_m[ispecies];
        }

        // Set the first contribution to the activity of water
        ln_a[iwater] = mSigma;

        // Loop over all charged species in the mixture
        for(Index i = 0; i < num_charged_species; ++i)
        {
            // The index of the current charged species
            const Index ispecies = icharged_species[i];

            // The molality of the charged species and its molar derivatives
            const auto mi = m[ispecies];

            // The electrical charge of the charged species
            const auto z = charges[i];

            // Update the Lambda parameter of the Debye-Huckel activity coefficient model
            Lambda = 1.0 + aions[i]*B*sqrtI;

			// Update the sigma parameter of the current ion
            if(aions[i] != 0.0) sigma = 3.0*pow(Lambda - 1, -3) * ((Lambda - 1)*(Lambda - 3) + 2*log(Lambda));
            else                sigma = 2.0;

            // Calculate the ln activity coefficient of the current charged species
            ln_g[ispecies] = ln10 * (-A*z*z*sqrtI/Lambda + bions[i]*I);

            // Calculate the ln activity of the current charged species
            ln_a[ispecies] = ln_g[ispecies] + ln_m[ispecies];

            // Calculate the contribution of current ion to the ln activity of water
			ln_a[iwater] += mi*ln_g[ispecies] + sigmacoeff*sigma*ln10 - I2*bions[i]/(z*z)*ln10;
        }

        // Finalize the computation of the activity of water (in mole fraction scale)
        ln_a[iwater] *= -1.0/nwo;

        // Set the activity coefficient of water (mole fraction scale)
        ln_g[iwater] = ln_a[iwater] - ln_xw;
    };

    return model;
}

struct DebyeHuckelParams::Impl
{
    /// The default value of the `a` parameter for ionic species.
    double aiondefault = 0.0;

    /// The default value of the `b` parameter for ionic species.
    double biondefault = 0.0;

    /// The default value of the `b` parameter for neutral species.
    double bneutraldefault = 0.1;

    /// The parameters `a` of the ionic species.
    std::map<std::string, double> aion;

    /// The parameters `b` of the ionic species.
    std::map<std::string, double> bion;

    /// The parameters `b` of the neutral species.
    std::map<std::string, double> bneutral;

    /// The Debye--Hückel parameter `å` used in PHREEQC v3 (Parkhurst and Appelo, 2013)
    const std::map<std::string, double> aion_phreeqc = {{"Al(OH)2+", 5.4}, {"Al(OH)4-", 4.5}, {"Al(SO4)2-", 4.5}, {"Al+++", 9}, {"AlF++", 5.4}, {"AlF2+", 5.4}, {"AlF4-", 4.5}, {"AlOH++", 5.4}, {"AlSO4+", 4.5}, {"Ba++", 4}, {"BaOH+", 5}, {"Br-", 3}, {"CO3--", 5.4}, {"Ca++", 5}, {"CaH2PO4+", 5.4}, {"CaHCO3+", 6}, {"CaPO4-", 5.4}, {"Cl-", 3.63}, {"Cu+", 2.5}, {"Cu++", 6}, {"CuCl+", 4}, {"CuCl2-", 4}, {"CuCl3-", 4}, {"CuCl3--", 5}, {"CuCl4--", 5}, {"CuOH+", 4}, {"F-", 3.5}, {"Fe(OH)2+", 5.4}, {"Fe(OH)3-", 5}, {"Fe(OH)4-", 5.4}, {"Fe++", 6}, {"Fe+++", 9}, {"FeCl++", 5}, {"FeCl2+", 5}, {"FeF++", 5}, {"FeF2+", 5}, {"FeH2PO4+", 5.4}, {"FeH2PO4++", 5.4}, {"FeHPO4+", 5}, {"FeOH+", 5}, {"FeOH++", 5}, {"FeSO4+", 5}, {"H+", 9}, {"H2PO4-", 5.4}, {"H2SiO4--", 5.4}, {"H3SiO4-", 4}, {"HCO3-", 5.4}, {"HPO4--", 5}, {"HS-", 3.5}, {"K+", 3.5}, {"KHPO4-", 5.4}, {"KSO4-", 5.4}, {"Li+", 6}, {"LiSO4-", 5}, {"Mg++", 5.5}, {"MgF+", 4.5}, {"MgH2PO4+", 5.4}, {"MgHCO3+", 4}, {"MgOH+", 6.5}, {"MgPO4-", 5.4}, {"Mn(OH)3-", 5}, {"Mn++", 6}, {"Mn+++", 9}, {"MnCl+", 5}, {"MnCl3-", 5}, {"MnF+", 5}, {"MnHCO3+", 5}, {"MnOH+", 5}, {"NH4+", 2.5}, {"NO2-", 3}, {"NO3-", 3}, {"Na+", 4.08}, {"NaHPO4-", 5.4}, {"NaSO4-", 5.4}, {"OH-", 3.5}, {"PO4---", 4}, {"S--", 5}, {"SO4--", 5}, {"SiF6--", 5}, {"Sr++", 5.26}, {"SrHCO3+", 5.4}, {"SrOH+", 5}, {"Zn++", 5}, {"ZnCl+", 4}, {"ZnCl3-", 4}, {"ZnCl4--", 5}};

    /// The Debye--Hückel parameter `b` used in PHREEQC v3 (Parkhurst and Appelo, 2013)
    const std::map<std::string, double> bion_phreeqc = {{"Ba++", 0.153}, {"Ca++", 0.165}, {"Cl-", 0.017}, {"K+", 0.015}, {"Mg++", 0.2}, {"Na+", 0.082}, {"SO4--", -0.04}, {"Sr++", 0.121}};

    /// The Debye--Hückel parameter `å` used in WATEQ4F (Ball and Nordstrom 1991, Truesdell and Jones 1974)
    const std::map<std::string, double> aion_wateq4f = {{"Ca++", 5.0}, {"Mg++", 5.5}, {"Na+", 4.0}, {"K+", 3.5}, {"Cl-", 3.5}, {"SO4--", 5.0}, {"HCO3-", 5.4}, {"CO3--", 5.4}, {"Sr++", 5.26}, {"H+", 9.0}, {"OH-", 3.5}, {"SrHCO3+", 5.4}, {"SrOH+", 5.0}, {"Cu(S4)2---", 23.0}, {"CuS4S5---", 25.0}, {"S2--", 6.5}, {"S3--", 8.0}, {"S4--", 10.0}, {"S5--", 12.0}, {"S6--", 14.0}, {"Ag(S4)2---", 22.0}, {"AgS4S5---", 24.0}, {"Ag(HS)S4--", 15.0}};

    /// The Debye--Hückel parameter `b` used in WATEQ4F (Ball and Nordstrom 1991, Truesdell and Jones 1974)
    const std::map<std::string, double> bion_wateq4f = {{"Ca++", 0.165}, {"Mg++", 0.20}, {"Na+", 0.075}, {"K+", 0.015}, {"Cl-", 0.015}, {"SO4--", -0.04}, {"HCO3-", 0.0}, {"CO3--", 0.0}, {"H2CO3(aq)", 0.0}, {"Sr++", 0.121}};

    /// The Debye--Hückel parameter `å` from Kielland (1937).
    const std::map<std::string, double> aion_kielland = {{"H+" , 9.0}, {"Li+" , 6.0}, {"Rb+" , 2.5}, {"Cs+" , 2.5}, {"NH4+" , 2.5}, {"Tl+" , 2.5}, {"Ag+" , 2.5}, {"K+" , 3.0}, {"Cl-" , 3.0}, {"Br-" , 3.0}, {"I-" , 3.0}, {"CN-" , 3.0}, {"NO2-" , 3.0}, {"NO3-" , 3.0}, {"OH-" , 3.5}, {"F-" , 3.5}, {"NCS-" , 3.5}, {"NCO-" , 3.5}, {"HS-" , 3.5}, {"ClO3-" , 3.5}, {"ClO4-" , 3.5}, {"BrO3-" , 3.5}, {"IO4-" , 3.5}, {"MnO4-" , 3.5}, {"Na+" , 4.0}, {"CdCl+" , 4.0}, {"ClO2-" , 4.0}, {"IO3-" , 4.0}, {"HCO3-" , 4.0}, {"H2PO4-" , 4.0}, {"HSO3-" , 4.0}, {"H2AsO4-" , 4.0}, {"Co(NH3)4(NO2)2+" , 4.0}, {"Hg2++" , 4.0}, {"SO4--" , 4.0}, {"S2O3--" , 4.0}, {"S2O6--" , 4.0}, {"S2O8--" , 4.0}, {"SeO4--" , 4.0}, {"CrO4--" , 4.0}, {"HPO4--" , 4.0}, {"Pb++" , 4.5}, {"CO3--" , 4.5}, {"SO3--" , 4.5}, {"MoO4--" , 4.5}, {"Co(NH3)5Cl++" , 4.5}, {"Fe(CN)5NO--" , 4.5}, {"Sr++" , 5.0}, {"Ba++" , 5.0}, {"Ra++" , 5.0}, {"Cd++" , 5.0}, {"Hg++" , 5.0}, {"S--" , 5.0}, {"S2O4--" , 5.0}, {"WO4--" , 5.0}, {"Ca++" , 6.0}, {"Cu++" , 6.0}, {"Zn++" , 6.0}, {"Sn++" , 6.0}, {"Mn++" , 6.0}, {"Fe++" , 6.0}, {"Ni++" , 6.0}, {"Co++" , 6.0}, {"Mg++" , 8.0}, {"Be++" , 8.0}, {"PO4---" , 4.0}, {"Fe(CN)6---" , 4.0}, {"Cr(NH3)6+++" , 4.0}, {"Co(NH3)6+++" , 4.0}, {"Co(NH3)5H2O+++" , 4.0}, {"Al+++" , 9.0}, {"Fe+++" , 9.0}, {"Cr+++" , 9.0}, {"Sc+++" , 9.0}, {"Y+++" , 9.0}, {"La+++" , 9.0}, {"In+++" , 9.0}, {"Ce+++" , 9.0}, {"Pr+++" , 9.0}, {"Nd+++" , 9.0}, {"Sm+++" , 9.0}, {"Fe(CN)6----" , 5.0}, {"Co(S2O3)(CN)5----" , 6.0}, {"Th++++" , 11.0}, {"Zn++++" , 11.0}, {"Ce++++" , 11.0}, {"Sn++++" , 11.0}, {"Co(SO3)2(CN)4-----" , 9.0}};

    Impl()
    {
    }
};

DebyeHuckelParams::DebyeHuckelParams()
: pimpl(new Impl())
{
	setPHREEQC();
}

auto DebyeHuckelParams::aiondefault(double value) -> void
{
    pimpl->aiondefault = value;
}

auto DebyeHuckelParams::aiondefault() const -> double
{
    return pimpl->aiondefault;
}

auto DebyeHuckelParams::aion(std::string name, double value) -> void
{
    pimpl->aion[name] = value;
}

auto DebyeHuckelParams::aion(const std::map<std::string, double>& pairs) -> void
{
    for(auto pair : pairs)
        aion(pair.first, pair.second);
}

auto DebyeHuckelParams::aion(double value) -> void
{
    for(auto& pair : pimpl->aion)
        pair.second = value;
    aiondefault(value);
}

auto DebyeHuckelParams::aion(std::string name) const -> double
{
    auto it = pimpl->aion.find(name);
    return it != pimpl->aion.end() ? it->second : aiondefault();
}

auto DebyeHuckelParams::biondefault(double value) -> void
{
    pimpl->biondefault = value;
}

auto DebyeHuckelParams::biondefault() const -> double
{
    return pimpl->biondefault;
}

auto DebyeHuckelParams::bion(std::string name, double value) -> void
{
    pimpl->bion[name] = value;
}

auto DebyeHuckelParams::bion(const std::map<std::string, double>& pairs) -> void
{
    for(auto pair : pairs)
        bion(pair.first, pair.second);
}

auto DebyeHuckelParams::bion(double value) -> void
{
    for(auto& pair : pimpl->bion)
        pair.second = value;
    biondefault(value);
}

auto DebyeHuckelParams::bion(std::string name) const -> double
{
    auto it = pimpl->bion.find(name);
    return it != pimpl->bion.end() ? it->second : biondefault();
}

auto DebyeHuckelParams::bneutraldefault(double value) -> void
{
    pimpl->bneutraldefault = value;
}

auto DebyeHuckelParams::bneutraldefault() const -> double
{
    return pimpl->bneutraldefault;
}

auto DebyeHuckelParams::bneutral(std::string name, double value) -> void
{
    pimpl->bneutral[name] = value;
}

auto DebyeHuckelParams::bneutral(const std::map<std::string, double>& pairs) -> void
{
    for(auto pair : pairs)
        bneutral(pair.first, pair.second);
}

auto DebyeHuckelParams::bneutral(double value) -> void
{
    for(auto& pair : pimpl->bneutral)
        pair.second = value;
    bneutraldefault(value);
}

auto DebyeHuckelParams::bneutral(std::string name) const -> double
{
    auto it = pimpl->bneutral.find(name);
    return it != pimpl->bneutral.end() ? it->second : bneutraldefault();
}

auto DebyeHuckelParams::setLimitingLaw() -> void
{
    aion(0.0);
    bion(0.0);
}

auto DebyeHuckelParams::setKielland1937() -> void
{
    aion(pimpl->aion_kielland);
}

auto DebyeHuckelParams::setWATEQ4F() -> void
{
    aion(pimpl->aion_wateq4f);
    bion(pimpl->bion_wateq4f);
}

auto DebyeHuckelParams::setPHREEQC() -> void
{
    aion(pimpl->aion_phreeqc);
    bion(pimpl->bion_phreeqc);
    bneutraldefault(0.1);
}

} // namespace Reaktoro
