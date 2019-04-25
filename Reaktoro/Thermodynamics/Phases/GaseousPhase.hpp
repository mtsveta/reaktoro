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

#pragma once

// C++ includes
#include <memory>

// Reaktoro includes
#include <Reaktoro/Core/Phase.hpp>

namespace Reaktoro {

// Forward declarations
class GaseousMixture;

/// Class that defines a gaseous phase
class GaseousPhase : public Phase
{
public:
    /// Construct a default GaseousPhase instance.
    GaseousPhase();

    /// Construct an GaseousPhase instance with given gaseous mixture.
    /// The Peng-Robinson equation of state is chosen by default to calculate the
    /// thermodynamic and chemical properties of this GaseousPhase object.
    explicit GaseousPhase(const GaseousMixture& mixture);

    /// Set the chemical model of the phase with the ideal gas equation of state.
    auto setChemicalModelIdeal() -> GaseousPhase&;

    /// Set the chemical model of the phase with the van der Waals equation of state.
    ///
    /// Reference: *van der Waals, J.D. (1910). The equation of state for gases and liquids. Nobel Lectures in Physics. pp. 254-265*.
    auto setChemicalModelVanDerWaals() -> GaseousPhase&;

    /// Set the chemical model of the phase with the Redlich-Kwong equation of state.
    ///
    /// Reference: *Redlich, O., Kwong, J.N.S. (1949). On The Thermodynamics of Solutions. Chem. Rev. 44(1) 233–244*.
    auto setChemicalModelRedlichKwong() -> GaseousPhase&;

    /// Set the chemical model of the phase with the Soave-Redlich-Kwong equation of state.
    ///
    /// Reference: *Soave, G. (1972). Equilibrium constants from a modified Redlich-Kwong equation of state, Chem. Eng. Sci., 27, 1197-1203*.
    auto setChemicalModelSoaveRedlichKwong() -> GaseousPhase&;

    /// Set the chemical model of the phase with the Peng-Robinson equation of state.
    ///
    /// Reference: *Peng, D.Y., Robinson, D.B. (1976). A New Two-Constant Equation of State. Industrial and Engineering Chemistry: Fundamentals 15: 59–64*.
    auto setChemicalModelPengRobinson() -> GaseousPhase&;

    /// Set the chemical model of the phase with the Spycher et al. (2003) equation of state.
    /// This model only supports the gaseous species `H2O(g)` and `CO2(g)`. Any other species
    /// will result in a runtime error.
    ///
    /// Reference: *Spycher, N., Pruess, K., Ennis-King, J. (2003). CO2-H2O mixtures in the
    /// geological sequestration of CO2. I. Assessment and calculation of mutual solubilities from 12 to 100°C
    /// and up to 600 bar. Geochimica et Cosmochimica Acta, 67(16), 3015–3031*.
    auto setChemicalModelSpycherPruessEnnis() -> GaseousPhase&;

    /// Set the chemical model of the phase with the Spycher and Reed (1988) equation of state.
    /// This model only supports the gaseous species `H2O(g)`, `CO2(g)`, and CH4(g). Any other
    /// species will result in a runtime error.
    ///
    /// Reference: *Spycher, N., Reed, M. (1988). Fugacity coefficients of H2, CO2,
    /// CH4, H2O and of H2O--CO2--CH4 mixtures: A virial equation treatment for
    /// moderate pressures and temperatures applicable to calculations of
    /// hydrothermal boiling. Geochimica et Cosmochimica Acta, 52(3), 739–749*.
    auto setChemicalModelSpycherReed() -> GaseousPhase&;

    /// Return the GaseousMixture instance
    auto mixture() const -> const GaseousMixture&;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktoro
