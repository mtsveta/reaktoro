// Reaktoro is a unified framework for modeling chemically reactive systems.
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

// Reaktoro includes
#include <Reaktoro/Kinetics/ODEState.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Common/ChemicalVector.hpp>

namespace Reaktoro {

/// The record of the knowledge database containing input, output, and derivatives data.
struct KineticRecord
{
    /// The temperature of the equilibrium state (in units of K).
    double T;

    /// The pressure of the equilibrium state (in units of Pa).
    double P;

    /// The amounts of elements in the equilibrium state (in units of mol).
    Vector be;

    /// The calculated equilibrium state at `T`, `P`, `be`.
    ChemicalState state;

    /// The chemical properties at the calculated equilibrium state.
    ChemicalProperties properties;

    /// The sensitivity derivatives at the calculated equilibrium state.
    EquilibriumSensitivity sensitivity;

    /// The matrix used to compute relative change of chemical potentials due to change in `be`.
    Matrix Mbe;

    /// Calculated kinetic state, containing initial and final value of [be, nk], initial and final time, and sensitivities
    ODEState ode_state;

    // Chemical kinetic rates
    ChemicalVector rates;
};

} // namespace Reaktoro