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
#include <Reaktoro/Kinetics/KineticOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumOptions.hpp>

namespace Reaktoro {

/// The options for the smart kinetics calculations.
/// @see SmartKineticSolver
struct SmartKineticOptions
{
    /// The boolean flag that indicates whether smart equilibrium solver should be used.
    bool use_smart_equilibrium_solver = false;

    /// The options for the equilibrium solver.
    EquilibriumOptions equilibrium;

    /// The options for the smart equilibrium solver.
    SmartEquilibriumOptions smart_equilibrium;

    /// The kinetic options to be used during learning operations.
    KineticOptions learning;

    /// The relative tolerance for estimated species mole amounts.
    double reltol = 1.0;

    /// The absolute tolerance for estimated species mole amounts.
    double abstol = 1e-14;

    /// The cutoff value for normalized species and element amounts.
    /// This parameter is used to ignore certain species/elements with tiny amounts
    /// during the error test for the equilibrium species.
    /// The species and elements with normalized amounts,
    /// \eq{n_{i}/\text{sum}(n)} and \eq{b_{j}/\text{sum}(b)} respectively,
    /// below this threshold will be discarded from the error control.
    double cutoff = -1e-5;

    /// The cutoff value for species mole fractions.
    /// This parameter is used to ignore certain species with tiny mole fractions
    /// during the error test. Those with mole fractions below this threshold
    /// will be discarded from the error control.
    double mole_fraction_cutoff = 1.0e-6;


};

} // namespace Reaktoro
