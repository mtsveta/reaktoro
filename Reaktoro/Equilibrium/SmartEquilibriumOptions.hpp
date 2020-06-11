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

// Reaktoro includes
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>

namespace Reaktoro {

/// The options for the smart equilibrium calculations.
/// @see SmartEquilibriumSolver
struct SmartEquilibriumOptions
{
    /// The equilibrium options to be used during learning operations.
    EquilibriumOptions learning;

    /// The cutoff value for species mole fractions.
    /// This parameter is used to ignore certain species with tiny mole fractions
    /// during the error test. Those with mole fractions below this threshold
    /// will be discarded from the error control.
    double mole_fraction_cutoff = 1.0e-14;

    /// The cutoff value for normalized species amounts.
    /// This parameter is used to ignore certain species with tiny amounts
    /// during the error test. Those with normalized amounts, \eq{n_{i}/\text{sum}(n)},
    /// below this threshold will be discarded from the error control.
    double amount_fraction_cutoff = 1.0e-14;

    /// The weighted residual tolerance used in the acceptance test.
    double tol = 1.0e-1;

    /// The relative tolerance used in the acceptance test.
    double reltol = 0.1;

    /// The absolute tolerance used in the acceptance test.
    double abstol = 1e-8;

    /// The small negative cutoff value for estimated species mole amounts.
    double cutoff = -1e-5;
};

} // namespace Reaktoro
