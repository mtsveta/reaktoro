// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018
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

namespace Reaktoro {

/// The options for the smart kinetics calculations.
/// @see SmartKineticSolver
struct SmartKineticOptions
{
    /// The kinetic options to be used during learning operations.
    KineticOptions learning;

    /// The relative tolerance for estimated species mole amounts.
    double reltol = 1.0;

    /// The absolute tolerance for estimated species mole amounts.
    double abstol = 1e-14;
};

} // namespace Reaktoro
