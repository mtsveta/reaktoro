// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2019
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

namespace Reaktoro {

/// Timing information of the operations during an kinetic calculation.
struct SmartKineticTiming {

    /// The time spent for solving the chemical kinetic problem.
    double solve = 0.0;

    /// The time spent for learning the chemical kinetic path
    double learn = 0.0;

    /// The time spent for estimating the chemical kinetic path
    double estimate = 0.0;

    /// Self addition assignment to accumulate kinetic timing.
    auto operator+=(const SmartKineticTiming& other) -> SmartKineticTiming&;

};
/// A type used to describe the result of an kinetic calculation.
struct SmartKineticResult
{
    /// The timing information of the operations during an equilibrium calculation.
    SmartKineticTiming timing;

    /// Self addition assignment to accumulate results.
    auto operator+=(const SmartKineticResult& other) -> SmartKineticResult&;
};

} // namespace Reaktoro
