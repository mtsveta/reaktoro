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
#include <Reaktoro/Optimization/OptimumSolverBase.hpp>

namespace Reaktoro {

// Forward declarations
struct OptimumOptions;
struct OptimumProblem;
struct OptimumResult;
struct OptimumState;

/// The class that implements the IpAction algorithm using an interior-point method
/// with a null-space KKT formulation.
class OptimumSolverIpAction : public OptimumSolverBase
{
public:
    /// Construct a default OptimumSolverIpAction instance.
    OptimumSolverIpAction();

    /// Construct a copy of an OptimumSolverIpAction instance.
    OptimumSolverIpAction(const OptimumSolverIpAction& other);

    /// Destroy this OptimumSolverIpAction instance.
    virtual ~OptimumSolverIpAction();

    /// Assign an OptimumSolverIpAction instance to this.
    auto operator=(OptimumSolverIpAction other) -> OptimumSolverIpAction&;

    /// Solve an optimisation problem.
    /// @param problem The definition of the optimisation problem
    /// @param state[in,out] The initial guess and the final state of the optimisation calculation
    virtual auto solve(const OptimumProblem& problem, OptimumState& state) -> OptimumResult;

    /// Solve an optimisation problem with given options.
    /// @param problem The definition of the optimisation problem
    /// @param state[in,out] The initial guess and the final state of the optimisation calculation
    /// @param options The options for the optimisation calculation
    virtual auto solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult;

    /// Return a clone of this instance.
    virtual auto clone() const -> OptimumSolverBase*;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
