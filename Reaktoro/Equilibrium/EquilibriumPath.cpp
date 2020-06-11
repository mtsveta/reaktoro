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

#include "EquilibriumPath.hpp"

// C++ includes
#include <list>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalOutput.hpp>
#include <Reaktoro/Core/ChemicalPlot.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Math/ODE.hpp>

namespace Reaktoro {

struct EquilibriumPath::Impl
{
    /// The chemical system
    ChemicalSystem system;

    /// The options for the equilibrium path calculation and visualization
    EquilibriumPathOptions options;

    /// The partition of the chemical system
    Partition partition;

    /// The output instance of the equilibrium path calculation
    ChemicalOutput output;

    /// The plots of the equilibrium path calculation
    std::vector<ChemicalPlot> plots;

    /// Construct a EquilibriumPath::Impl instance
    explicit Impl(const ChemicalSystem& system)
    : Impl(Partition(system))
    {}

    /// Construct a EquilibriumPath::Impl instance
    explicit Impl(const Partition& partition)
    : system(partition.system()), partition(partition)
    {}

    /// Set the options for the equilibrium path calculation and visualization
    auto setOptions(const EquilibriumPathOptions& options_) -> void
    {
        options = options_;
    }

    /// Solve the path of equilibrium states between two chemical states
    auto solve(const ChemicalState& state_i, const ChemicalState& state_f) -> EquilibriumPathResult
    {
        // The result of this equilibrium path calculation
        EquilibriumPathResult result;

        // The number of equilibrium species
        const unsigned Ne = partition.numEquilibriumSpecies();

        // The indices of species and elements in the equilibrium partition
        const auto& ies = partition.indicesEquilibriumSpecies();
        const auto& iee = partition.indicesEquilibriumElements();

        /// The temperatures at the initial and final chemical states
        const double T_i = state_i.temperature();
        const double T_f = state_f.temperature();

        /// The pressures at the initial and final chemical states
        const double P_i = state_i.pressure();
        const double P_f = state_f.pressure();

        /// The molar amounts of the elements in the equilibrium partition at the initial and final chemical states
        const Vector be_i = state_i.elementAmountsInSpecies(ies)(iee);
        const Vector be_f = state_f.elementAmountsInSpecies(ies)(iee);

        // The equilibrium solver
        EquilibriumSolver equilibrium(partition);
        equilibrium.setOptions(options.equilibrium);

        /// The sensitivity of the equilibrium state
        EquilibriumSensitivity sensitivity;

        // The chemical state updated throughout the path calculation
        ChemicalState state = state_i;

        // The ODE function describing the equilibrium path
        ODEFunction f = [&](double t, VectorConstRef ne, VectorRef res) -> int
        {
            static double tprev = t;

            // Skip if t is greater or equal than 1
            if(t > 0.0 && std::abs(t - tprev) >= options.maxstep) return 1;

            // Update tprev
            tprev = t;

            // Skip if t is greater or equal than 1
            if(t >= 1) return 0;

            // Calculate T, P, be at current t
            const double T  = T_i + t * (T_f - T_i);
            const double P  = P_i + t * (P_f - P_i);
            const Vector be = be_i + t * (be_f - be_i);

            // Perform the equilibrium calculation at T, P, be
            result.equilibrium += equilibrium.solve(state, T, P, be);

            // Calculate the sensitivity of the equilibrium state
            sensitivity = equilibrium.sensitivity();

            // The sensitivity derivatives with respect to equilibrium species
            const auto dnedT  = sensitivity.dndT(ies);
            const auto dnedP  = sensitivity.dndP(ies);
            const auto dnedbe = sensitivity.dndb(ies, iee);

            // Check if the calculation succeeded
            if(!result.equilibrium.optimum.succeeded) return 1;

            // Calculate the right-hand side vector of the ODE
            res = dnedT*(T_f - T_i) + dnedP*(P_f - P_i) + dnedbe*(be_f - be_i);

            return 0;
        };

        // Initialize the ODE problem
        ODEProblem problem;
        problem.setNumEquations(Ne);
        problem.setFunction(f);

        // The initial and final molar amounts of equilibrium species
        Vector ne_i = state_i.speciesAmounts()(ies);
        Vector ne_f = state_i.speciesAmounts()(ies);

        // Adjust the absolute tolerance parameters for each component
        options.ode.abstols = options.ode.abstol * ((ne_i + ne_f)/2.0 + 1.0);

        // Ensure the iteration algorithm is not Newton
        options.ode.iteration = ODEIterationMode::Functional;

        ODESolver ode;
        ode.setOptions(options.ode);
        ode.setProblem(problem);

        // Initialize initial conditions
        double t = 0.0;
        Vector& ne = ne_i;

        // Initialize the ODE solver
        ode.initialize(t, ne);

        // Initialize the output of the equilibrium path calculation
        if(output) output.open();

        // Initialize the plots of the equilibrium path calculation
        for(auto& plot : plots) plot.open();

        // Perform the integration from t = 0 to t = 1
        while(t < 1.0)
        {
            // Update the output with current state
            if(output) output.update(state, t);

            // Update the plots with current state
            for(auto& plot : plots) plot.update(state, t);

            // Integrate one time step only
            ode.integrate(t, ne);
        }

        // Update the output with the final state
        if(output) output.update(state_f, 1.0);

        // Update the plots with the final state
        for(auto& plot : plots) plot.update(state_f, 1.0);

        return result;
    }
};

EquilibriumPath::EquilibriumPath(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

EquilibriumPath::EquilibriumPath(const EquilibriumPath& other)
: pimpl(new Impl(*other.pimpl))
{}

EquilibriumPath::~EquilibriumPath()
{}

auto EquilibriumPath::operator=(EquilibriumPath other) -> EquilibriumPath&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumPath::setOptions(const EquilibriumPathOptions& options) -> void
{
    pimpl->setOptions(options);
}

auto EquilibriumPath::setPartition(const Partition& partition) -> void
{
    RuntimeError("Cannot proceed with EquilibriumPath::setPartition.",
        "EquilibriumPath::setPartition is deprecated. "
        "Use constructor EquilibriumPath(const Partition&) instead.");
}

auto EquilibriumPath::solve(const ChemicalState& state_i, const ChemicalState& state_f) -> EquilibriumPathResult
{
    return pimpl->solve(state_i, state_f);
}

auto EquilibriumPath::output() -> ChemicalOutput
{
    pimpl->output = ChemicalOutput(pimpl->system);
    return pimpl->output;
}

auto EquilibriumPath::plot() -> ChemicalPlot
{
    pimpl->plots.push_back(ChemicalPlot(pimpl->system));
    return pimpl->plots.back();
}

auto EquilibriumPath::plots(unsigned num) -> std::vector<ChemicalPlot>
{
    for(unsigned i = 0; i < num; ++i) plot();
    return pimpl->plots;
}

auto EquilibriumPath::system() const -> const ChemicalSystem&
{
    return pimpl->system;
}

auto EquilibriumPath::partition() const -> const Partition&
{
    return pimpl->partition;
}

} // namespace Reaktoro
