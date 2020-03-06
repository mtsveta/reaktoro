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

#include "KineticSolver.hpp"

// C++ includes
#include <functional>
using namespace std::placeholders;

// Reaktoro includes
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Profiling.hpp>
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumSolver.hpp>
#include <Reaktoro/Kinetics/KineticOptions.hpp>
#include <Reaktoro/Kinetics/KineticResult.hpp>
#include <Reaktoro/Kinetics/KineticProblem.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>

namespace Reaktoro {

struct KineticSolver::Impl
{
    /// The kinetically-controlled chemical reactions
    ReactionSystem reactions;

    /// The chemical system instance
    ChemicalSystem system;

    /// The partition of the species in the chemical system
    Partition partition;

    /// The options of the kinetic solver
    KineticOptions options;

    /// The result of the kinetic solver
    KineticResult result;

    /// The equilibrium solver instance
    EquilibriumSolver equilibrium;

    /// The equilibrium solver instance
    SmartEquilibriumSolver smart_equilibrium;

    /// The sensitivity of the equilibrium state
    EquilibriumSensitivity sensitivity;

    /// The ODE solver instance
    ODESolver ode;

    /// The indices of the equilibrium and kinetic species
    Indices ies, iks;

    /// The indices of the elements in the equilibrium and kinetic partition
    Indices iee, ike;

    /// The number of equilibrium and kinetic species
    Index Ne, Nk;

    /// The number of elements in the equilibrium and kinetic partition
    Index Ee, Ek;

    /// The formula matrix of the equilibrium species
    Matrix Ae, Ak;

    /// The stoichiometric matrix w.r.t. the equilibrium and kinetic species
    Matrix Se, Sk;

    /// The coefficient matrix `A` of the kinetic rates
    Matrix A;

    /// The coefficient matrix `B` of the source rates
    Matrix B;

    /// The temperature of the chemical system (in units of K)
    double T;

    /// The pressure of the chemical system (in units of Pa)
    double P;

    /// The molar composition of the equilibrium and kinetic species
    Vector ne, nk;

    /// The molar abundance of the elements in the equilibrium species
    Vector be;

    /// The combined vector of elemental molar abundance and composition of kinetic species [be nk]
    Vector benk;

    /// The chemical properties of the system
    ChemicalProperties properties;

    /// The vector with the values of the reaction rates
    ChemicalVector r;

    /// The partial derivatives of the reaction rates `r` w.r.t. to `be`, `ne`, `nk`, `and `u = [be nk]`
    Matrix drdbe, drdne, drdnk, drdu;

    /// The source term
    ChemicalVector q;

    /// The partial derivatives of the source rates `q` w.r.t. to `be`, `ne`, `nk`, `and `u = [be nk]`
    Matrix dqdbe, dqdne, dqdnk, dqdu;

    /// The function that calculates the source term in the problem
    std::function<ChemicalVector(const ChemicalProperties&)> source_fn;

    /// Construct an instance of KineticSolver::Impl
    Impl(const ReactionSystem& reactions, const Partition& partition)
    : reactions(reactions),
      system(partition.system()),
      partition(partition),
      equilibrium(partition)
    {
        // Set the indices of the equilibrium and kinetic species
        ies = partition.indicesEquilibriumSpecies();
        iks = partition.indicesKineticSpecies();

        // Set the indices of the equilibrium and kinetic elements
        iee = partition.indicesEquilibriumElements();
        ike = partition.indicesKineticElements();

        // Set the number of equilibrium and kinetic species
        Ne = ies.size();
        Nk = iks.size();

        // Set the number of equilibrium and kinetic elements
        Ee = iee.size();
        Ek = ike.size();

        // Initialise the formula matrix of the equilibrium partition
        Ae = partition.formulaMatrixEquilibriumPartition();

        // Initialise the formula matrix of the kinetic partition
        // Ak_tmp is of the size (Indices of elements that are included in kinetic species) x Nk
        Matrix Ak_tmp = partition.formulaMatrixKineticPartition();
        if(Nk)
        {
            // Initialize formula matrix of the kinetic partition
            // Ak is of the size N x Nk
            Ak = Matrix::Zero(system.numElements(), Nk);

            for(Index i = 0; i < ike.size(); ++i)
            {
                // Copy the rows of Ak_tmp to the positions of the kinetic elements
                Ak.row(ike[i]) << Ak_tmp.row(i);
            }
        }

        // Initialise the stoichiometric matrices w.r.t. the equilibrium and kinetic species
        Se = cols(reactions.stoichiometricMatrix(), ies);
        Sk = cols(reactions.stoichiometricMatrix(), iks);

        // Initialise the coefficient matrix `A` of the kinetic rates
        A.resize(Ee + Nk, reactions.numReactions());
        A.topRows(Ee) = Ae * tr(Se);
        A.bottomRows(Nk) = tr(Sk);

        // Auxiliary identity matrix
        const Matrix I = identity(Ne + Nk, Ne + Nk);

        // Auxiliary selected equilibrium and kinetic rows of the identity matrix
        const Matrix Ie = rows(I, ies);
        const Matrix Ik = rows(I, iks);

        // Initialise the coefficient matrix `B` of the source rates
        B = zeros(Ee + Nk, system.numSpecies());
        B.topRows(Ee) = Ae * Ie;
        B.bottomRows(Nk) = Ik;

        // Allocate memory for the partial derivatives of the reaction rates `r` w.r.t. to `u = [be nk]`
        drdu.resize(reactions.numReactions(), Ee + Nk);

        // Allocate memory for the partial derivatives of the source rates `q` w.r.t. to `u = [be nk]`
        dqdu.resize(system.numSpecies(), Ee + Nk);

        if (Nk)
        {
            /*
            //* TODO: remove when debugging is finished
            std::cout << "ies :"; for (auto elem : ies) std::cout << elem << " "; std::cout << std::endl;
            std::cout << "iks :"; for (auto elem : iks) std::cout << elem << " "; std::cout << std::endl;


            std::cout << "setPartition function :" << std::endl;

            std::cout << "iee :"; for (auto elem : iee) std::cout << elem << " "; std::cout << std::endl;
            std::cout << "ike :"; for (auto elem : ike) std::cout << elem << " "; std::cout << std::endl;

            std::cout << "S \n" << reactions.stoichiometricMatrix() << std::endl;
            std::cout << "Se \n" << Se << std::endl;
            std::cout << "Sk \n" << Sk << std::endl;

            std::cout << "A \n" << A << std::endl;
            std::cout << "Ae \n" << Ae << std::endl;
            std::cout << "Ak \n" << Ak << std::endl;

            std::cout << "B \n" << B << std::endl;

            std::cout << "Ee " << Ee << std::endl;
            std::cout << "Nk " << Nk << std::endl;
            */
        }
    }

    auto setOptions(const KineticOptions& options_) -> void
    {
        // Initialise the options of the kinetic solver
        options = options_;
    }

    auto addSource(ChemicalState state, double volumerate, std::string units) -> void
    {
        const Index num_species = system.numSpecies();
        const double volume = units::convert(volumerate, units, "m3/s");
        state.scaleVolume(volume);
        const Vector n = state.speciesAmounts();
        auto old_source_fn = source_fn;

        source_fn = [=](const ChemicalProperties& properties)
        {
            ChemicalVector q(num_species);
            q.val = n;
            if(old_source_fn)
                q += old_source_fn(properties);
            return q;
        };
    }

    auto addPhaseSink(std::string phase, double volumerate, std::string units) -> void
    {
        const double volume = units::convert(volumerate, units, "m3/s");
        const Index iphase = system.indexPhaseWithError(phase);
        const Index ifirst = system.indexFirstSpeciesInPhase(iphase);
        const Index size = system.numSpeciesInPhase(iphase);
        auto old_source_fn = source_fn;
        ChemicalScalar phasevolume;
        ChemicalVector q(size);

        source_fn = [=](const ChemicalProperties& properties) mutable
        {
            const auto n = properties.composition();
            const auto np = rows(n, ifirst, size);
            auto qp = rows(q, ifirst, size);
            phasevolume = properties.phaseVolumes()[iphase];
            qp = -volume*np/phasevolume;
            if(old_source_fn)
                q += old_source_fn(properties);
            return q;
        };
    }

    auto addFluidSink(double volumerate, std::string units) -> void
    {
        const double volume = units::convert(volumerate, units, "m3/s");
        const Indices& isolid_species = partition.indicesSolidSpecies();
        auto old_source_fn = source_fn;
        ChemicalScalar fluidvolume;
        ChemicalVector q;

        source_fn = [=](const ChemicalProperties& properties) mutable
        {
            const auto n = properties.composition();
            fluidvolume = properties.fluidVolume();
            q = -volume*n/fluidvolume;
            rows(q, isolid_species).fill(0.0);
            if(old_source_fn)
                q += old_source_fn(properties);
            return q;
        };
    }

    auto addSolidSink(double volumerate, std::string units) -> void
    {
        const double volume = units::convert(volumerate, units, "m3/s");
        const Indices& ifluid_species = partition.indicesFluidSpecies();
        auto old_source_fn = source_fn;
        ChemicalScalar solidvolume;
        ChemicalVector q;

        source_fn = [=](const ChemicalProperties& properties) mutable
        {
            const auto n = properties.composition();
            solidvolume = properties.solidVolume();
            q = -volume*n/solidvolume;
            rows(q, ifluid_species).fill(0.0);
            if(old_source_fn)
                q += old_source_fn(properties);
            return q;
        };
    }

    auto initialize(ChemicalState& state, double tstart) -> void
    {
        // Initialise the temperature and pressure variables
        T = state.temperature();
        P = state.pressure();

        // Extract the composition of the equilibrium and kinetic species
        const auto& n = state.speciesAmounts();
        ne = n(ies);
        nk = n(iks);

        // Assemble the vector benk = [be nk]
        benk.resize(Ee + Nk);
        benk.head(Ee) = Ae * ne;
        benk.tail(Nk) = nk;

        // Initialize
        initialize(state, tstart, benk);

    }
    auto initialize(ChemicalState& state, double tstart, VectorConstRef benk) -> void
    {
        // Initialise the temperature and pressure variables
        T = state.temperature();
        P = state.pressure();

        // Define the ODE function
        ODEFunction ode_function = [&](double t, VectorConstRef u, VectorRef res)
        {
            return function(state, t, u, res);
        };

        // Define the jacobian of the ODE function
        ODEJacobian ode_jacobian = [&](double t, VectorConstRef u, MatrixRef res)
        {
            return jacobian(state, t, u, res);
        };

        // Initialise the ODE problem
        ODEProblem problem;
        problem.setNumEquations(Ee + Nk);
        problem.setFunction(ode_function);
        problem.setJacobian(ode_jacobian);

        // Define the options for the ODE solver
        ODEOptions options_ode = options.ode;

        // Set the ODE problem and initialize the ODE solver
        ode.setProblem(problem);
        ode.setOptions(options_ode);
        ode.initialize(tstart, benk);

        // Set the options of the equilibrium solver
        equilibrium.setOptions(options.equilibrium);
        smart_equilibrium.setOptions(options.smart_equilibrium);
    }

    auto step(ChemicalState& state, double t) -> double
    {
        const double tfinal = Index(-1);
        step(state, t, tfinal);
        return t;
    }

    auto step(ChemicalState& state, double t, double tfinal) -> double
    {
        // Extract the composition vector of the equilibrium and kinetic species
        const auto& n = state.speciesAmounts();
        ne = n(ies);
        nk = n(iks);

        // Assemble the vector benk = [be nk]
        benk.head(Ee) = Ae * ne;
        benk.tail(Nk) = nk;

        // Perform one ODE step integration
        ode.integrate(t, benk, tfinal);

        // Extract the `be` and `nk` entries of the vector `benk`
        be = benk.head(Ee);
        nk = benk.tail(Nk);

        // Update the composition of the kinetic species
        state.setSpeciesAmounts(nk, iks);

        // Update the composition of the equilibrium species
        if(options.use_smart_equilibrium_solver) smart_equilibrium.solve(state, T, P, be);
        else equilibrium.solve(state, T, P, be);

        return t;
    }

    auto solve(ChemicalState& state, double t, double dt, VectorConstRef b) -> void
    {
        // Initialize result structure
        result = {};

        tic(0);

        // Extract the composition of the kinetic species from the state
        const auto &n = state.speciesAmounts();
        nk = n(iks);

        // Assemble the vector benk = [be nk]
        benk.resize(Ee + Nk);
        benk.head(Ee) = b - Ak * nk;
        benk.tail(Nk) = nk;

        // Initialise the chemical kinetics solver
        timeit(initialize(state, t, benk), result.timing.initialize=);

        //std::cout << "benk before  ..." << tr(benk) << std::endl;

        // Integrate the chemical kinetics ODE from `t` to `t + dt`
        timeit(ode.solve(t, dt, benk), result.timing.integrate=);
        //timeit(ode.solve_explicit_1st_order(t, dt, benk), result.timing.integrate=);
        //timeit(ode.solve_explicit_2nd_order(t, dt, benk), result.timing.integrate=);
        //timeit(ode.solve_implicit_1st_order(t, dt, benk), result.timing.integrate=);

        // Extract the `be` and `nk` entries of the vector `benk`
        be = benk.head(Ee);
        nk = benk.tail(Nk);

        //std::cout << "state.setSpeciesAmounts(nk, iks)  ..." << std::endl;

        // Update the composition of the kinetic species
        state.setSpeciesAmounts(nk, iks);

        // Equilibrate equilibrium species
        tic(1);

        if(options.use_smart_equilibrium_solver){
            SmartEquilibriumResult res = {};
            res += smart_equilibrium.solve(state, T, P, be);
            result.smart_equilibrium += res;
        }
        else{
            EquilibriumResult res = {};
            res += equilibrium.solve(state, T, P, be);
            result.equilibrium += res;
        }

        toc(1, result.timing.equilibrate);

        toc(0, result.timing.solve);

    }
    auto solve(ChemicalState& state, double t, double dt, VectorConstRef b, Index step, Index icell) -> void
    {
        // Initialize result structure
        result = {};

        tic(0);

        // Extract the composition of the kinetic species from the state
        const auto &n = state.speciesAmounts();
        nk = n(iks);

        // Assemble the vector benk = [be nk]
        benk.resize(Ee + Nk);
        benk.head(Ee) = b - Ak * nk;
        benk.tail(Nk) = nk;

        // Initialise the chemical kinetics solver
        timeit(initialize(state, t, benk), result.timing.initialize=);

        /*
        if(step >= 178)
        {

            std::cout << "step " << step << ", icell " << icell << std::endl;
            std::cout << "benk before  ..." << tr(benk) << std::endl;
        }
        */

        // Integrate the chemical kinetics ODE from `t` to `t + dt`
        //timeit(ode.solve(t, dt, benk), result.timing.integrate=);
        //timeit(ode.solve_explicit_1st_order(t, dt, benk), result.timing.integrate=);
        //timeit(ode.solve_explicit_2nd_order(t, dt, benk), result.timing.integrate=);
        timeit(ode.solve_implicit_1st_order(t, dt, benk), result.timing.integrate=);

        /*
        if(step > 178)
        {
            std::cout << "benk after  ..." << tr(benk) << std::endl;
        }
        */
        // Extract the `be` and `nk` entries of the vector `benk`
        be = benk.head(Ee);
        nk = benk.tail(Nk);

        //std::cout << "state.setSpeciesAmounts(nk, iks)  ..." << std::endl;

        // Update the composition of the kinetic species
        state.setSpeciesAmounts(nk, iks);

        // Equilibrate equilibrium species
        tic(1);

        if(options.use_smart_equilibrium_solver){
            SmartEquilibriumResult res = {};
            res += smart_equilibrium.solve(state, T, P, be);
            result.smart_equilibrium += res;
        }
        else{
            EquilibriumResult res = {};
            res += equilibrium.solve(state, T, P, be);
            result.equilibrium += res;
        }

        toc(1, result.timing.equilibrate);

        toc(0, result.timing.solve);

    }

    auto function(ChemicalState& state, double t, VectorConstRef u, VectorRef res) -> int
    {
        // Extract the `be` and `nk` entries of the vector [be, nk]
        be = u.head(Ee);
        nk = u.tail(Nk);

        // Check for non-finite values in the vector `benk`
        for(Index i = 0; i < u.rows(); ++i)
            if(!std::isfinite(u[i]))
                return 1; // ensure the ode solver will reduce the time step

        // Update the composition of the kinetic species in the member `state`
        state.setSpeciesAmounts(nk, iks);

        // Solve the equilibrium problem using the elemental molar abundance `be`
        //if(options.use_smart_equilibrium_solver) // using smart equilibrium solver
        if(false) // block using smart equilibrium solver (for now)
        {
            SmartEquilibriumResult res = {};

            // smart_equilibrium_result
            timeit(res += smart_equilibrium.solve(state, T, P, be), result.timing.integrate_equilibration+=);

            /*
            std::cout << " - learn    : " << res.timing.learn << std::endl;
            std::cout << " - estimate : " << res.timing.estimate << std::endl;
            std::cout << "   - search : " << res.timing.estimate_search << std::endl;
            getchar();
            */
            result.smart_equilibrium += res;

            // If smart calculation failed, use cold-start
            if(!res.estimate.accepted && !res.learning.gibbs_energy_minimization.optimum.succeeded)
            {
                //std::cout << "restart smart_equilibrium: " << res.learning.gibbs_energy_minimization.optimum.succeeded << std::endl;

                state.setSpeciesAmounts(0.0);
                timeit( res = smart_equilibrium.solve(state, T, P, be), result.timing.integrate_equilibration+=);

            }

            // Assert the smart equilibrium calculation did not fail
            Assert(res.estimate.accepted || res.learning.gibbs_energy_minimization.optimum.succeeded,
                   "Could not calculate the rates of the species.",
                   "The smart equilibrium calculation failed.");

        }
        else  // using conventional equilibrium solver
        {
            EquilibriumResult res = {};

            timeit(res += equilibrium.solve(state, T, P, be), result.timing.integrate_equilibration +=);

            result.equilibrium += res;

            // Check if the calculation failed, if so, use cold-start
            if (!res.optimum.succeeded) {
                //std::cout << "t : " << t << std::endl;
                //std::cout << "n : " << tr(state.speciesAmounts()) << std::endl;
                state.setSpeciesAmounts(0.0);
                timeit(res = equilibrium.solve(state, T, P, be), result.timing.integrate_equilibration +=);

            }
            if (!res.optimum.succeeded) {
                std::cout << "t : " << t << std::endl;
                std::cout << "n : " << tr(state.speciesAmounts()) << std::endl;
                //std::cout << "n : " << state.speciesAmount("H+") << std::endl;
                //std::cout << "error : " << res.optimum.error << std::endl;
                //std::cout << "iterations : " << res.optimum.iterations << std::endl;
                return 1; // ensure the ode solver will reduce the time step
            }
            // Assert the equilibrium calculation did not fail
            Assert(res.optimum.succeeded,
                   "Could not calculate the rates of the species.",
                   "The equilibrium calculation failed.");
        }


        //std::cout << "properties ..." << std::endl;
        // Update the chemical properties of the system
        timeit(properties = state.properties(), result.timing.integrate_chemical_properties+=);

        //std::cout << "reactions ..." << std::endl;
        // Calculate the kinetic rates of the reactions
        timeit(r = reactions.rates(properties), result.timing.integrate_reaction_rates+=);

        // Calculate the right-hand side function of the ODE
        res = A * r.val;

        // Add the function contribution from the source rates
        if(source_fn)
        {
            // Evaluate the source function
            q = source_fn(properties);

            // Add the contribution of the source rates
            res += B * q.val;
        }

        return 0;
    }

    auto jacobian(ChemicalState& state, double t, VectorConstRef u, MatrixRef res) -> int
    {
        // Calculate the sensitivity of the equilibrium state
        sensitivity = equilibrium.sensitivity();

        // Extract the columns of the kinetic rates derivatives w.r.t. the equilibrium and kinetic species
        drdne = cols(r.ddn, ies);
        drdnk = cols(r.ddn, iks);

        // The sensitivity derivatives with respect to equilibrium species
        const auto dnedbe = sensitivity.dndb(ies, iee);

        // Calculate the derivatives of `r` w.r.t. `be` using the equilibrium sensitivity
        drdbe = drdne * dnedbe;

        // Assemble the partial derivatives of the reaction rates `r` w.r.t. to `u = [be nk]`
        drdu << drdbe, drdnk;

        // Calculate the Jacobian matrix of the ODE function
        res = A * drdu;

        // Add the Jacobian contribution from the source rates
        if(source_fn)
        {
            // Extract the columns of the source rates derivatives w.r.t. the equilibrium and kinetic species
            dqdne = cols(q.ddn, ies);
            dqdnk = cols(q.ddn, iks);

            // Calculate the derivatives of `q` w.r.t. `be` using the equilibrium sensitivity
            dqdbe = dqdne * dnedbe;

            // Assemble the partial derivatives of the source rates `q` w.r.t. to `u = [be nk]`
            dqdu << dqdbe, dqdnk;

            // Add the contribution of the source rates
            res += B * dqdu;
        }

        return 0;
    }

};

KineticSolver::KineticSolver()
{
    RuntimeError("Cannot proceed with KineticSolver().",
        "KineticSolver() constructor is deprecated. "
        "Use constructor KineticSolver(const ReactionSystem&, const Partition&) instead.");
}

KineticSolver::KineticSolver(const KineticSolver& other)
: pimpl(new Impl(*other.pimpl))
{}

KineticSolver::KineticSolver(const ReactionSystem& reactions)
{
    RuntimeError("Cannot proceed with KineticSolver(const ReactionSystem&).",
        "KineticSolver(const ReactionSystem&) constructor is deprecated. "
        "Use constructor KineticSolver(const ReactionSystem&, const Partition&) instead.");
}

KineticSolver::KineticSolver(const ReactionSystem& reactions, const Partition& partition)
: pimpl(new Impl(reactions, partition))
{}

KineticSolver::~KineticSolver()
{}

auto KineticSolver::operator=(KineticSolver other) -> KineticSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto KineticSolver::setOptions(const KineticOptions& options) -> void
{
    pimpl->setOptions(options);
}

auto KineticSolver::setPartition(const Partition& partition) -> void
{
    RuntimeError("Cannot proceed with KineticSolver::setPartition.",
        "KineticSolver::setPartition is deprecated. "
        "Use constructor KineticSolver(const ReactionSystem&, const Partition&) instead.");
}

auto KineticSolver::addSource(const ChemicalState& state, double volumerate, std::string units) -> void
{
    pimpl->addSource(state, volumerate, units);
}

auto KineticSolver::addPhaseSink(std::string phase, double volumerate, std::string units) -> void
{
    pimpl->addPhaseSink(phase, volumerate, units);
}

auto KineticSolver::addFluidSink(double volumerate, std::string units) -> void
{
    pimpl->addFluidSink(volumerate, units);
}

auto KineticSolver::addSolidSink(double volumerate, std::string units) -> void
{
    pimpl->addSolidSink(volumerate, units);
}

auto KineticSolver::initialize(ChemicalState& state, double tstart) -> void
{
    pimpl->initialize(state, tstart);
}

auto KineticSolver::initialize(ChemicalState& state, double tstart, VectorConstRef benk) -> void
{
    pimpl->initialize(state, tstart, benk);
}

auto KineticSolver::step(ChemicalState& state, double t) -> double
{
    return pimpl->step(state, t);
}

auto KineticSolver::step(ChemicalState& state, double t, double dt) -> double
{
    return pimpl->step(state, t, dt);
}

auto KineticSolver::solve(ChemicalState& state, double t, double dt, VectorConstRef b) -> void
{
    pimpl->solve(state, t, dt, b);
}

auto KineticSolver::solve(ChemicalState& state, double t, double dt, VectorConstRef b, Index step, Index icell) -> void{
    pimpl->solve(state, t, dt, b, step, icell);
}

auto KineticSolver::properties() const -> const ChemicalProperties&
{
    return pimpl->properties;
}

auto KineticSolver::result() const -> const KineticResult&
{
    return pimpl->result;
}


} // namespace Reaktoro
