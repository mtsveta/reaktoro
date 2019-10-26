// This file is part of Reaktoro (https://reaktoro.org).
//
// Reaktoro is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// Reaktoro is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include "ReactiveTransportSolver.hpp"

// C++ includes
#include <algorithm>

// Reaktoro includes
#include <Reaktoro/Common/Profiling.hpp>

namespace Reaktoro {

struct ReactiveTransportSolver::Impl
{
    /// The chemical system common to all degrees of freedom in the chemical field.
    ChemicalSystem system;

    /// The solver for solving the transport equations
    TransportSolver transport_solver;

    /// The options for the reactive transport calculations.
    ReactiveTransportOptions options;

    /// The result information of the last reactive transport time step calculation.
    ReactiveTransportResult result;

    /// The equilibrium solver using conventional Gibbs energy minimization approach.
    EquilibriumSolver equilibrium_solver;

    /// The equilibrium solver using a smart on-demand learning strategy.
    SmartEquilibriumSolver smart_equilibrium_solver;

    /// The equilibrium solver using conventional Gibbs energy minimization approach.
    KineticSolver kinetic_solver;

    /// The equilibrium solver using a smart on-demand learning strategy.
    SmartKineticSolver smart_kinetic_solver;

    /// The list of chemical output objects
    std::vector<ChemicalOutput> outputs;

    /// The amounts of fluid elements on the boundary.
    Vector bbc;

    /// The amounts of a fluid element on each cell of the mesh.
    Matrix bf;

    /// The amounts of a solid element on each cell of the mesh.
    Matrix bs;

    /// The amounts of an element on each cell of the mesh.
    Matrix b;

    /// The current number of steps in the solution of the reactive transport equations.
    Index steps = 0;

    /// Construct a ReactiveTransportSolver::Impl instance.
    Impl(const ChemicalSystem& system)
    : system(system), equilibrium_solver(system), smart_equilibrium_solver(system)
    {
        setBoundaryState(ChemicalState(system));
    }

    Impl(const ChemicalSystem& system, const ReactionSystem& reactions, const Partition& partition)
    : system(system), equilibrium_solver(system), smart_equilibrium_solver(system),
    kinetic_solver(reactions), smart_kinetic_solver(reactions)
    {
        setBoundaryState(ChemicalState(system));
        kinetic_solver.setPartition(partition);
        smart_kinetic_solver.setPartition(partition);
        equilibrium_solver.setPartition(partition);
        smart_equilibrium_solver.setPartition(partition);
    }

    /// Set the options for the reactive transport calculations.
    auto setOptions(const ReactiveTransportOptions& options) -> void
    {
        this->options = options;

        // Set options of equilibrium solvers
        equilibrium_solver.setOptions(options.equilibrium);
        smart_equilibrium_solver.setOptions(options.smart_equilibrium);

        // Set options of kinetic solvers
        kinetic_solver.setOptions(options.kinetics);
        this->options.kinetics.use_smart_equilibrium_solver = options.use_smart_equilibrium_solver; // TODO: make sure the flag is set

        // Set options of smart kinetic solvers
        smart_kinetic_solver.setOptions(options.smart_kinetics);
        this->options.smart_kinetics.use_smart_equilibrium_solver = options.use_smart_equilibrium_solver; // TODO: make sure the flag is set
    }

    /// Initialize the mesh discretizing the computational domain for reactive transport.
    auto setMesh(const Mesh& mesh) -> void
    {
        transport_solver.setMesh(mesh);
    }

    /// Initialize the velocity of the reactive transport model.
    auto setVelocity(double val) -> void
    {
        transport_solver.setVelocity(val);
    }

    /// Initialize the diffusion of the reactive transport model.
    auto setDiffusionCoeff(double val) -> void
    {
        transport_solver.setDiffusionCoeff(val);
    }

    /// Initialize boundary conditions of the reactive transport model.
    auto setBoundaryState(const ChemicalState& state) -> void
    {
        bbc = state.elementAmounts();
    }

    /// Initialize time step of the reactive transport sequential algorithm.
    auto setTimeStep(double val) -> void
    {
        transport_solver.setTimeStep(val);
    }

    /// Add the output to the reactive transport modelling.
    auto output() -> ChemicalOutput
    {
        outputs.emplace_back(ChemicalOutput(system));
        return outputs.back();
    }

    /// Initialize the reactive transport solver.
    auto initialize() -> void
    {
        // Initialize mesh and corresponding amount e
        const Mesh& mesh = transport_solver.mesh();
        const Index num_elements = system.numElements();
        const Index num_cells = mesh.numCells();

        // Initialize amount of elements in fluid and solid phases
        bf.resize(num_cells, num_elements);
        bs.resize(num_cells, num_elements);
        b.resize(num_cells, num_elements);

        // Initialize equilibrium solver based on the parameter
        transport_solver.setOptions(options.transport);
        transport_solver.initialize();

    }

    /// Perform one time step of a reactive transport calculation.
    auto step(ChemicalField& field) -> ReactiveTransportResult
    {
        // The result of the reactive transport step
        ReactiveTransportResult rt_result = {};

        // Auxiliary variables
        const auto& mesh = transport_solver.mesh();
        const auto& num_elements = system.numElements();
        const auto& num_cells = mesh.numCells();
        const auto& ifs = system.indicesFluidSpecies();
        const auto& iss = system.indicesSolidSpecies();
        auto& states = field.states();
        auto& properties = field.properties();

        // Open the the file for outputting chemical states
        for(auto output : outputs)
        {
            output.suffix("-" + std::to_string(steps));
            output.open();
        }

        //---------------------------------------------------------------------------
        // Step 1: Perform a time step transport calculation for each fluid element
        //---------------------------------------------------------------------------
        tic(0);

        // Collect the amounts of elements in the solid and fluid species
        for(Index icell = 0; icell < num_cells; ++icell)
        {
            bf.row(icell) = states[icell].elementAmountsInSpecies(ifs);
            bs.row(icell) = states[icell].elementAmountsInSpecies(iss);
        }

        // Left boundary condition cell
        Index icell_bc = 0;

        // Get porosity of the left boundary cell
        const auto phi_bc = properties[icell_bc].fluidVolume().val /  properties[icell_bc].volume().val;

        // Ensure the result of each fluid element transport calculation can be saved
        result.transport_of_element.resize(num_elements);

        // Transport the elements in the fluid species
        for(Index ielement = 0; ielement < num_elements; ++ielement)
        {
            // Scale BC with a porosity of the boundary cell
            transport_solver.setBoundaryValue(phi_bc * bbc[ielement]);
            transport_solver.step(bf.col(ielement));

            // Save the result of this element transport calculation.
            result.transport_of_element[ielement] = transport_solver.result();
        }

        // Sum the amounts of elements distributed among fluid and solid species
        b.noalias() = bf + bs;

        toc(0, result.timing.transport);

        //---------------------------------------------------------------------------
        // Step 2: Perform a time step equilibrium calculation for each cell
        //---------------------------------------------------------------------------
        tic(1);

        if(options.use_smart_equilibrium_solver)
        {
            // Ensure the result of each cell's smart equilibrium calculation can be saved
            result.smart_equilibrium_at_cell.resize(num_cells);

            for(Index icell = 0; icell < num_cells; ++icell)
            {
                const auto T = states[icell].temperature();
                const auto P = states[icell].pressure();

                // Solve with a smart equilibrium solver
                smart_equilibrium_solver.solve(states[icell], T, P, b.row(icell));
                //if(!smart_equilibrium_solver.result().estimate.accepted) std::cout << " on (step, cell) : (" << steps << ", " << icell << ")" << std::endl;

                // Update chemical properties of the field
                properties[icell] = smart_equilibrium_solver.properties();

                // Save the result of this cell's smart equilibrium calculation
                result.smart_equilibrium_at_cell[icell] = smart_equilibrium_solver.result();
            }
        }
        else
        {
            // Ensure the result of each cell's equilibrium calculation can be saved.
            result.equilibrium_at_cell.resize(num_cells);

            for(Index icell = 0; icell < num_cells; ++icell)
            {
                const auto T = states[icell].temperature();
                const auto P = states[icell].pressure();

                // Solve with a conventional equilibrium solver
                equilibrium_solver.solve(states[icell], T, P, b.row(icell));

                // Update chemical properties of the field
                properties[icell] = equilibrium_solver.properties();

                // Save the result of this cell's smart equilibrium calculation.
                result.equilibrium_at_cell[icell] = equilibrium_solver.result();
            }
        }

        toc(1, result.timing.equilibrium);

        // Update the output files with the chemical state of every cell
        for(Index icell = 0; icell < num_cells; ++icell)
            for(auto output : outputs) {
                output.update(states[icell], icell);
                //output.update(states[icell], properties[icell], icell);
            }

        // Output chemical states in the output files
        for(auto output : outputs)
            output.close();

        // Increment the current number of reactive transport steps
        ++steps;

        return rt_result;
    }

    /// Perform one time step of a reactive transport calculation.
    auto stepKinetics(ChemicalField& field) -> ReactiveTransportResult
    {
        // The result of the reactive transport step
        ReactiveTransportResult rt_result = {};

        // Auxiliary variables
        const auto& mesh = transport_solver.mesh();
        const auto& num_elements = system.numElements();
        const auto& num_cells = mesh.numCells();
        const auto& ifs = system.indicesFluidSpecies();
        const auto& iss = system.indicesSolidSpecies();
        const auto& dt = transport_solver.timeStep();
        const auto& t_start = steps * dt;

        auto& states = field.states();
        auto& properties = field.properties();

        // Open the the file for outputting chemical states
        for(auto output : outputs)
        {
            output.suffix("-" + std::to_string(steps));
            output.open();
        }

        //---------------------------------------------------------------------------
        // Step 1: Perform a time step transport calculation for each fluid element
        //---------------------------------------------------------------------------
        tic(0);

        // Collect the amounts of elements in the solid and fluid species
        for(Index icell = 0; icell < num_cells; ++icell)
        {
            bf.row(icell) = states[icell].elementAmountsInSpecies(ifs);
            bs.row(icell) = states[icell].elementAmountsInSpecies(iss);
        }

        // Left boundary condition cell
        Index icell_bc = 0;
        // Get porosity of the left boundary cell
        const auto phi_bc = properties[icell_bc].fluidVolume().val /  properties[icell_bc].volume().val;

        // Ensure the result of each fluid element transport calculation can be saved
        result.transport_of_element.resize(num_elements);

        // Transport the elements in the fluid species
        for(Index ielement = 0; ielement < num_elements; ++ielement)
        {
            // Scale BC with a porosity of the boundary cell
            transport_solver.setBoundaryValue(phi_bc * bbc[ielement]);
            transport_solver.step(bf.col(ielement));

            // Save the result of this element transport calculation.
            result.transport_of_element[ielement] = transport_solver.result();
        }

        // Sum the amounts of elements distributed among fluid and solid species
        b.noalias() = bf + bs;

        toc(0, result.timing.transport);

        //---------------------------------------------------------------------------
        // Step 2: Perform a time step equilibrium calculation for each cell
        //---------------------------------------------------------------------------
        tic(1);

        if(options.use_smart_kinetic_solver)
        {
            // Ensure the result of each cell's smart kinetic calculation can be saved
            result.smart_kinetics_at_cell.resize(num_cells);
            result.smart_equilibrium_at_cell.resize(num_cells);
            result.equilibrium_at_cell.resize(num_cells);

            for(Index icell = 0; icell < num_cells; ++icell)
            {
                const auto T = states[icell].temperature();
                const auto P = states[icell].pressure();

                // Solve with a smart kinetic solver
                smart_kinetic_solver.solve(states[icell], t_start, dt, b.row(icell));

                // Update chemical properties of the field
                properties[icell] = smart_kinetic_solver.properties();
                /*
                std::cout << "icell        : " << icell << std::endl;
                if(steps > 200 && icell >- 0 && icell < 20) {
                    std::cout << " step, cell : " << steps << " " << icell;
                    std::cout << " Calcite : " << states[icell].speciesAmount("Calcite") ;
                    std::cout << " Dolomite : " << states[icell].speciesAmount("Dolomite") << std::endl;
                    getchar();
                }
                */
                // Save the result of this cell's smart equilibrium calculation.
                result.smart_kinetics_at_cell[icell] = smart_kinetic_solver.result();
                if(options.use_smart_equilibrium_solver)    result.smart_equilibrium_at_cell[icell] = smart_kinetic_solver.result().smart_equilibrium;
                else                                        result.equilibrium_at_cell[icell] = smart_kinetic_solver.result().equilibrium;
                /*
                // * TODO: debugging code to remove

                std::cout << " - solve                     : " << result.smart_kinetics_at_cell[icell].timing.solve << std::endl;
                std::cout << "   - learn : " << result.smart_kinetics_at_cell[icell].timing.learn << " (" << result.smart_kinetics_at_cell[icell].timing.learn / result.smart_kinetics_at_cell[icell].timing.solve * 100 << " %)" << std::endl;
                std::cout << "     - store                 : " << result.smart_kinetics_at_cell[icell].timing.learn_storage << " (" << result.smart_kinetics_at_cell[icell].timing.learn_storage / result.smart_kinetics_at_cell[icell].timing.solve * 100 << " %)" << std::endl;
                std::cout << "     - integrate             : " << result.smart_kinetics_at_cell[icell].timing.learn_integration << " (" << result.smart_kinetics_at_cell[icell].timing.learn_storage / result.smart_kinetics_at_cell[icell].timing.solve * 100 << " %)" << std::endl;
                std::cout << "       - chemical properties : " << result.smart_kinetics_at_cell[icell].timing.learn_chemical_properties << " (" << result.smart_kinetics_at_cell[icell].timing.learn_chemical_properties / result.smart_kinetics_at_cell[icell].timing.solve * 100 << " %)" << std::endl;
                std::cout << "       - integrate_reaction_rates      : " << result.smart_kinetics_at_cell[icell].timing.learn_reaction_rates << " (" << result.smart_kinetics_at_cell[icell].timing.learn_reaction_rates / result.smart_kinetics_at_cell[icell].timing.solve * 100 << " %)" << std::endl;
                std::cout << "       - sensitivity         : " << result.smart_kinetics_at_cell[icell].timing.learn_sensitivity << " (" << result.smart_kinetics_at_cell[icell].timing.learn_sensitivity / result.smart_kinetics_at_cell[icell].timing.solve * 100 << " %)" << std::endl;
                std::cout << "       - equilibration       : " << result.smart_kinetics_at_cell[icell].timing.learn_equilibration << " (" << result.smart_kinetics_at_cell[icell].timing.learn_equilibration / result.smart_kinetics_at_cell[icell].timing.solve * 100 << " %)" << std::endl;
                std::cout << "   - estimate              : " << result.smart_kinetics_at_cell[icell].timing.estimate << " (" << result.smart_kinetics_at_cell[icell].timing.estimate / result.smart_kinetics_at_cell[icell].timing.solve * 100 << " %)" << std::endl;
                std::cout << "     - estimate            : " << result.smart_kinetics_at_cell[icell].timing.estimate_search << " (" << result.smart_kinetics_at_cell[icell].timing.estimate_search / result.smart_kinetics_at_cell[icell].timing.solve * 100 << " %)" << std::endl;
                std::cout << "     - acceptance          : " << result.smart_kinetics_at_cell[icell].timing.estimate_acceptance << " (" << result.smart_kinetics_at_cell[icell].timing.estimate_acceptance / result.smart_kinetics_at_cell[icell].timing.solve * 100 << " %)" << std::endl;
                */
            }
        }
        else
        {
            // Ensure the result of each cell's kinetic calculation can be saved.
            result.kinetics_at_cell.resize(num_cells);
            result.smart_equilibrium_at_cell.resize(num_cells);
            result.equilibrium_at_cell.resize(num_cells);

            for(Index icell = 0; icell < num_cells; ++icell)
            {
                const auto T = states[icell].temperature();
                const auto P = states[icell].pressure();

                // Solve with a conventional kinetic solver
                kinetic_solver.solve(states[icell], t_start, dt, b.row(icell));

                // Update chemical properties of the field
                properties[icell] = kinetic_solver.properties();

                // Save the result of this cell's smart equilibrium calculation.
                result.kinetics_at_cell[icell] = kinetic_solver.result();
                if(options.use_smart_equilibrium_solver)    result.smart_equilibrium_at_cell[icell] = kinetic_solver.result().smart_equilibrium;
                else                                        result.equilibrium_at_cell[icell] = kinetic_solver.result().equilibrium;

                /*
                // * TODO: debugging code to remove
                std::cout << " - solve                 : " << result.kinetics_at_cell[icell].timing.solve << std::endl;
                std::cout << "   - chemical properties : " << result.kinetics_at_cell[icell].timing.chemical_properties << " (" << result.kinetics_at_cell[icell].timing.chemical_properties / result.kinetics_at_cell[icell].timing.solve * 100 << " %)" << std::endl;
                std::cout << "   - integrate_reaction_rates      : " << result.kinetics_at_cell[icell].timing.integrate_reaction_rates << " (" << result.kinetics_at_cell[icell].timing.integrate_reaction_rates / result.kinetics_at_cell[icell].timing.solve * 100 << " %)" << std::endl;
                std::cout << "   - sensitivity         : " << result.kinetics_at_cell[icell].timing.sensitivity << " (" << result.kinetics_at_cell[icell].timing.sensitivity / result.kinetics_at_cell[icell].timing.solve * 100 << " %)" << std::endl;
                std::cout << "   - equilibration       : " << result.kinetics_at_cell[icell].timing.equilibrate << " (" << result.kinetics_at_cell[icell].timing.equilibrate / result.kinetics_at_cell[icell].timing.solve * 100 << " %)" << std::endl;
                if(options.use_smart_equilibrium_solver)
                {
                    std::cout << "     - learning       : " << result.kinetics_at_cell[icell].smart_equilibrium.timing.learn << " (" << result.kinetics_at_cell[icell].smart_equilibrium.timing.learn / result.kinetics_at_cell[icell].timing.solve * 100 << " %)" << std::endl;
                    std::cout << "     - estimation     : " << result.kinetics_at_cell[icell].smart_equilibrium.timing.estimate << " (" << result.kinetics_at_cell[icell].smart_equilibrium.timing.estimate / result.kinetics_at_cell[icell].timing.solve * 100 << " %)" << std::endl;
                }
                getchar();
                */
            }
        }

        toc(1, result.timing.kinetics);

        /*
         * // TODO: remove, left from debugging
        double solve_time = 0.0;
        double initialize = 0.0;
        double integrate = 0.0;
        double chemical_properties = 0.0;
        double reaction_rates = 0.0;
        double sensitivity = 0.0;
        double equilibration = 0.0;
        double estimation = 0.0;
        double estimation_search = 0.0;
        double learning = 0.0;
        double kin_learn = 0.0;
        double kin_estimate = 0.0;
        double kin_estimate_search = 0.0;
        double kin_estimate_acceptance = 0.0;

        auto sum = [&](const ReactiveTransportResult& res){
            for(unsigned int i=0; i < num_cells; ++i)
            {
                if(options.use_smart_kinetic_solver){
                    initialize += res.smart_kinetics_at_cell.at(i).timing.initialize;
                    solve_time += res.smart_kinetics_at_cell.at(i).timing.solve;
                    integrate += res.smart_kinetics_at_cell.at(i).timing.learn_integration;

                    chemical_properties += res.smart_kinetics_at_cell.at(i).timing.learn_chemical_properties;
                    reaction_rates += res.smart_kinetics_at_cell.at(i).timing.learn_reaction_rates;
                    sensitivity += res.smart_kinetics_at_cell.at(i).timing.learn_sensitivity;
                    equilibration += res.smart_kinetics_at_cell.at(i).timing.learn_equilibration;
                    kin_learn += res.smart_kinetics_at_cell.at(i).timing.learn;
                    kin_estimate += res.smart_kinetics_at_cell.at(i).timing.estimate;
                    kin_estimate_search += res.smart_kinetics_at_cell.at(i).timing.estimate_search;
                    kin_estimate_acceptance += res.smart_kinetics_at_cell.at(i).timing.estimate_acceptance;

                }else{
                    initialize += res.kinetics_at_cell.at(i).timing.initialize;
                    solve_time += res.kinetics_at_cell.at(i).timing.solve;
                    integrate += res.kinetics_at_cell.at(i).timing.integrate;

                    chemical_properties += res.kinetics_at_cell.at(i).timing.integrate_chemical_properties;
                    reaction_rates += res.kinetics_at_cell.at(i).timing.integrate_reaction_rates;
                    sensitivity += res.kinetics_at_cell.at(i).timing.integrate_sensitivity;
                    equilibration += res.kinetics_at_cell.at(i).timing.integrate_equilibration;
                }
                if(options.use_smart_equilibrium_solver) {

                    learning += res.smart_equilibrium_at_cell.at(i).timing.learn;
                    estimation += res.smart_equilibrium_at_cell.at(i).timing.estimate;
                    estimation_search += res.smart_equilibrium_at_cell.at(i).timing.estimate_search;
                }
            }
        };
        sum(result);


        std::cout << "transport : " << result.timing.transport << std::endl;
        std::cout << "kinetics  : " << result.timing.kinetics << std::endl;
        std::cout << " - solve                   : " << solve_time << std::endl;
        std::cout << "   - learn                 : " << kin_learn << " (" << kin_learn / solve_time * 100 << " %)" << std::endl;
        //std::cout << "     - initialize            : " << initialize << " (" << initialize / solve_time * 100 << " %)" << std::endl;
        std::cout << "     - integrate             : " << integrate << " (" << integrate / solve_time * 100 << " %)" << std::endl;
        std::cout << "       - chemical properties : " << chemical_properties << " (" << chemical_properties / solve_time * 100 << " %)" << std::endl;
        std::cout << "       - sensitivity         : " << sensitivity << " (" << sensitivity / solve_time * 100 << " %)" << std::endl;
        std::cout << "       - integrate_reaction_rates      : " << integrate_reaction_rates << " (" << integrate_reaction_rates / solve_time * 100 << " %)" << std::endl;
        std::cout << "       - equilibration       : " << equilibration << " (" << equilibration / solve_time * 100 << " %)" << std::endl;
        if (options.use_smart_kinetic_solver) {
            std::cout << "   - estimate                 : " << kin_estimate << " (" << kin_estimate / solve_time * 100 << " %)" << std::endl;
            std::cout << "     - search                 : " << kin_estimate_search << " (" << kin_estimate_search / solve_time * 100 << " %)" << std::endl;
            std::cout << "     - acceptance             : " << kin_estimate_acceptance << " (" << kin_estimate_acceptance / solve_time * 100 << " %)" << std::endl;
        }
        int number_learnings;
        auto sum_learnings = [&](const ReactiveTransportResult& res){
            unsigned int number_learnings = 0;
            for(unsigned int i=0; i < num_cells; ++i)
                number_learnings += static_cast<int>(res.smart_kinetics_at_cell.at(i).estimate.accepted);
            return number_learnings;
        };
        std::cout << "------------------------------------------------" << std::endl;
        //getchar();
        */
        /*
        if(options.use_smart_kinetic_solver) {
            auto sum_learnings = [&](const ReactiveTransportResult &res) {
                unsigned int number_learnings = 0;
                for (unsigned int i = 0; i < num_cells; ++i)
                    number_learnings += static_cast<int>(res.smart_kinetics_at_cell.at(i).estimate.accepted);
                return number_learnings;
            };
            std::cout << "number of learnings on step " << steps << " : " << num_cells - sum_learnings(result)
                      << " out of " << num_cells << std::endl;
            if (!(steps % 10)) getchar();
        }
         */

        /*
        if(steps == 49) {
            //TODO: code for debugging
            std::cout << "transport : " << result.timing.transport << std::endl;
            std::cout << "kinetics  : " << result.timing.kinetics << std::endl;
            std::cout << " - solve                   : " << solve_time << std::endl;
            std::cout << "   - initialize            : " << initialize << " (" << initialize / solve_time * 100 << " %)"
                      << std::endl;
            std::cout << "   - integrate             : " << integrate << " (" << integrate / solve_time * 100 << " %)"
                      << std::endl;
            std::cout << "     - chemical properties : " << chemical_properties << " ("
                      << chemical_properties / solve_time * 100 << " %)" << std::endl;
            std::cout << "     - sensitivity : " << sensitivity << " ("
                      << sensitivity / solve_time * 100 << " %)" << std::endl;
           std::cout << "     - integrate_reaction_rates      : " << integrate_reaction_rates << " (" << integrate_reaction_rates / solve_time * 100
                      << " %)" << std::endl;
            std::cout << "     - equilibration       : " << equilibration << " (" << equilibration / solve_time * 100
                      << " %)" << std::endl;
            if (options.use_smart_equilibrium_solver) {
                std::cout << "       - learning          : " << learning << " (" << learning / solve_time * 100 << " %)"
                          << std::endl;
                std::cout << "       - estimation        : " << estimation << " (" << estimation / solve_time * 100
                          << " %)" << std::endl;
                std::cout << "         - estimation_search : " << estimation_search << " ("
                          << estimation_search / solve_time * 100 << " %)" << std::endl;
            }
            std::cout << "------------------------------------------------" << std::endl;
            getchar();
        }
        */

        // Update the output files with the chemical state of every cell
        for(Index icell = 0; icell < num_cells; ++icell)
            for(auto output : outputs){
                output.update(states[icell], icell);
                //output.update(states[icell], properties[icell], icell);
            }

        // Output chemical states in the output files
        for(auto output : outputs)
            output.close();

        // Increment the current number of reactive transport steps
        ++steps;

         return rt_result;
    }
};

ReactiveTransportSolver::ReactiveTransportSolver(const ChemicalSystem& system)
: pimpl(new Impl(system))
{
}

ReactiveTransportSolver::ReactiveTransportSolver(const ChemicalSystem& system, const ReactionSystem& reactions, const Partition& partition)
: pimpl(new Impl(system, reactions, partition))
{
}

ReactiveTransportSolver::ReactiveTransportSolver(const ReactiveTransportSolver& other)
: pimpl(new Impl(*other.pimpl))
{
}

ReactiveTransportSolver::~ReactiveTransportSolver()
{}

auto ReactiveTransportSolver::operator=(ReactiveTransportSolver other) -> ReactiveTransportSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto ReactiveTransportSolver::setOptions(const ReactiveTransportOptions& options) -> void
{
    pimpl->setOptions(options);
}

auto ReactiveTransportSolver::setMesh(const Mesh& mesh) -> void
{
    pimpl->setMesh(mesh);
}

auto ReactiveTransportSolver::setVelocity(double val) -> void
{
    pimpl->setVelocity(val);
}

auto ReactiveTransportSolver::setDiffusionCoeff(double val) -> void
{
    pimpl->setDiffusionCoeff(val);
}

auto ReactiveTransportSolver::setBoundaryState(const ChemicalState& state) -> void
{
    pimpl->setBoundaryState(state);
}

auto ReactiveTransportSolver::setTimeStep(double val) -> void
{
    pimpl->setTimeStep(val);
}

auto ReactiveTransportSolver::output() -> ChemicalOutput
{
    return pimpl->output();
}

auto ReactiveTransportSolver::initialize() -> void
{
    pimpl->initialize();
}

auto ReactiveTransportSolver::step(ChemicalField& field) -> ReactiveTransportResult
{
    return pimpl->step(field);
}

auto ReactiveTransportSolver::stepKinetics(ChemicalField& field) -> ReactiveTransportResult
{
    return pimpl->stepKinetics(field);
}

auto ReactiveTransportSolver::result() const -> const ReactiveTransportResult&
{
    return pimpl->result;
}

auto ReactiveTransportSolver::system() const -> const ChemicalSystem&
{
    return pimpl->system;
}

auto ReactiveTransportSolver::timeStep() const -> double
{
    return pimpl->transport_solver.timeStep();
}

} // namespace Reaktoro
