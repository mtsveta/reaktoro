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

#include "SmartKineticSolver.hpp"

// C++ includes
#include <functional>
using namespace std::placeholders;

// Reaktoro includes
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/Exception.hpp>
//#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Common/Profiling.hpp>
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
#include <Reaktoro/Equilibrium/SmartEquilibriumSolver.hpp>
#include <Reaktoro/Kinetics/SmartKineticOptions.hpp>
#include <Reaktoro/Kinetics/SmartKineticResult.hpp>
#include <Reaktoro/Kinetics/KineticState.hpp>

namespace Reaktoro {

/// A class used to store the node of tree for smart equilibrium calculations.
struct TreeNode{

    KineticState state;
    ChemicalState chemical_state;
    ChemicalProperties properties;
    EquilibriumSensitivity sensitivity;
    ChemicalVector rates;

    TreeNode(KineticState state, ChemicalState chemical_state, ChemicalProperties properties,
             EquilibriumSensitivity sensitivity, ChemicalVector rates)
            : state(state), chemical_state(chemical_state), properties(properties), sensitivity(sensitivity), rates(rates){}

    TreeNode& operator=(const TreeNode& other){

        this->state = other.state;
        this->chemical_state = other.chemical_state;
        this->properties = other.properties;
        this->sensitivity = other.sensitivity;
        this->rates = other.rates;

        return *this;
    }
};

struct SmartKineticSolver::Impl
{
    /// The kinetically-controlled chemical reactions
    ReactionSystem reactions;

    /// The chemical system instance
    ChemicalSystem system;

    /// The partition of the species in the chemical system
    Partition partition;

    /// The options of the kinetic solver
    SmartKineticOptions options;

    /// The result of the smart kinetic solver
    SmartKineticResult result;

    // The solver for solving the equilibrium equations using classical approach
    EquilibriumSolver equilibrium;

    /// The solver for solving the equilibrium equations using smart on-demand learning algorithm
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

    /// The number of elements in the equilibrium partition
    Index Ee;

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

     /// The tree used to save the calculated equilibrium states and respective sensitivities
    std::list<TreeNode> tree;

    Impl()
    {}

    Impl(const ReactionSystem& reactions)
    : reactions(reactions), system(reactions.system()), equilibrium(system), smart_equilibrium(system)
    {
        setPartition(Partition(system));
    }

    /*
    auto summarize() -> void {

        std::cout << "***********************************************************************************" << std::endl;
        std::cout << "***********************************************************************************" << std::endl;
        std::cout << "total time of SmartKineticSolver : " << stats.time_total << std::endl;
        std::cout << " - time kinetics             : " << stats.kinetics_time_total <<  " (" << stats.kinetics_time_total / stats.time_total * 100 << " %)" << std::endl;
        std::cout << "   - time estimaton          : " << stats.estimation_time_total <<  " (" << stats.estimation_time_total / stats.kinetics_time_total * 100 << " %)" << std::endl;
        std::cout << "     - time acceptance       : " << stats.confidence_time_total <<  " (" << stats.confidence_time_total / stats.estimation_time_total * 100 << " %)" << std::endl;
        std::cout << "     - time search           : " << stats.search_time_total     <<  " (" << stats.search_time_total / stats.estimation_time_total * 100 << " %)" << std::endl;
        std::cout << "   - time integration        : " << stats.integraion_time_total <<  " (" << stats.integraion_time_total / stats.kinetics_time_total * 100 << " %)" << std::endl;
        std::cout << "     - incr time prop. eval. : " << stats.prop_eval_time_total << " (" << stats.prop_eval_time_total / stats.integraion_time_total * 100 << " %)" << std::endl;
        std::cout << "     - incr time eq calls     : " << stats.eq_calls_time_total <<  " (" << stats.eq_calls_time_total / stats.integraion_time_total * 100 << " %)" << std::endl;
        if (options.learning.use_smart_equilibrium_solver) std::cout << "       - incr time eq search  : " << stats.eq_search_time_total <<  " (" << stats.eq_search_time / stats.integraion_time_total * 100 << " %)" << std::endl;
        std::cout << " - time equilibration        : " << stats.equilibration_time_total <<  " (" << stats.equilibration_time_total / stats.time_total * 100 << " %)" << std::endl;
        std::cout << " - # of integration cell     : " << stats.conv_counter << " out of " << stats.conv_counter + stats.smart_counter
                  << " (" << float(stats.conv_counter) / float(stats.conv_counter + stats.smart_counter) * 100 << ") % " << std::endl;
        std::cout << "total time of SmartKineticSolver (without prop. evals) : " << stats.time_total - stats.prop_eval_time_total << std::endl;
        // <<  " -> speed-up of " << stats.time_total / (stats.time_total - stats.prop_eval_time_total) << std::endl;
        std::cout << "total time of SmartKineticSolver (without prop. evals and search) : "
                  << stats.time_total - stats.prop_eval_time_total - stats.search_time_total << std::endl;
        //<<  " -> speed-up of " << stats.time_total / (stats.time_total - stats.prop_eval_time_total - stats.search_time_total) << std::endl;
        if (options.learning.use_smart_equilibrium_solver) {
            std::cout << "total time of SmartKineticSolver (without prop. evals, search, and eq search) : "
                      << stats.time_total - stats.prop_eval_time_total - stats.search_time_total - stats.eq_search_time_total << std::endl;
            //<<  " -> speed-up of " << stats.time_total / (stats.time_total - stats.prop_eval_time_total - stats.search_time_total - stats.eq_search_time_total) << std::endl;
        }
        std::cout << "***********************************************************************************" << std::endl;
        std::cout << "***********************************************************************************" << std::endl;

    }
    auto console() -> void {
        //std::cout << "\nsmart_counter : " << stats.smart_counter << std::endl;
        //std::cout << "conv_counter  : " << stats.conv_counter << std::endl;

        double total_time_kinetics = stats.estimation_time + stats.confidence_time + stats.integraion_time;
        double total_time = stats.estimation_time + stats.confidence_time + stats.integraion_time + stats.equilibration_time;

        std::cout << "incr time total              : " << total_time << std::endl;
        std::cout << " - incr time kinetics        : " << total_time_kinetics <<  " (" << total_time_kinetics / total_time * 100 << " %)" << std::endl;
        std::cout << "   - incr time estimaton     : " << stats.estimation_time <<  " (" << stats.estimation_time / total_time_kinetics * 100 << " %)" << std::endl;
        std::cout << "     - incr time acceptance  : " << stats.confidence_time <<  " (" << stats.confidence_time / stats.estimation_time * 100 << " %)" << std::endl;
        std::cout << "     - incr time search      : " << stats.search_time     <<  " (" << stats.search_time / stats.estimation_time * 100 << " %)" << std::endl;
        std::cout << "   - incr time integration   : " << stats.integraion_time <<  " (" << stats.integraion_time / total_time_kinetics * 100 << " %)" << std::endl;
        std::cout << "     - incr time prop. eval. : " << stats.prop_eval_time << " (" << stats.prop_eval_time / stats.integraion_time * 100 << " %)" << std::endl;
        std::cout << "     - incr time eq calls    : " << stats.eq_calls_time <<  " (" << stats.eq_calls_time / stats.integraion_time * 100 << " %)" << std::endl;
        if (options.learning.use_smart_equilibrium_solver)  std::cout << "       - incr time eq search : " << stats.eq_search_time <<  " (" << stats.eq_search_time / total_time_kinetics * 100 << " %)" << std::endl;
        std::cout << " - incr time equilibration   : " << stats.equilibration_time <<  " (" << stats.equilibration_time / total_time * 100 << " %)" << std::endl;
        std::cout << " - # of integration cell     : " << stats.conv_counter << " out of " << stats.conv_counter + stats.smart_counter
                  << " (" << float(stats.conv_counter) / float(stats.conv_counter + stats.smart_counter) * 100 << ") % " << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
    }
    auto refresh() -> void {

        // Add the smart counter to the vector of smart counters on each step
        stats.smart_numbers.emplace_back(stats.smart_counter);

        stats.estimation_time_total += stats.estimation_time;
        stats.confidence_time_total += stats.confidence_time;
        stats.search_time_total += stats.search_time;
        stats.integraion_time_total += stats.integraion_time;
        stats.prop_eval_time_total     += stats.prop_eval_time;

        stats.eq_calls_time_total  += stats.eq_calls_time;
        stats.eq_search_time_total += stats.eq_search_time;

        stats.equilibration_time_total += stats.equilibration_time;
        stats.kinetics_time_total += stats.integraion_time + stats.estimation_time;

        stats.time_total = stats.kinetics_time_total + stats.equilibration_time_total;

        // Update the counter of smartly predicted cells per step
        //stats.smart_counter = 0;

        // Update the counter of smartly predicted cells per step
        //stats.conv_counter = 0;

        // Update the time measured for the current step
        stats.estimation_time = 0.0;
        stats.confidence_time = 0.0;
        stats.integraion_time = 0.0;
        stats.search_time = 0.0;
        stats.equilibration_time = 0.0;
        stats.prop_eval_time     = 0.0;
        stats.kinetics_time      = 0.0;

        stats.eq_calls_time = 0.0;
        stats.eq_search_time = 0.0;

        // Empty the tree for the next step
        //tree.clear();
    }
    */

    auto setOptions(const SmartKineticOptions& options) -> void
    {
        // Initialise the options of the kinetic solver
        this->options = options;

        // Set options for the equilibrium calculations
        equilibrium.setOptions(options.equilibrium);
        smart_equilibrium.setOptions(options.smart_equilibrium);

    }

    auto setPartition(const Partition& partition_) -> void
    {
        // Initialise the partition member
        partition = partition_;

        // Set the partition of the equilibrium solver
        smart_equilibrium.setPartition(partition);
        equilibrium.setPartition(partition);

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

        {
            /*
            * TODO: remove when debugging is finished
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
            */
        }
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
        benk.head(Ee) = be;
        benk.tail(Nk) = nk;

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
        ODEOptions options_ode = options.learning.ode;

        // Set the ODE problem and initialize the ODE solver
        ode.setProblem(problem);
        ode.setOptions(options_ode);
        ode.initialize(tstart, benk);

        // Set the options of the equilibrium solver
        smart_equilibrium.setOptions(options.learning.smart_equilibrium);
        equilibrium.setOptions(options.learning.equilibrium);

    }

    auto setElementsAmountsPerCell(const ChemicalState& state, VectorConstRef b) -> void
    {
        // Fetch the amounts of spacies from the chemical state
        const auto& n = state.speciesAmounts();
        nk = n(iks);
        ne = n(ies);

        // Update amounts of elements
        be = b - Ak * nk;

    }

    auto learn(ChemicalState& state, double t, double dt) -> void
    {
        SmartKineticResult res = {};

        /// The initial state, rhf, Jacobian, and sensitivity matrix of
        /// combined vector of elemental molar abundance and composition of kinetic species [be nk]
        Vector benk0, benk_f, benk_J, benk_S;

        // Initialize the initial state benk0
        benk0.resize(Ee + Nk);
        benk0 = benk;

        // Allocate the memory for the rhs, Jacobian, and sensitivity matrix
        benk_f.resize(Ee + Nk);
        benk_J.resize(Ee + Nk, Ee + Nk);
        benk_S.resize(Ee + Nk, Ee + Nk);

        // Initialize sensitivity matrix by the identity matrix
        benk_S.setIdentity();

        // Initialize the kinetics state with the data at times t0 and tk
        KineticState kin_state;
        kin_state.n0 = benk0;
        kin_state.nk = benk;
        kin_state.tk = t;

        // Reconstract `benk` by the conventional numerical integration
        timeit(ode.solve(t, dt, benk, benk_f, benk_S), result.timing.learn_integrate +=);

        // Update the composition of the kinetic species
        nk = benk.tail(Nk);
        state.setSpeciesAmounts(nk, iks);

        // Save the sensitivity values, the result time, and the obtain species' amount
        kin_state.tk1 = t;
        kin_state.nk1 = benk;
        kin_state.dndn0 = benk_S;
        kin_state.dndt = benk_f;
        kin_state.dfdn = benk_J;

        // Evaluate chemical properties, reaction rate, and sensitivities
        timeit(properties = state.properties(), result.timing.learn_chemical_properties +=);
        timeit(r = reactions.rates(properties), result.timing.learn_reaction_rates+=);
        timeit(sensitivity = options.use_smart_equilibrium_solver ? smart_equilibrium.sensitivity() : equilibrium.sensitivity(),
                result.timing.learn_sensitivity +=);

        // Save the kinetic state to the tree of learned states
        timeit(tree.emplace_back(TreeNode(kin_state, state, properties, sensitivity, r)), result.timing.learn_store+=);
    }
    /*
    auto stepSmart(ChemicalState& state, double t, double tfinal) -> double
    {
        // TODO: rewrite this method.
        /*
        //Vector benk0 = benk;
        TreeNode ref_node;
        Vector benk_est;
        benk_est.resize(Ee + Nk, 1);

        double dt = tfinal - t;

        Time confidence_time = time();
        auto confident = isInitialGuessConfident(benk0, t, ref_node); // considering n0 as the initial condition

        Time start = time();
        //if(false){
        if(confident) { // Try estimating using sensitivities

            // Algorithm:
            // the reference state: n, n0, S = dn/dn0, t, f = dn/dt
            // new initial values :  n0_tilde, t_tilde
            // the predicted state: n_tilde = n + S * (n0_tilde - n0) + f * (t_tilde - t)

            // Clarification:
            // benk0 is new initial condition (n0_tilde)
            // t     is new time (t_tilde)
            // ref_node.state.nk1   is already calculated by integration vector  (n)
            // ref_node.state.dndn0 is the sensitivity w.r.t. the initial condition (S)
            // ref_node.state.n0    is the initial condition of reference vector (n0)
            // ref_node.state.dndt  is the sensitivity w.r.t. the change in time (S)
            // ref_node.state.tk1   is the time, at which reference state is constructed (t)

            // Linear approximation of the new kinetics path value
            dt = ref_node.state.tk1 - ref_node.state.tk;
            benk = ref_node.state.nk1 + ref_node.state.dndn0 * (benk0 - ref_node.state.n0)
                   + ref_node.state.dndt * (t + dt - ref_node.state.tk1);

            //std::cout << "tk     : " << ref_node.state.tk << std::endl;
            //std::cout << "tk1    : " << ref_node.state.tk1 - ref_node.state.tk << std::endl;
            //std::cout << "tfinal : " << tfinal << std::endl;
            //std::cout << "benk (after estimation)  : " << tr(benk) << std::endl;
            //std::cout << "time for estimating      : " << elapsed(start) << std::endl;
            t += dt;
        }
        else { // Perform one ODE step integration with forward sensitivities

            // Initialize the kinetics state with the data at times t0 and tk
            KineticState kin_state;
            kin_state.n0 = benk0;
            kin_state.nk = benk;
            kin_state.tk = t;

            // Integrate ODE within the time interval [t, tfinal]
            ode.integrate(t, benk, tfinal, benk_S, benk_J, benk_f);

            // Increase the counter of the conventionally integrated cells
            stats.conv_counter += 1;

            // Save the sensitivity values, the result time, and the obtain species' amount
            kin_state.tk1 = t;
            kin_state.nk1 = benk;
            kin_state.dndn0 = benk_S;
            kin_state.dndt = benk_f;
            //state.dfdn = benk_J;

            // Save the kinetic state as the tree node
            TreeNode node;
            node.state = kin_state;
            tree.emplace_back(node);

            //std::cout << "tk    : " << state.tk << std::endl;
            //std::cout << "tk1   : " << state.tk1 << std::endl;
            //std::cout << "tfinal: " << tfinal << std::endl;
            //std::cout << "dt    : " << state.tk1 - state.tk << std::endl;
            //std::cout << "benk (after integration) : " << tr(benk) << std::endl;
            //std::cout << "time for integrating     : " << elapsed(start) << std::endl;//if (confident) std::cout << "error (estim - integr)   : " << (benk_est - benk).squaredNorm() << std::endl;
        }

        return t;
    }


    auto solveSmart(ChemicalState& state, double& t, double tfinal, Index step, Index cell, KineticResult& res) -> void
    {
        Time start = time();
        bool result = estimate(t, res);
        stats.estimation_time += elapsed(start);
        //res.smart.stats.time_estimate = stats.estimation_time;

        start = time();
        if (result)
            return;
        else
            learn(state, t, tfinal, step, cell, res);
        stats.integraion_time += elapsed(start);
        //res.smart.stats.time_estimate = stats.integraion_time;

    }
    */
    auto estimate(ChemicalState& state, double& t, VectorConstRef benk0) -> void
    {
        // Skip estimation if no previous full computation has been done
        if(tree.empty())
            return;

        // Relative and absolute tolerance parameters
        // TODO: should they be different scale?
        auto kinetics_reltol = options.reltol;
        auto kinetics_abstol = options.abstol;
        auto equilibrium_reltol = options.reltol;
        auto equilibrium_abstol = options.abstol;
        auto equilibrium_cutoff = options.cutoff;

        // Comparison function based on the Euclidean distance
        auto distancefn = [&](const TreeNode& a, const TreeNode& b)
        {
            const auto& nk_a = a.state.nk;
            const auto& nk_b = b.state.nk;

            double dist_to_a = std::pow((nk_a - benk0).squaredNorm(), 0.5);
            double dist_to_b = std::pow((nk_b - benk0).squaredNorm(), 0.5);

            return dist_to_a < dist_to_b;
        };

        //---------------------------------------------------------------------------------------
        // Step 1: Search for the reference element (closest to the new state input conditions)
        //---------------------------------------------------------------------------------------
        tic(0);

        // Find the reference element (nearest to the new state benk)
        auto it = std::min_element(tree.begin(), tree.end(), distancefn);

        toc(0, result.timing.estimate_search);


        //----------------------------------------------------------------------------
        // Step 2: Calculate predicted state with a first-order Taylor approximation
        //----------------------------------------------------------------------------
        tic(1);

        // Get all the data stored in the reference element
        const auto& benk0_ref = it->state.n0;
        const auto& benkk_ref = it->state.nk;
        const auto& benkk1_ref = it->state.nk1;
        const auto& tk_ref = it->state.tk;
        const auto& tk1_ref = it->state.tk1;
        const auto& dndn0_ref = it->state.dndn0;
        const auto& dndt_ref = it->state.dndt;
        const auto& state_ref = it->chemical_state;
        const auto& properties_ref = it->properties;
        const auto& sensitivity_ref = it->sensitivity;
        const auto& r_ref = it->rates;

        const auto& N = ies.size() + iks.size();
        const auto& be_ref = benkk_ref.head(Ee);
        const auto& n_ref = state_ref.speciesAmounts();
        const auto& ne_ref = n_ref(ies);
        const auto& nk_ref = n_ref(iks);

        // Algorithm:
        // the reference state: n, n0, S = dn/dn0, t, f = dn/dt
        // new initial values :  n0_tilde, t_tilde
        // the predicted state: n_tilde = n + S * (n0_tilde - n0) + f * (t_tilde - t)

        // Clarification:
        // benk0 is new initial condition (n0_tilde)
        // t     is new time (t_tilde)
        // ref_node.state.nk1   is already calculated by integration vector  (n)
        // ref_node.state.dndn0 is the sensitivity w.r.t. the initial condition (S)
        // ref_node.state.n0    is the initial condition of reference vector (n0)
        // ref_node.state.dndt  is the sensitivity w.r.t. the change in time (S)
        // ref_node.state.tk1   is the time, at which reference state is constructed (t)

        // Perform smart estimation of benk
        double dt = tk1_ref - tk_ref;
        benk = benkk1_ref + dndn0_ref * (benk0 - benk0_ref) + dndt_ref * (t - tk_ref);

        // Fetch the be and nk unknows from vector benk = [be; nk]
        VectorConstRef be = benk.head(Ee);
        VectorConstRef nk = benk.tail(Nk);

        // Update the composition of the kinetic species
        state.setSpeciesAmounts(nk, iks);

        // Get the sensitivity derivatives dln(a) / dn
        MatrixConstRef dlnadn_ref = properties_ref.lnActivities().ddn; // TODO: this line is assuming all species are equilibrium specie! get the rows and columns corresponding to equilibrium species
        VectorConstRef lna_ref = properties_ref.lnActivities().val;

        // Get the sensitivity derivatives w.r.t. the equilibrium species dln(a) / dne
        MatrixConstRef dlnadne_ref = dlnadn_ref(ies, ies);
        VectorConstRef lnae_ref = lna_ref(ies);

        properties = state.properties();
        r = reactions.rates(properties);

        VectorConstRef lna = properties.lnActivities().val;
        VectorConstRef lnae = lna(ies);

        //VectorConstRef diff_lnae = lnae_ref - lnae_new;
        //double dist_lnae = std::pow((lnae_ref - lnae_new).squaredNorm(), 2.0);
        //std::cout << "lnae_new: " << tr(lnae_new) << std::endl;
        //std::cout << "lnae_ref: " << tr(lnae_ref) << std::endl;
        //std::cout << "diff_lnae: " << tr(diff_lnae) << std::endl;
        //std::cout << "dist(lnae_new, lnae_ref): " << dist_lnae  << std::endl;

        // Define vectors of equilibirium species, delta(ne), delta(lnae)
        Vector ne, delta_ne, delta_lnae;

        delta_ne.noalias() = sensitivity_ref.dndb * (be - be_ref);    // delta(n) = dn/db * (b - b0)
        ne.noalias() = ne_ref + delta_ne;                             // n = n0 + delta(n)
        delta_lnae.noalias() = dlnadne_ref * delta_ne;                // delta(delta_lnae) = d (ln ae) / d(ne)

        //std::cout << "delta_ne: " << tr(delta_ne) << std::endl;
        //std::cout << "delta_lna: " << tr(delta_lna) << std::endl;
        //std::cout << "ne_new: " << tr(ne_new) << std::endl;

        toc(1, result.timing.estimate_mat_vec_mul);

        //----------------------------------------------
        // Step 3: Checking the acceptance criterion
        //----------------------------------------------
        tic(2);


        //const bool amount_check = ne.minCoeff() > equilibrium_cutoff;
        const auto& x = properties_ref.moleFractions().val;
        const auto& xe = x(ies);


        // Check the variations of equilibrium species
        bool equilibrium_amount_check = true;
        bool equilibrium_variation_check = true;

        const double fraction_tol = equilibrium_abstol * 1e-2;
        for(int i = 0; i < ne.size(); ++i){
            // If the fraction is too small, skip the variational check
            if(x[i] < fraction_tol) continue;

            ///std::cout << "abs(delta_lna[i]): " << std::abs(delta_lna[i]) << std::endl;
            //std::cout << "abstol + reltol * std::abs(lnae_ref[i])): " << abstol + reltol * std::abs(lnae_ref[i]) << std::endl;

            // Perform the variational check
            if(std::abs(delta_lnae[i]) > equilibrium_abstol + equilibrium_reltol * std::abs(lnae_ref[i])) {
                equilibrium_variation_check = false;
                result.estimate.failed_with_species = system.species(i).name();
                result.estimate.failed_with_amount = ne[i];
                result.estimate.failed_with_chemical_potential = lnae_ref[i] + delta_lnae[i]; // TODO: change to the chemical potential

                //result.estimate.accepted = false;    // variation test failed
                //i = ne.size();        // finish running through the loop
            }
            if(ne[i] < equilibrium_cutoff){
                equilibrium_amount_check = false;
                //std::cout << "neg amount in i: " << i << std::endl;
                //std::cout << "xe(i) : " << xe(i, 0) << " < fraction_tol : " << fraction_tol << std::endl;
                //std::cout << "ne[i] : " << ne[i] << std::endl;
                //std::cout << "Species: " << system.species(i).name() << std::endl;
                //getchar();
            }
        }
        /*
        for(int i = 0; i < ne_new.size(); ++i){
            // If the fraction is too small, skip the variational check
            if(x[i] < fraction_tol) continue;

            // Perform the variational check
            if(std::abs(diff_lnae[i]) > abstol + reltol * std::abs(lnae_ref[i])) {
                acceptance_eq = false;    // variation test failed
            }
        }
        */
        // Check the variations of kinetics' rates
        // ---------------------------------------------------------------------------------------
        //reltol = options.reltol;
        //abstol = options.abstol;

        // Initialize delta_n = [delta_ne; delta_nk]
        Vector delta_n;
        delta_n.resize(N);
        delta_n(ies) << delta_ne;
        delta_n(iks) << nk - nk_ref;
        //std::cout << "delta_ne: " << tr(delta_ne) << std::endl;
        //std::cout << "nk_new - nk_ref: " << nk_new - nk_ref << std::endl;
        //std::cout << "delta_n: " << tr(delta_n) << std::endl;

        // Compute delta_rk = drdn * delta_n
        Vector delta_rk;
        delta_rk.noalias() = r_ref.ddn * delta_n;
        //std::cout << "delta_rk.array().abs(): " << delta_rk.array().abs() << std::endl;
        //std::cout << "abstol + reltol * r_ref.val.array().abs(): " << abstol + reltol * r_ref.val.array().abs() << std::endl;
        const bool kinetics_variation_check
        = (delta_rk.array().abs() <= kinetics_abstol + kinetics_reltol * r_ref.val.array().abs()).all();

        //std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << tr(benk) << std::endl;
        //bool confident = true;
        if (equilibrium_variation_check && equilibrium_amount_check && kinetics_variation_check)
        {
            t += dt;
            return ;
        }

        result.estimate.accepted = true;
    }

    // learn(ChemicalState& state, double t, double dt)
    // auto estimate(ChemicalState& state, double& t, VectorConstRef benk0) -> SmartKineticResult&
    auto solve(ChemicalState& state, double t, double dt) -> SmartKineticResult
    {
        // Initialise the chemical kinetics solver
        initialize(state, t);

        tic(0);

        // Reset the result of the last smart equilibrium calculation
        result = {};

        // Perform a smart estimate for the chemical state
        timeit(estimate(state, t), result.timing.estimate =);

        // Perform a learning step if the smart prediction is not satisfactory
        if (!result.estimate.accepted)
            timeit(learn(state, t, dt), result.timing.learn =);


        toc(0, result.timing.solve);

        //std::cout << "benk: " << tr(benk) << std::endl;
        //std::cout << "t: " << t << std::endl;

        // Extract the `be` and `nk` entries of the vector `benk`
        be = benk.head(Ee);
        nk = benk.tail(Nk);

        // Update the composition of the kinetic species
        state.setSpeciesAmounts(nk, iks);

        return result;

        // Update the composition of the equilibrium species
        //equilibrium->solve(state, T, P, be);
        //Time start = time();
        //EquilibriumResult res;
        //if (options.learning.use_smart_equilibrium_solver) smart_equilibrium->solve(state, T, P, be);
        //else                                    equilibrium->solve(state, T, P, be);
        //stats.eq_calls_time += elapsed(start);

    }
    /*
    auto solveSmartPostprocess(ChemicalState& state, double& t, double tfinal, Index step, Index cell, KineticResult& res) -> void
    {
        Time start = time();
        bool success = estimateUsingResult(state, t, res, step, cell);
        stats.estimation_time += elapsed(start);
        //res.smart.stats.time_estimate = stats.estimation_time;

        start = time();
        if (success)
            return;
        else
            learn(state, t, tfinal, step, cell, res);
        stats.integraion_time += elapsed(start);
        //res.smart.stats.time_integrate = stats.integraion_time;
    }
    auto solvePath(ChemicalState& state, double t0, double t1, std::string units, KineticResult& result) -> void
    {
        t0 = units::convert(t0, units, "s");
        t1 = units::convert(t1, units, "s");

        initialize(state, t0);

        double t = t0;

        // Initialize the output of the equilibrium path calculation
        if(output) output.open();

        // Initialize the plots of the equilibrium path calculation
        for(auto& plot : plots) plot.open();

        while(t < t1)
        {
            // Update the output with current state
            if(output) output.update(state, t);

            // Update the plots with current state
            for(auto& plot : plots) plot.update(state, t);

            // Integration time step
            t = stepSmart(state, t, t1);
        }

        // Update the output with the final state
        if(output) output.update(state, t1);

        // Update the plots with the final state
        for(auto& plot : plots) plot.update(state, t1);


    }

    auto solveKineticsEquilibrium(ChemicalState &state, double t, double dt, VectorConstRef b, Index step, Index cell, KineticResult& res) -> void
    {

        postequilibrate = (cell == 0) ? true : postequilibrate;
        // Initialise the chemical kinetics solver
        initialize(state, t, b);

        // Integration time step
        solveSmartPostprocess(state, t, t + dt, step, cell, res);

        // Extract the `be` and `nk` entries of the vector `benk`
        be = benk.head(Ee);
        nk = benk.tail(Nk);

        // Update the composition of the kinetic species
        state.setSpeciesAmounts(nk, iks);
        Vector n_kin;
        n_kin.noalias() = state.speciesAmounts();

        // Update the composition of the equilibrium species
        Time start = time();
        //equilibrium->solve(state, T, P, be);
        if (postequilibrate) {
            // Update the composition of the equilibrium species
            if (options.learning.use_smart_equilibrium_solver) smart_equilibrium->solve(state, T, P, be);
            else equilibrium->solve(state, T, P, be);

            stats.equilibration_time += elapsed(start);
            //res.smart.stats.time_equilibrate = stats.equilibration_time;

            // Update the learning indices for the smart equilibrium

        }

        if (postequilibrate) {
            Vector n_eq;
            n_eq.noalias() = state.speciesAmounts();
            double rel_error =  std::sqrt((n_kin - n_eq).squaredNorm()) / std::sqrt((n_kin).squaredNorm());
            //std::cout << "(" << step << ", " << cell << ") : || n_kin - n_eq || / || n_eq ||  = " << rel_error << std::endl;
            //postequilibrate = (rel_error <= 1e-10) ? false : true;
        }
    }
    auto solveEquilibriumKinetics(ChemicalState& state, double t, double dt, VectorConstRef b, Index step, Index cell, KineticResult& res) -> void
    {

        // Initialise the temperature and pressure variables
        T = state.temperature();
        P = state.pressure();

        // Extract the composition of the equilibrium and kinetic species
        const auto& n = state.speciesAmounts();
        nk = n(iks);

        // Let kinetics know about the change in transport
        be = b - Ak * nk;

        // Set the options of the equilibrium solver
        if (options.learning.use_smart_equilibrium_solver)
            smart_equilibrium->setOptions(options.learning.smart_equilibrium);
        else
            equilibrium->setOptions(options.learning.equilibrium);

        Time start = time();

        // Update the composition of the equilibrium species
        //equilibrium->solve(state, T, P, be);

        if (options.learning.use_smart_equilibrium_solver)
            smart_equilibrium->solve(state, T, P, be);
        else
            equilibrium->solve(state, T, P, be);

        stats.equilibration_time += elapsed(start);
        //res.smart.stats.time_equilibrate = stats.equilibration_time;

        // Initialise the chemical kinetics solver
        initialize(state, t, b);

        // Integration time step
        // solveSmart(state, t, t + dt, step, cell, res);
        solveSmartPostprocess(state, t, t + dt, step, cell, res);

        // Extract the `be` and `nk` entries of the vector `benk`
        be = benk.head(Ee);
        nk = benk.tail(Nk);

        // Update the composition of the kinetic species
        state.setSpeciesAmounts(nk, iks);
    }
    */
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
        if(options.use_smart_equilibrium_solver) // using smart equilibrium solver
        {
            SmartEquilibriumResult res = {};

            // smart_equilibrium_result
            timeit(res += smart_equilibrium.solve(state, T, P, be), result.timing.learn_equilibration+=);

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
                timeit( res = smart_equilibrium.solve(state, T, P, be), result.timing.learn_equilibration+=);

            }

            // Assert the smart equilibrium calculation did not fail
            Assert(res.estimate.accepted || res.learning.gibbs_energy_minimization.optimum.succeeded,
                   "Could not calculate the rates of the species.",
                   "The smart equilibrium calculation failed.");

        }
        else  // using conventional equilibrium solver
        {
            EquilibriumResult res = {};

            timeit(res += equilibrium.solve(state, T, P, be), result.timing.learn_equilibration +=);

            result.equilibrium += res;

            // Check if the calculation failed, if so, use cold-start
            if (!res.optimum.succeeded) {
                state.setSpeciesAmounts(0.0);
                timeit(res = equilibrium.solve(state, T, P, be), result.timing.learn_equilibration +=);

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
        timeit(properties = state.properties(), result.timing.learn_chemical_properties+=);

        //std::cout << "reactions ..." << std::endl;
        // Calculate the kinetic rates of the reactions
        timeit(r = reactions.rates(properties), result.timing.learn_reaction_rates+=);

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

        //std::cout << "end of func u: " << tr(u) << std::endl;

        return 0;
    }

    auto jacobian(ChemicalState& state, double t, VectorConstRef u, MatrixRef res) -> int
    {
        // Calculate the sensitivity of the equilibrium state
        //sensitivity = equilibrium->sensitivity();

        if (options.learning.use_smart_equilibrium_solver)   sensitivity = smart_equilibrium.sensitivity();
        else    sensitivity = equilibrium.sensitivity();

        // Extract the columns of the kinetic rates derivatives w.r.t. the equilibrium and kinetic species
        drdne = cols(r.ddn, ies);
        drdnk = cols(r.ddn, iks);

        // Calculate the derivatives of `r` w.r.t. `be` using the equilibrium sensitivity
        drdbe = drdne * sensitivity.dndb;

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
            dqdbe = dqdne * sensitivity.dndb;

            // Assemble the partial derivatives of the source rates `q` w.r.t. to `u = [be nk]`
            dqdu << dqdbe, dqdnk;

            // Add the contribution of the source rates
            res += B * dqdu;
        }
        //std::cout << "jac res :\n" << res << std::endl;

        return 0;
    }
    /*
    auto showTree(const Index& step) -> void{
        //if (options.learning.use_smart_equilibrium_solver) smart_equilibrium->showTree(step);
    }
    */
};

SmartKineticSolver::SmartKineticSolver()
: pimpl(new Impl())
{}

SmartKineticSolver::SmartKineticSolver(const SmartKineticSolver& other)
: pimpl(new Impl(*other.pimpl))
{}

SmartKineticSolver::SmartKineticSolver(const ReactionSystem& reactions)
        : pimpl(new Impl(reactions))
{}

SmartKineticSolver::~SmartKineticSolver()
{}

auto SmartKineticSolver::operator=(SmartKineticSolver other) -> SmartKineticSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto SmartKineticSolver::setOptions(const SmartKineticOptions& options) -> void
{
    pimpl->setOptions(options);
}

auto SmartKineticSolver::setPartition(const Partition& partition) -> void
{
    pimpl->setPartition(partition);
}

auto SmartKineticSolver::addSource(const ChemicalState& state, double volumerate, std::string units) -> void
{
    pimpl->addSource(state, volumerate, units);
}

auto SmartKineticSolver::addPhaseSink(std::string phase, double volumerate, std::string units) -> void
{
    pimpl->addPhaseSink(phase, volumerate, units);
}

auto SmartKineticSolver::addFluidSink(double volumerate, std::string units) -> void
{
    pimpl->addFluidSink(volumerate, units);
}

auto SmartKineticSolver::addSolidSink(double volumerate, std::string units) -> void
{
    pimpl->addSolidSink(volumerate, units);
}

auto SmartKineticSolver::initialize(ChemicalState& state, double tstart) -> void
{
    pimpl->initialize(state, tstart);
}

} // namespace Reaktoro
