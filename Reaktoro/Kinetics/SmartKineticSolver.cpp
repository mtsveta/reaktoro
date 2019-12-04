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
#include <deque>

using namespace std::placeholders;

// Reaktoro includes
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/Exception.hpp>
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
struct SmartKineticNode{

    KineticState state;
    ChemicalState chemical_state;
    ChemicalProperties properties;
    EquilibriumSensitivity sensitivity;
    ChemicalVector rates;

    // The count of how often the element was successfully used for estimation
    unsigned int usage_count = 0;
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

    /// The chemical properties of the chemical system
    ChemicalProperties properties;

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
    std::deque<SmartKineticNode> tree;

    /// The sensitivity matrix of combined vector of elemental molar abundance and composition of kinetic species [be nk]
    Matrix benk_S;

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
        if(options.learning.use_smart_equilibrium_solver) std::cout << "       - incr time eq search  : " << stats.eq_search_time_total <<  " (" << stats.eq_search_time / stats.integraion_time_total * 100 << " %)" << std::endl;
        std::cout << " - time equilibration        : " << stats.equilibration_time_total <<  " (" << stats.equilibration_time_total / stats.time_total * 100 << " %)" << std::endl;
        std::cout << " - # of integration cell     : " << stats.conv_counter << " out of " << stats.conv_counter + stats.smart_counter
                  << " (" << float(stats.conv_counter) / float(stats.conv_counter + stats.smart_counter) * 100 << ") % " << std::endl;
        std::cout << "total time of SmartKineticSolver (without prop. evals) : " << stats.time_total - stats.prop_eval_time_total << std::endl;
        // <<  " -> speed-up of " << stats.time_total / (stats.time_total - stats.prop_eval_time_total) << std::endl;
        std::cout << "total time of SmartKineticSolver (without prop. evals and search) : "
                  << stats.time_total - stats.prop_eval_time_total - stats.search_time_total << std::endl;
        //<<  " -> speed-up of " << stats.time_total / (stats.time_total - stats.prop_eval_time_total - stats.search_time_total) << std::endl;
        if(options.learning.use_smart_equilibrium_solver) {
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
        if(options.learning.use_smart_equilibrium_solver)  std::cout << "       - incr time eq search : " << stats.eq_search_time <<  " (" << stats.eq_search_time / total_time_kinetics * 100 << " %)" << std::endl;
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

        // Allocate the memory for the sensitivity matrix
        benk_S.resize(Ee + Nk, Ee + Nk);

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


    auto initialize(ChemicalState& state, double tstart) -> void {

        // Extract the composition of the equilibrium and kinetic species
        const auto &n = state.speciesAmounts();
        nk = n(iks);
        ne = n(ies);

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
        ODEOptions options_ode = options.learning.ode;

        // Set the ODE problem and initialize the ODE solver
        ode.setProblem(problem);
        ode.setOptions(options_ode);
        ode.initialize(tstart, benk);

        // Set the options of the equilibrium solver
        smart_equilibrium.setOptions(options.learning.smart_equilibrium);
        equilibrium.setOptions(options.learning.equilibrium);

    }

    auto learn(ChemicalState& state, double t, double dt) -> void
    {
        SmartKineticResult res = {};

        // Initialize sensitivity matrix by the identity matrix
        benk_S.setIdentity();

        // Initialize the kinetics state with the data at times t0 and t0
        KineticState kin_state;
        kin_state.u0 = benk;
        kin_state.t0 = t;

        // Reconstract `benk` by the conventional numerical integration
        timeit(ode.solve(t, dt, benk, benk_S), result.timing.learn_integration +=);
        //timeit(ode.integrate(t, benk, t + dt, benk_S), result.timing.learn_integration +=); // This approach produces quite delayed estimations, figure out why?

        // Save the sensitivity values, the result time, and the obtain species' amount
        kin_state.t = t;
        kin_state.u = benk;
        kin_state.dudu0 = benk_S;

        // Update the composition of the amounts of equilibrium elements and kinetic species
        be = benk.head(Ee);
        nk = benk.tail(Nk);
        state.setSpeciesAmounts(nk, iks);

        // Evaluate chemical properties, reaction rate, and sensitivities (not needed, sinse it is evaluated in ode.solve())
        /*
        timeit(properties = state.properties() , result.timing.learn_chemical_properties +=);
        timeit(r = reactions.rates(properties), result.timing.learn_reaction_rates+=);
        timeit(sensitivity = equilibrium.sensitivity(), result.timing.learn_sensitivity +=);
        */

        // Save the kinetic state to the tree of learned states
        timeit(tree.emplace_back(SmartKineticNode{kin_state, state, properties, sensitivity, r}), result.timing.learn_storage+=);

        /*
        // TODO: remove, left from debugging
        std::cout << "t0     : " << kin_state.t0 << std::endl;
        std::cout << "t      : " << kin_state.t << std::endl;
        std::cout << "u0     : " << tr(kin_state.u0) << std::endl;
        std::cout << "u      : " << tr(kin_state.u) << std::endl;
        std::cout << "dndn0 : \n" << kin_state.dudu0 << std::endl;
        std::cout << "properties.lnActivities().val : \n" << tr(properties.lnActivities().val) << std::endl;
        std::cout << "r.val : \n" << r.val << std::endl;
        std::cout << "sensitivity.dndb : \n" << sensitivity.dndb << std::end;
        */
    }

    auto estimate(ChemicalState& state, double& t) -> void
    {
        // Skip estimation if no previous full computation has been done
        if(tree.empty())
            return;

        // Relative and absolute tolerance parameters
        auto kinetics_reltol = options.reltol;
        auto kinetics_abstol = options.abstol;
        auto kinetics_cutoff = options.cutoff;
        auto fraction_tol = kinetics_abstol * 1e-2; // TODO: is it important to multiply by 1e-2

        // Define initial state of the problem
        Vector benk0 = benk;

        // Comparison function based on the Euclidean distance
        auto distancefn = [&](const SmartKineticNode& a, const SmartKineticNode& b)
        {
            // Fetch benk0 from each TreeNode saved as u0 in ChemicalState
            const auto& benk0_a = a.state.u0;
            const auto& benk0_b = b.state.u0;

            return (benk0_a - benk0).squaredNorm() < (benk0_b - benk0).squaredNorm();  // TODO: We need to extend this later with T and P contributions too
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

        // Fetch the data stored in the reference element
        const auto& benk0_ref = it->state.u0;
        const auto& benk_ref = it->state.u;
        const auto& t0_ref = it->state.t0;
        const auto& t_ref = it->state.t;
        const auto& dndn0_ref = it->state.dudu0;

        // Algorithm:
        // the reference state : u, u0, S = du/du0, t, f = du/dt
        // new initial values  : u0_new
        // the predicted state : u_new = u + S * (u0_new - u0)

        // Clarification:
        // benk0 is new initial condition (u0_tilde)
        // it.state.u0    is the initial condition of reference vector (u0)
        // it.state.u    is already calculated by integration vector  (u)
        // it.state.dudu0 is the sensitivity w.r.t. the initial condition (S)

        // Perform smart estimation of benk
        Vector benk_new;
        benk_new.noalias() = benk_ref + dndn0_ref * (benk0 - benk0_ref);

        toc(1, result.timing.estimate_mat_vec_mul);

        //----------------------------------------------
        // Step 3: Checking the acceptance criterion
        //----------------------------------------------
        tic(2);

        // Fetch the be and nk unknowns from vector benk = [be; nk]
        VectorConstRef be_new = benk_new.head(Ee);
        VectorConstRef nk_new = benk_new.tail(Nk);

        // Fetch properties of the reference state
        const auto& state_ref = it->chemical_state;
        const auto& properties_ref = it->properties;
        const auto& sensitivity_ref = it->sensitivity;
        const auto& r_ref = it->rates;

        const auto& N = ies.size() + iks.size();
        const auto& be_ref = benk_ref.head(Ee);
        const auto& nk_ref = benk_ref.tail(Nk);
        const auto& n_ref = state_ref.speciesAmounts();
        const auto& ne_ref = n_ref(ies);

        // Get the sensitivity derivatives dln(a) / dn
        MatrixConstRef dlnadn_ref = properties_ref.lnActivities().ddn; // TODO: this line is assuming all species are equilibrium specie! get the rows and columns corresponding to equilibrium species
        VectorConstRef lna_ref = properties_ref.lnActivities().val;

        // Get the sensitivity derivatives w.r.t. the equilibrium species dln(a) / dne
        MatrixConstRef dlnaedne_ref = dlnadn_ref(ies, ies);
        VectorConstRef lnae_ref = lna_ref(ies);

        // Define vectors of equilibrium species ne, delta(ne), delta(lnae)
        Vector ne, delta_ne, delta_lnae;

        delta_ne.noalias() = sensitivity_ref.dndb * (be_new - be_ref);    // delta(ne) = dn/db * (be - be0)
        ne.noalias() = ne_ref + delta_ne;                             // ne = ne_ref + delta(ne)
        delta_lnae.noalias() = dlnaedne_ref * delta_ne;               // delta(delta_lnae) = d (ln ae) / d(ne)

        //std::cout << "delta_ne: " << tr(delta_ne) << std::endl;
        //std::cout << "delta_lna: " << tr(delta_lnae) << std::endl;
        //std::cout << "ne: " << tr(ne) << std::endl;

        // Fetch mole fractions
        const auto& x_ref = properties_ref.moleFractions().val;
        const auto& xe_ref = x_ref(ies);
        const auto& xk_ref = x_ref(iks);

        // -------------------------------------------------------------------------------------------------------------
        // Check the variations of equilibrium species
        // -------------------------------------------------------------------------------------------------------------
        bool equilibrium_variation_check = true;
        bool equilibrium_neg_amount_check = ne.minCoeff() > kinetics_cutoff;

        // Loop via equilibrium species and consider only big enough fractions
        for(int i = 0; i < xe_ref.size(); ++i){

            // If the fraction is too small, skip the variational check
            if(xe_ref(i, 0) < fraction_tol)
                continue;

            // Perform the variational check
            if(std::abs(delta_lnae[i]) > kinetics_abstol + kinetics_reltol * std::abs(lnae_ref[i])) {
                equilibrium_variation_check = false;
                result.estimate.failed_with_species = system.species(ies[i]).name();
                result.estimate.failed_with_amount = ne[i];
                result.estimate.failed_with_chemical_potential = lnae_ref[i] + delta_lnae[i]; // TODO: change to the chemical potential

                /*
                std::cout << "ne[i]     : " << ne[i] << std::endl;
                std::cout << "system.species(ies[i]).name()     : " << system.species(ies[i]).name() << std::endl;
                std::cout << "lnae_ref[i] + delta_lnae[i]     : " << lnae_ref[i] + delta_lnae[i] << std::endl;
                std::cout << "|delta_lna[i]| : " << std::abs(delta_lnae[i]) << std::endl;
                std::cout << "|lnae_ref[i]|  : " << std::abs(lnae_ref[i]) << std::endl;
                */
                break;
            }
        }

        /// -------------------------------------------------------------------------------------------------------------
        // Check the variations of kinetic species
        // -------------------------------------------------------------------------------------------------------------

        // Initialize delta_n = [delta_ne; delta_nk]
        Vector delta_nk;
        delta_nk.noalias() = nk_new - nk_ref;
        Vector delta_n;
        delta_n.resize(N);
        delta_n(ies) << delta_ne;
        delta_n(iks) << delta_nk;

        // Initialize reaction rates
        Vector delta_r = r_ref.ddn * delta_n;

        // TODO: loop for more then one kinetic species
        bool kinetics_r_variation_check = true;
        for(Index i = 0; i < xk_ref.size(); ++i){
            // If the fraction is too small, skip the variational check
            if(xk_ref(i, 0) < fraction_tol)
                continue;

            if(std::abs(delta_r.array()[i]) > kinetics_abstol + kinetics_reltol * std::abs(r_ref.val.array()[i])){
                /*
                std::cout << "var in i: " << iks[i] << std::endl;
                std::cout << "Species: " << system.species(iks[i]).name() << std::endl;
                std::cout << "|delta_r|  : " << delta_r.array().abs() << std::endl;
                std::cout << "|r_ref|    : " << r_ref.val.array().abs() << std::endl;
                std::cout << "|delta_r|  : " << delta_r.array().abs() << " < ";
                std::cout << "abstol + reltol * |r_ref|    : " << kinetics_abstol + kinetics_reltol * std::abs(r_ref.val.array()[i]) << std::endl;
                */
                kinetics_r_variation_check = false;
                break;
            }

        }

        if(!equilibrium_variation_check || !equilibrium_neg_amount_check || !kinetics_r_variation_check)
        //if(!equilibrium_variation_check || !kinetics_r_variation_check)
        {
            /*
            if(step>=1000) {
                //std::cout << "|xk_ref|  : " << TR(xk_ref)  << std::endl;
                //std::cout << fraction_tol  : " << fraction_tol  << std::endl;
                //std::cout << "|delta_r| : " << delta_r.array().abs() << std::endl;
                //std::cout << "|r_ref|   : " << r_ref.val.array().abs() << std::endl;
                //std::cout << "abstol + reltol * |r_ref|    : " << kinetics_abstol + kinetics_reltol * std::abs(r_ref.val.array()[i]) << std::endl;

                std::cout << "lnae_var_check     : " << equilibrium_variation_check;
                std::cout << ", eq_amount_check  : " << equilibrium_neg_amount_check;
                std::cout << ", r_var_check    : " << kinetics_r_variation_check << std::endl;

                getchar();
            }
            */
            /*
            std::cout << "lnae_var_check     : " << equilibrium_variation_check;
            std::cout << ", eq_amount_check  : " << equilibrium_neg_amount_check;
            std::cout << ", r_var_check    : " << kinetics_r_variation_check << std::endl;
            getchar();
            */
            return;
        }

        // Mark estimated result as accepted
        result.estimate.accepted = true;

        // Remember index of reference element used
        result.estimate.reference_state_index = std::distance(tree.begin(), it);

        // Assign small values to all the amount in the interval [cutoff, 0] (instead of mirroring above)
        //for(unsigned int i = 0; i < benk_new.size(); ++i) if(benk_new[i] >= kinetics_cutoff && benk_new[i] < 0) benk_new[i] = options.equilibrium.epsilon;
        for(unsigned int i = 0; i < benk_new.size(); ++i) if(benk_new[i] < 0) benk_new[i] = options.equilibrium.epsilon;

        // Update the solution of kinetic problem by new estimated value
        benk = benk_new;

        // Update properties by the reference one
        // properties = properties_ref;

        toc(2, result.timing.estimate_acceptance);

        // Increase the count of successfully used (for estimation) element
        it->usage_count++;

    }

    auto solve(ChemicalState& state, double t, double dt, VectorConstRef b) -> void
    {

        // Reset the result of the last smart equilibrium calculation
        result = {};

        tic(0);

        // ----------------------------------------------------------------------------------
        // Initialize kinetic problem
        // ----------------------------------------------------------------------------------

        // Extract the composition of the kinetic species from the state
        const auto &n = state.speciesAmounts();
        nk = n(iks);

        // Assemble the vector benk = [be nk]
        benk.resize(Ee + Nk);
        be = b - Ak * nk;
        benk.head(Ee) = be;
        benk.tail(Nk) = nk;

        // Initialise the chemical kinetics solver
        timeit(initialize(state, t, benk), result.timing.initialize=);

        // ----------------------------------------------------------------------------------
        // Smartly solve kinetic species
        // ----------------------------------------------------------------------------------

        // Perform a smart estimate for the chemical state
        timeit(estimate(state, t), result.timing.estimate =);

        // Perform a learning step if the smart prediction is not satisfactory
        if(!result.estimate.accepted)
            timeit(learn(state, t, dt), result.timing.learn =);

        // Extract the `be` and `nk` entries of the vector `benk`
        be = benk.head(Ee);
        nk = benk.tail(Nk);

        // Update the composition of the kinetic species
        state.setSpeciesAmounts(nk, iks);

        // ----------------------------------------------------------------------------------
        // Equilibrate equilibrium species
        // ----------------------------------------------------------------------------------
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

        /*
        std::cout << "smart equilibrium accepted? " << result.smart_equilibrium.estimate.accepted << ", index of ref elem in equilibrium " << result.smart_equilibrium.estimate.reference_state_index << std::endl; // << " out " << smart_equilibrium.tree().size()
        std::cout << "smart kinetics accepted   ? " << result.estimate.accepted << ", index of ref elem in kinetics " << result.estimate.reference_state_index << std::endl; // << " out of " << tree.size() <<
        getchar();
        */
        // If tree is of a certain size, do the cleanup
        if (tree.size() == 100){
            auto it = tree.begin();
            auto counter = 0;
            for (auto node : tree) {
                if (!node.usage_count) {
                    it = tree.erase(it); // return the iterator pointing on the elemnt after the removed one
                }
                else
                    it++; // increase the iterator
            }
        }
        // Print the tree
        if (tree.size() == 99 or tree.size() == 100 or tree.size() == 101 or t / dt = 49) {
            std::cout << "\nTree with reference elements for smart kinetics of the size " << tree.size() << std::endl;
            unsigned int i = 0;
            for (auto node : tree) {
                std::cout << "node " << i << " used " << node.usage_count
                          << std::endl;
                i += 1;
            }
            getchar();
        }

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
        if(false) // NOTE: for now, using equilibrium solver
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
            if(!res.optimum.succeeded) {
                state.setSpeciesAmounts(0.0);
                timeit(res = equilibrium.solve(state, T, P, be), result.timing.learn_equilibration +=);

            }
            if(!res.optimum.succeeded) {
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
        //timeit( sensitivity = options.use_smart_equilibrium_solver ? smart_equilibrium.sensitivity() : equilibrium.sensitivity(),
        //        result.timing.learn_sensitivity+=);
        sensitivity = equilibrium.sensitivity();

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
        //if(options.learning.use_smart_equilibrium_solver) smart_equilibrium->showTree(step);
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

auto SmartKineticSolver::initialize(ChemicalState& state, double tstart, VectorConstRef benk) -> void
{
    pimpl->initialize(state, tstart, benk);
}

auto SmartKineticSolver::result() const -> const SmartKineticResult&
{
    return pimpl->result;
}

auto SmartKineticSolver::properties() const -> const ChemicalProperties&
{
    return pimpl->properties;
}

auto SmartKineticSolver::solve(ChemicalState& state, double t, double dt, VectorConstRef b) -> void
{
    pimpl->solve(state, t, dt, b);
}

} // namespace Reaktoro
