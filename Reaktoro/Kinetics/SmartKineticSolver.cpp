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
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Common/TimeUtils.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Core/ChemicalOutput.hpp>
#include <Reaktoro/Core/ChemicalPlot.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumSolver.hpp>
#include <Reaktoro/Kinetics/KineticOptions.hpp>
#include <Reaktoro/Kinetics/KineticProblem.hpp>
#include <Reaktoro/Kinetics/KineticState.hpp>
#include <Reaktoro/Kinetics/KineticResult.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>
#include <Reaktoro/Transport.hpp>

namespace Reaktoro {

struct SmartKineticSolver::Impl
{
    /// The kinetically-controlled chemical reactions
    ReactionSystem reactions;

    /// The chemical system instance
    ChemicalSystem system;

    /// The partition of the species in the chemical system
    Partition partition;

    /// The options of the kinetic solver
    KineticOptions options;

    /// The equilibrium solver instance
    // EquilibriumSolver equilibrium;

    // The solver for solving the equilibrium equations using classical approach
    std::unique_ptr<EquilibriumSolver> equilibrium;

    /// The solver for solving the equilibrium equations using smart on-demand learning algorithm
    std::unique_ptr<SmartEquilibriumSolver> smart_equilibrium;

    /// The sensitivity of the equilibrium state
    EquilibriumSensitivity sensitivity;

    /// The sensitivity of the equilibrium state
    KineticState forward_sensitivity;

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
    Matrix Ae;

    /// The formula matrix of the kinetic species
    Matrix Ak;

    /// The stoichiometric matrix w.r.t. the equilibrium species
    Matrix Se;

    /// The stoichiometric matrix w.r.t. the kinetic species
    Matrix Sk;

    /// The coefficient matrix `A` of the kinetic rates
    Matrix A;

    /// The coefficient matrix `B` of the source rates
    Matrix B;

    /// The temperature of the chemical system (in units of K)
    double T;

    /// The pressure of the chemical system (in units of Pa)
    double P;

    /// The molar composition of the equilibrium species
    Vector ne;

    /// The molar composition of the kinetic species
    Vector nk;

    /// The molar abundance of the elements in the equilibrium species
    Vector be;

    /// The combined vector of elemental molar abundance and composition of kinetic species [be nk]
    Vector benk;

    /// The initial state of combined vector of elemental molar abundance and composition of kinetic species [be nk]
    Vector benk0;

    /// The vector of right-hand side  values at benk
    Vector benk_f;

    /// The matrix of sensitivities of combined vector [be nk]
    Matrix benk_S;

    /// The Jacobian matrix of combined vector [be nk]
    Matrix benk_J;

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

    /// The output instance of the kinetic path calculation
    ChemicalOutput output;

    /// The plots of the kinetic path calculation
    std::vector<ChemicalPlot> plots;

    double postequilibrate;

    /// A class used to store the node of tree for smart equilibrium calculations.
    struct TreeNode{

        KineticState state;
        ChemicalState chemical_state;
        ChemicalProperties properties;
        EquilibriumSensitivity sensitivity;
        ChemicalVector rates;

        Index num_predicted = 0;
        Index step = 0;
        Index cell = 0;

        TreeNode(KineticState state,
                 ChemicalState chemical_state,
                 ChemicalProperties properties,
                 EquilibriumSensitivity sensitivity,
                 ChemicalVector rates,
                 Index num_predicted, Index step, Index cell)
                : state(state), chemical_state(chemical_state),
                  properties(properties), sensitivity(sensitivity), rates(rates),
                  num_predicted(num_predicted), step(step), cell(cell){}

        TreeNode& operator=(const TreeNode& other){
            this->state = other.state;
            this->chemical_state = other.chemical_state;
            this->properties = other.properties;
            this->sensitivity = other.sensitivity;
            this->num_predicted = other.num_predicted;
            this->step = other.step;
            this->cell = other.cell;

            return *this;
        }


    };

    double dt_in = 0;

    struct Statistics{

        // Vector that stores number of smartly predicted on each step
        std::vector<int> smart_numbers;

        // Counter of smartly predicted cells on the current step
        int smart_counter = 0;

        // Counter of conventionaly integrated cells on the current step
        int conv_counter = 0;

        // Total time spent to search for the reference state and access its confidence
        double confidence_time = 0.0;

        // Total time spent to estimate new states
        double estimation_time = 0.0;

        // Total time spent to search for the nearest states
        double search_time = 0.0;

        // Total time spent to integrate new states
        double integraion_time = 0.0;

        // Time for properties evals
        double prop_eval_time = 0.0;

        double eq_calls_time  = 0.0;
        double eq_search_time = 0.0;


        // Total time spent to integrate new states
        double equilibration_time = 0.0;

        // Total time spent to integrate new states
        double kinetics_time = 0.0;

        // Total time spent to search for the reference state and access its confidence
        double confidence_time_total = 0.0;

        // Total time spent to estimate new states
        double estimation_time_total = 0.0;

        // Total time spent to search for the nearest states
        double search_time_total = 0.0;

        // Total time spent to integrate new states
        double integraion_time_total = 0.0;

        // Total time spent to integrate new states
        double equilibration_time_total = 0.0;

        // Total time for properties evals
        double prop_eval_time_total = 0.0;

        double eq_calls_time_total  = 0.0;
        double eq_search_time_total  = 0.0;

        // Total time spent to integrate new states
        double kinetics_time_total = 0.0;

        // Total time spent to integrate new states
        double time_total = 0.0;

    };

    Statistics stats;

    /// The tree used to save the calculated equilibrium states and respective sensitivities
    std::list<TreeNode> tree;

    Impl()
    {}

    Impl(const ReactionSystem& reactions, const KineticOptions& kin_options)
            : reactions(reactions), system(reactions.system())
    {
        // Initialize equilibrium solver based on the parameter
        // TODO: write it so that not both are initilized
        if (kin_options.equilibrium.smart.is_smart) {
            smart_equilibrium = std::make_unique<SmartEquilibriumSolver>(system);
            smart_equilibrium->setOptions(kin_options.equilibrium);
        }
        else {
            equilibrium = std::make_unique<EquilibriumSolver>(system);
            equilibrium->setOptions(kin_options.equilibrium);
        }
        setOptions(kin_options);

        //setPartition(Partition(system));

    }

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
        if (options.equilibrium.smart.is_smart) std::cout << "       - incr time eq search  : " << stats.eq_search_time_total <<  " (" << stats.eq_search_time / stats.integraion_time_total * 100 << " %)" << std::endl;
        std::cout << " - time equilibration        : " << stats.equilibration_time_total <<  " (" << stats.equilibration_time_total / stats.time_total * 100 << " %)" << std::endl;
        std::cout << " - # of integration cell     : " << stats.conv_counter << " out of " << stats.conv_counter + stats.smart_counter
                  << " (" << float(stats.conv_counter) / float(stats.conv_counter + stats.smart_counter) * 100 << ") % " << std::endl;
        std::cout << "total time of SmartKineticSolver (without prop. evals) : " << stats.time_total - stats.prop_eval_time_total << std::endl;
        // <<  " -> speed-up of " << stats.time_total / (stats.time_total - stats.prop_eval_time_total) << std::endl;
        std::cout << "total time of SmartKineticSolver (without prop. evals and search) : "
                  << stats.time_total - stats.prop_eval_time_total - stats.search_time_total << std::endl;
        //<<  " -> speed-up of " << stats.time_total / (stats.time_total - stats.prop_eval_time_total - stats.search_time_total) << std::endl;
        if (options.equilibrium.smart.is_smart) {
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
        if (options.equilibrium.smart.is_smart)  std::cout << "       - incr time eq search : " << stats.eq_search_time <<  " (" << stats.eq_search_time / total_time_kinetics * 100 << " %)" << std::endl;
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
    auto setOptions(const KineticOptions& options_) -> void
    {
        // Initialise the options of the kinetic solver
        options = options_;
        /*
        options.smart.is_smart = opts.smart.is_smart;
        options.smart.pertubation_tol = opts.smart.pertubation_tol;
        options.smart.abstol = opts.smart.abstol;
        options.smart.reltol = opts.smart.reltol;

        options.equilibrium.smart.is_smart = opts.equilibrium.smart.is_smart;
        options.equilibrium.smart.reltol = opts.equilibrium.smart.reltol;
        options.equilibrium.smart.abstol = opts.equilibrium.smart.abstol;
        options.equilibrium.track_statistics = opts.equilibrium.track_statistics;
        */

    }

    auto setPartition(const Partition& partition_) -> void
    {
        // Initialise the partition member
        partition = partition_;

        // Set the partition of the equilibrium solver
        if (options.equilibrium.smart.is_smart)
            smart_equilibrium->setPartition(partition);
        else
            equilibrium->setPartition(partition);

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

        // Initialise the formula matrix of the kinetic partitions
        Indices ike = partition.indicesKineticElements();
        Matrix Ak_ = partition.formulaMatrixKineticPartition();

        if (Nk) {
            //std::cout << "Ae  :" << Ae << std::endl;
            //std::cout << "A   :" << system.formulaMatrix() << std::endl;
            //std::cout << "ike :"; for (auto elem : ike) std::cout << elem << " "; std::cout << std::endl;
            //std::cout << "Ak_ :" << Ak_ << std::endl;
            //std::cout << "Ak :" << Ak << std::endl;

            Ak = Matrix::Zero(system.numElements(), Nk);
            for (unsigned i = 0; i < ike.size(); i++)   Ak.row(ike[i]) << Ak_.row(i);
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
            std::cout << "setPartition function of SmartKineticSolver:" << std::endl;
            std::cout << "ies :"; for (auto elem : ies) std::cout << elem << " "; std::cout << std::endl;
            std::cout << "iks :"; for (auto elem : iks) std::cout << elem << " "; std::cout << std::endl;
            std::cout << "S \n" << reactions.stoichiometricMatrix() << std::endl;
            std::cout << "Se \n" << Se << std::endl;
            std::cout << "Sk \n" << Sk << std::endl;

            std::cout << "A \n" << A << std::endl;
            std::cout << "Ae \n" << Ae << std::endl;
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

    auto initialize(ChemicalState& state, double tstart, VectorConstRef b) -> void
    {
        // Initialise the temperature and pressure variables
        T = state.temperature();
        P = state.pressure();

        // Extract the composition of the equilibrium and kinetic species
        const auto& n = state.speciesAmounts();
        ne = n(ies);
        nk = n(iks);

        // Let kinetics know about the change in transport
        /*
        std::cout << "b   : " << tr(b) << std::endl;
        std::cout << "Ae   : " << Ae << std::endl;
        std::cout << "ne   : " << ne << std::endl;
        std::cout << "Ak   : " << Ak << std::endl;
        std::cout << "nk   : " << nk << std::endl;
        std::cout << "Ak * nk   : " << tr(Ak * nk) << std::endl;
        std::cout << "Ae * ne   : " << tr(Ae * ne) << std::endl;
        std::cout << "iks :"; for (auto elem : iks) std::cout << elem << " "; std::cout << std::endl;
        std::cout << "ies :"; for (auto elem : ies) std::cout << elem << " "; std::cout << std::endl;
        */
        be = b - Ak * nk;
        //std::cout << "b   : " << tr(b) << std::endl;
        //std::cout << "be : " << tr(be) << std::endl;
        //std::cout << "bk  : " << tr(Ak * nk) << std::endl;

        // Assemble the vector benk = [be nk]
        benk.resize(Ee + Nk);
        benk.head(Ee) = be;
        benk.tail(Nk) = nk;

        // Initialize the initial state benk0
        benk0.resize(Ee + Nk);
        benk0 = benk;

        // Allocate the memory for the sensitivity matrix
        benk_S.resize(Ee + Nk, Ee + Nk);
        benk_S.setIdentity();

        // Allocate the memory for the sensitivity matrix
        benk_J.resize(Ee + Nk, Ee + Nk);

        // Allocate the memory for the right-hand side
        benk_f.resize(Ee + Nk);

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
        if (options.equilibrium.smart.is_smart) smart_equilibrium->setOptions(options.equilibrium);
        else                                    equilibrium->setOptions(options.equilibrium);

        // Initialize Jacobian and the right-hand side at the t0
        //problem.jacobian(tstart, benk0, benk_J);
        //problem.function(tstart, benk0, benk_f);

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

        // Initialize the initial state benk0
        benk0.resize(Ee + Nk);
        benk0 = benk;

        // Allocate the memory for the sensitivity matrix
        benk_S.resize(Ee + Nk, Ee + Nk);
        benk_S.setIdentity();

        // Allocate the memory for the sensitivity matrix
        benk_J.resize(Ee + Nk, Ee + Nk);

        // Allocate the memory for the right-hand side
        benk_f.resize(Ee + Nk);

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
        //equilibrium->setOptions(options.equilibrium);
        if (options.equilibrium.smart.is_smart) smart_equilibrium->setOptions(options.equilibrium);
        else                                    equilibrium->setOptions(options.equilibrium);

        // Initialize Jacobian and the right-hand side
        // problem.jacobian(tstart, benk0, benk_J);
        // problem.function(tstart, benk0, benk_f);

    }

    auto estimate(double& t, KineticResult& res) -> bool
    {
        // If the tree is empty abort estimation
        if(tree.empty())
            return false;

        // Comparison function based on the Euclidean distance
        auto compare_n0_tk = [&](const TreeNode& a, const TreeNode& b)
        {
            const auto& n0_a = a.state.n0;
            const auto& n0_b = b.state.n0;

            const auto& tk_a = a.state.tk;
            const auto& tk_b = b.state.tk;

            double dist_to_a = std::pow((n0_a - benk0).squaredNorm(), 2.0) + std::pow(t - tk_a, 2);
            double dist_to_b = std::pow((n0_b - benk0).squaredNorm(), 2.0) + std::pow(t - tk_b, 2);

            return dist_to_a < dist_to_b;
        };

        // Find the reference element (closest to the new state be)
        auto closest_n0 = std::min_element(tree.begin(), tree.end(), compare_n0_tk);

        // Get all the data stored in the reference element
        const auto& n0_ref = closest_n0->state.n0;
        const auto& tk_ref = closest_n0->state.tk;

        double distance = std::sqrt((benk0 - n0_ref).squaredNorm() + std::pow(t - tk_ref, 2));
        auto confident = distance < options.smart.pertubation_tol;
        //std::cout << "\nn0     : " << tr(n0) << " " << tk << std::endl;
        //std::cout << "n0_ref : " << tr(n0_ref) << " " << tk_ref << std::endl;
        //std::cout << "\ndis = " << distance << " > tol = " << options.smart.pertubation_tol << std::endl;

        if (confident) // Try estimating using sensitivities
        {
            Time start = time();
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
            const auto& nk1_ref = closest_n0->state.nk1;
            const auto& tk1_ref = closest_n0->state.tk1;
            const auto& dndn0 = closest_n0->state.dndn0;
            const auto& dndt = closest_n0->state.dndt;
            double dt = tk1_ref - tk_ref;

            benk = nk1_ref + dndn0 * (benk0 - n0_ref) + dndt * (t + tk_ref);

            //std::cout << "tk    : " << ref_node.state.tk << std::endl;
            //std::cout << "tk1    : " << ref_node.state.tk1 << std::endl;
            //std::cout << "tfinal : " << tfinal << std::endl;
            //std::cout << "dt     : " << dt << std::endl;
            //std::cout << "benk               : " << tr(benk) << std::endl;
            //std::cout << "benk (estimation)  : " << tr(benk) << std::endl;

            t += dt;

            // Profile estimation time of the cell
            stats.estimation_time += elapsed(start);
            res.smart.stats.time_estimate = stats.estimation_time;

            // Count the smartly predicted cell
            stats.smart_counter += 1;
            // Update the number of predicted cells
            closest_n0->num_predicted += 1;

            return true;
        }
        else{
            return false;
        }
    }

    auto learn(ChemicalState& state, double t, double tfinal, Index step, Index cell, KineticResult& res) -> bool
    {
        //Time start = time();

        // Initialize the kinetics state with the data at times t0 and tk
        KineticState kin_state;
        kin_state.n0 = benk0;
        kin_state.nk = benk;
        kin_state.tk = t;

        //ode.integrate(t, benk, tfinal, benk_S, benk_J, benk_f);
        ode.solve(t, tfinal - t, benk, benk_S, benk_J, benk_f);

        // Updating the chemical state
        //be = benk.head(Ee);
        nk = benk.tail(Nk);

        // Update the composition of the kinetic species
        state.setSpeciesAmounts(nk, iks);

        /*
        double t0 = t;
        double dt = 0.0, dt_in = 0.0;

        while (t <= tfinal){
            ode.integrate_(t, benk, tfinal, dt, benk_S, benk_J, benk_f);
            std::cout << "dt = " << dt << std:: endl;
            dt = (dt > dt_in) ? dt : dt_in;
        }
        */

        // Save the sensitivity values, the result time, and the obtain species' amount
        kin_state.tk1 = t;
        kin_state.nk1 = benk;
        kin_state.dndn0 = benk_S;
        kin_state.dndt = benk_f;
        //state.dfdn = benk_J;

        //std::cout << "\ntk     : " << state.tk << std::endl;
        //std::cout << "tk1    : " << state.tk1 << std::endl;
        //std::cout << "tfinal : " << tfinal << std::endl;

        // Save the kinetic state as the tree node
        ChemicalProperties prop;
        EquilibriumSensitivity eq_sens;

        Time start = time();
        // TODO: do we need to evaluate them here?
        if (options.equilibrium.smart.is_smart) {
            prop = smart_equilibrium->properties();
            eq_sens = smart_equilibrium->sensitivity();
        }
        else {
            prop = equilibrium->properties();
            eq_sens = equilibrium->sensitivity();
        }
        stats.prop_eval_time += elapsed(start);
        //ChemicalVector rates = reactions.rates(prop);
        //std::cout << "rates : " << rates.val << std::endl;
        //std::cout << "r     : " << r.val << std::endl;

        //node.properties = properties;
        //node.sensitivity = sensitivity;
        // TODO: can I initialize with state.properties() and state.sensitivities()?
        //tree.emplace_back(TreeNode(kin_state, state, properties, sensitivity, step, cell, 0));
        tree.emplace_back(TreeNode(kin_state, state, prop, eq_sens, r, step, cell, 0));
        //std::cout << "tk    : " << state.tk << std::endl;
        //std::cout << "tk1    : " << state.tk1 << std::endl;
        //std::cout << "tfinal : " << tfinal << std::endl;
        //std::cout << "dt     : " << state.tk1 - state.tk << std::endl;
        //std::cout << "benk (integration) : " << tr(benk) << std::endl;
        //std::cout << "time for integrating     : " << elapsed(start) << std::endl;
        //if (confident) std::cout << "error (estim - integr)   : " << (benk_est - benk).squaredNorm() << std::endl;

        // Profile estimation time of the cell
        //stats.integraion_time += elapsed(start);
        //res.smart.stats.time_integrate = stats.integraion_time;
        res.smart.addLearningIndex(cell);

        // Count the conventionally integrated cell
        stats.conv_counter += 1;
    }
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
        */
        return t;
    }


    auto solveSmart(ChemicalState& state, double& t, double tfinal, Index step, Index cell, KineticResult& res) -> void
    {
        Time start = time();
        bool result = estimate(t, res);
        stats.estimation_time += elapsed(start);
        res.smart.stats.time_estimate = stats.estimation_time;

        start = time();
        if (result)
            return;
        else
            learn(state, t, tfinal, step, cell, res);
        stats.integraion_time += elapsed(start);
        res.smart.stats.time_estimate = stats.integraion_time;

    }
    auto estimateUsingResult(ChemicalState& state, double& t, KineticResult& res, Index step, Index cell) -> bool
    {
        // If the tree is empty abort estimation
        if(tree.empty())
            return false;

        // Comparison function based on the Euclidean distance
        auto compare_nk = [&](const TreeNode& a, const TreeNode& b)
        {
            const auto& nk1_a = a.state.nk;
            const auto& nk1_b = b.state.nk;

            double dist_to_a = std::pow((nk1_a - benk).squaredNorm(), 0.5);
            double dist_to_b = std::pow((nk1_b - benk).squaredNorm(), 0.5);

            return dist_to_a < dist_to_b;
        };
        // Comparison function based on the Euclidean distance
        auto compare_nk_tk = [&](const TreeNode& a, const TreeNode& b)
        {
            const auto& nk1_a = a.state.nk;
            const auto& nk1_b = b.state.nk;

            const auto& tk_a = a.state.tk;
            const auto& tk_b = b.state.tk;

            double dist_to_a = std::sqrt(std::pow((nk1_a - benk).squaredNorm(), 2.0) + std::pow(t - tk_a, 2));
            double dist_to_b = std::sqrt(std::pow((nk1_a - benk).squaredNorm(), 2.0) + std::pow(t - tk_b, 2));

            return dist_to_a < dist_to_b;
        };
        // Comparison function based on the Euclidean distance
        auto compare_n0_tk = [&](const TreeNode& a, const TreeNode& b)
        {
            const auto& n0_a = a.state.n0;
            const auto& n0_b = b.state.n0;

            const auto& tk_a = a.state.tk;
            const auto& tk_b = b.state.tk;

            double dist_to_a = std::pow((n0_a - benk0).squaredNorm(), 2.0) + std::pow(t - tk_a, 2);
            double dist_to_b = std::pow((n0_b - benk0).squaredNorm(), 2.0) + std::pow(t - tk_b, 2);

            return dist_to_a < dist_to_b;
        };

        Time start = time();

        // Find the reference element (nearest to the new state benk)
        auto closest_n0 = std::min_element(tree.begin(), tree.end(), compare_nk);

        stats.search_time += elapsed(start);
        res.smart.stats.time_search = stats.search_time;

        start = time();

        // Get all the data stored in the reference element
        const auto& benk0_ref = closest_n0->state.n0;
        const auto& benkk1_ref = closest_n0->state.nk1;
        const auto& benkk_ref = closest_n0->state.nk;
        const auto& tk_ref = closest_n0->state.tk;
        const auto& tk1_ref = closest_n0->state.tk1;
        const auto& dndn0_ref = closest_n0->state.dndn0;
        const auto& dndt_ref = closest_n0->state.dndt;
        const ChemicalState& state_ref = closest_n0->chemical_state;
        const ChemicalProperties& prop_ref = closest_n0->properties;
        const EquilibriumSensitivity& eqsens_ref = closest_n0->sensitivity;
        const ChemicalVector& r_ref = closest_n0->rates;

        const auto& be_ref = benkk_ref.head(Ee);
        const auto& n_ref = state_ref.speciesAmounts();
        const auto& N = ies.size() + iks.size();
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

        double dt = tk1_ref - tk_ref;
        benk = benkk1_ref + dndn0_ref * (benk0 - benk0_ref) + dndt_ref * (t - tk_ref);
        //std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << tr(benk) << std::endl;

        VectorConstRef be_new = benk.head(Ee);
        VectorConstRef nk_new = benk.tail(Nk);

        double dist_nk = std::pow((nk_ref - nk_new).squaredNorm(), 0.5);
        //std::cout << "nk_new: " << tr(nk_new) << std::endl;
        //std::cout << "nk_ref: " << tr(nk_ref) << std::endl;
        //std::cout << "dist(nk_new, nk_ref): " << dist_nk  << std::endl;

        // Update the composition of the kinetic species
        ChemicalState& state_new = state;
        state_new.setSpeciesAmounts(nk_new, iks);
        ChemicalProperties prop_new;
        Time start_prop = time();
        if (options.equilibrium.smart.is_smart) prop_new = smart_equilibrium->properties();
        else                                    prop_new = equilibrium->properties();
        stats.prop_eval_time += elapsed(start_prop);
        VectorConstRef lna_new = prop_new.lnActivities().val;
        VectorConstRef lnae_new = lna_new(ies);
        ChemicalVector r_new = reactions.rates(prop_new);

        // Retreval of kinetic and equilibrium species out of
        // Get the sensitivity derivatives dln(a) / dn
        MatrixConstRef dlnadn_ref = prop_ref.lnActivities().ddn; // TODO: this line is assuming all species are equilibrium specie! get the rows and columns corresponding to equilibrium species
        VectorConstRef lna_ref = prop_ref.lnActivities().val;

        // Get the sensitivity derivatives w.r.t. the equilibrium species dln(a) / dne
        MatrixConstRef dlnadne_ref = dlnadn_ref(ies, ies);
        VectorConstRef lnae_ref = lna_ref(ies);

        //VectorConstRef diff_lnae = lnae_ref - lnae_new;
        //double dist_lnae = std::pow((lnae_ref - lnae_new).squaredNorm(), 2.0);
        //std::cout << "lnae_new: " << tr(lnae_new) << std::endl;
        //std::cout << "lnae_ref: " << tr(lnae_ref) << std::endl;
        //std::cout << "diff_lnae: " << tr(diff_lnae) << std::endl;
        //std::cout << "dist(lnae_new, lnae_ref): " << dist_lnae  << std::endl;

        Vector delta_ne;
        delta_ne.noalias() = eqsens_ref.dndb * (be_new - be_ref);    // delta(n) = dn/db * (b - b0)
        //std::cout << "delta_ne: " << tr(delta_ne) << std::endl;
        Vector ne_new;
        ne_new.noalias() = ne_ref + delta_ne;                             // n = n0 + delta(n)
        Vector delta_lna;
        delta_lna.noalias() = dlnadne_ref * delta_ne;
        //std::cout << "delta_lna: " << tr(delta_lna) << std::endl;
        //std::cout << "ne_new: " << tr(ne_new) << std::endl;

        // const ChemicalState& state0 = it->state;
        // const ChemicalProperties& properties0 = it->properties;
        // const EquilibriumSensitivity& sensitivity0 = it->sensitivity;

        auto reltol = options.smart.reltol;
        auto abstol = options.smart.abstol;
        //auto reltol = options.equilibrium.smart.reltol;
        //auto abstol = options.equilibrium.smart.abstol;
        //const bool acceptance_eq = (delta_lna.array().abs() <= abstol + reltol * lna_ref.array().abs()).all();

        const double cutoff = -1e-5;    // relaxation parameter
        const bool amount_check = ne_new.minCoeff() > cutoff;
        const auto& x = prop_ref.moleFractions();
        //std::cout << "x: " << x << std::endl;

        // Check the variations of equilibrium species
        bool acceptance_eq = true;
        const double fraction_tol = abstol * 1e-2;
        for(int i = 0; i < ne_new.size(); ++i){
            // If the fraction is too small, skip the variational check
            if(x[i] < fraction_tol) continue;

            ///std::cout << "abs(delta_lna[i]): " << std::abs(delta_lna[i]) << std::endl;
            //std::cout << "abstol + reltol * std::abs(lnae_ref[i])): " << abstol + reltol * std::abs(lnae_ref[i]) << std::endl;

            // Perform the variational check
            if(std::abs(delta_lna[i]) > abstol + reltol * std::abs(lnae_ref[i])) {
                acceptance_eq = false;    // variation test failed
                i = ne_new.size();        // finish running through the loop
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
        reltol = options.smart.reltol;
        abstol = options.smart.abstol;

        // Initialize delta_n = [delta_ne; delta_nk]
        Vector delta_n;
        delta_n.resize(N);
        //std::cout << "delta_ne: " << tr(delta_ne) << std::endl;
        //std::cout << "nk_new - nk_ref: " << nk_new - nk_ref << std::endl;
        delta_n(ies) << delta_ne;
        delta_n(iks) << nk_new - nk_ref;
        //std::cout << "delta_n: " << tr(delta_n) << std::endl;

        // Compute delta_rk = drdn * delta_n
        Vector delta_rk;
        delta_rk.noalias() = r_ref.ddn * delta_n;
        //std::cout << "delta_rk.array().abs(): " << delta_rk.array().abs() << std::endl;
        //std::cout << "abstol + reltol * r_ref.val.array().abs(): " << abstol + reltol * r_ref.val.array().abs() << std::endl;
        const bool accaptance_kin = (delta_rk.array().abs() <= abstol + reltol * r_ref.val.array().abs()).all();

        stats.confidence_time += elapsed(start);
        res.smart.stats.time_confidence = stats.confidence_time;

        //std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << tr(benk) << std::endl;
        //bool confident = true;
        if (acceptance_eq && amount_check && accaptance_kin) {

            //std::cout << "benk (estimation) : " << tr(benk) << std::endl;
            //benk = benk_new;
            // Count the smartly predicted cell
            stats.smart_counter += 1;
            t += dt;
            return true;
        }else{
            //std::cout << "(" << step << ", " << cell << ") : ";
            //std::cout << "accept_eq  = " << acceptance_eq << "\t";
            //std::cout << "accapt_kin = " << accaptance_kin << std::endl;
            return false;
        }
    }

    auto solveSmartPostprocess(ChemicalState& state, double& t, double tfinal, Index step, Index cell, KineticResult& res) -> void
    {
        Time start = time();
        bool success = estimateUsingResult(state, t, res, step, cell);
        stats.estimation_time += elapsed(start);
        res.smart.stats.time_estimate = stats.estimation_time;

        start = time();
        if (success)
            return;
        else
            learn(state, t, tfinal, step, cell, res);
        stats.integraion_time += elapsed(start);
        res.smart.stats.time_integrate = stats.integraion_time;
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

            /*
            result.times.emplace_back(t);
            result.calcite.emplace_back(state.speciesAmount("Calcite"));
            result.hco3.emplace_back(state.speciesAmount("HCO3-"));
            result.ca.emplace_back(state.speciesAmount("Ca++"));
            */

            // Integration time step
            t = stepSmart(state, t, t1);
        }

        // Update the output with the final state
        if(output) output.update(state, t1);

        // Update the plots with the final state
        for(auto& plot : plots) plot.update(state, t1);

        /*
        result.times.emplace_back(t);
        result.calcite.emplace_back(state.speciesAmount("Calcite"));
        result.hco3.emplace_back(state.speciesAmount("HCO3-"));
        result.ca.emplace_back(state.speciesAmount("Ca++"));
        */
    }

    auto solve(ChemicalState& state, double t, double dt) -> void
    {
        // Initialise the chemical kinetics solver
        initialize(state, t);

        std::cout << "\nbenk0: " << tr(benk) << std::endl;
        std::cout << "t0: " << t << std::endl;

        // Integrate the chemical kinetics ODE from `t` to `t + dt`
        ode.solve(t, dt, benk);

        std::cout << "benk: " << tr(benk) << std::endl;
        std::cout << "t: " << t << std::endl;

        // Extract the `be` and `nk` entries of the vector `benk`
        be = benk.head(Ee);
        nk = benk.tail(Nk);

        // Update the composition of the kinetic species
        state.setSpeciesAmounts(nk, iks);

        // Update the composition of the equilibrium species
        //equilibrium->solve(state, T, P, be);
        //Time start = time();
        EquilibriumResult res;
        if (options.equilibrium.smart.is_smart) res += smart_equilibrium->solve(state, T, P, be);
        else                                    equilibrium->solve(state, T, P, be);
        //stats.eq_calls_time += elapsed(start);

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
            if (options.equilibrium.smart.is_smart) smart_equilibrium->solve(state, T, P, be, res.equilibrium, step, cell);
            else res.equilibrium += equilibrium->solve(state, T, P, be);

            stats.equilibration_time += elapsed(start);
            res.smart.stats.time_equilibrate = stats.equilibration_time;

            // Update the learning indices for the smart equilibrium
            if (options.equilibrium.smart.is_smart && !res.equilibrium.smart.succeeded){
                res.equilibrium.smart.addLearningIndex(cell);
                res.equilibrium.smart.tree_size++;
            }
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
        if (options.equilibrium.smart.is_smart)
            smart_equilibrium->setOptions(options.equilibrium);
        else
            equilibrium->setOptions(options.equilibrium);

        Time start = time();

        // Update the composition of the equilibrium species
        //equilibrium->solve(state, T, P, be);
        /*
        if (options.equilibrium.smart.is_smart)
            smart_equilibrium->solve(state, T, P, be, step, cell);
        else
            equilibrium->solve(state, T, P, be);
        */
        if (options.equilibrium.smart.is_smart)
            res.equilibrium += smart_equilibrium->solve(state, T, P, be, step, cell);
        else
            res.equilibrium += equilibrium->solve(state, T, P, be);

        stats.equilibration_time += elapsed(start);
        res.smart.stats.time_equilibrate = stats.equilibration_time;

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

    auto function(ChemicalState& state, double t, VectorConstRef u, VectorRef res) -> int
    {
        // Extract the `be` and `nk` entries of the vector [be, nk]
        be = u.head(Ee);
        nk = u.tail(Nk);

        // Check for non-finite values in the vector `benk`
        for(int i = 0; i < u.rows(); ++i)
            if(!std::isfinite(u[i]))
                return 1; // ensure the ode solver will reduce the time step

        // Update the composition of the kinetic species in the member `state`
        state.setSpeciesAmounts(nk, iks);

        // Solve the equilibrium problem using the elemental molar abundance `be`
        //auto result = equilibrium->solve(state, T, P, be);
        Time start = time();
        EquilibriumResult result;
        if (options.equilibrium.smart.is_smart)
            result = smart_equilibrium->solve(state, T, P, be);
        else
            result = equilibrium->solve(state, T, P, be);

        // Check if the calculation failed, if so, use cold-start
        if(!result.optimum.succeeded)
        {
            state.setSpeciesAmounts(0.0);
            //result = equilibrium->solve(state, T, P, be);
            if (options.equilibrium.smart.is_smart) result = smart_equilibrium->solve(state, T, P, be);
            else                                    result = equilibrium->solve(state, T, P, be);

        }
        stats.eq_calls_time += elapsed(start);
        stats.eq_search_time += result.smart.estimate_stats.time_search;

        // Assert the equilibrium calculation did not fail
        Assert(result.optimum.succeeded,
               "Could not calculate the rates of the species.",
               "The equilibrium calculation failed.");

        start = time();
        // Update the chemical properties of the system
        properties = state.properties();
        stats.prop_eval_time += elapsed(start);

        // Calculate the kinetic rates of the reactions
        r = reactions.rates(properties);

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
        //std::cout << "fun res :" << tr(res) << std::endl;

        return 0;
    }

    auto jacobian(ChemicalState& state, double t, VectorConstRef u, MatrixRef res) -> int
    {
        // Calculate the sensitivity of the equilibrium state
        //sensitivity = equilibrium->sensitivity();

        if (options.equilibrium.smart.is_smart)
            sensitivity = smart_equilibrium->sensitivity();
        else
            sensitivity = equilibrium->sensitivity();

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
    auto showTree(const Index& step) -> void{
        if (options.equilibrium.smart.is_smart) smart_equilibrium->showTree(step);
    }
};

SmartKineticSolver::SmartKineticSolver()
        : pimpl(new Impl())
{}

SmartKineticSolver::SmartKineticSolver(const ReactionSystem& reactions, const KineticOptions& kin_options)
        : pimpl(new Impl(reactions, kin_options))
{}

SmartKineticSolver::~SmartKineticSolver()
{}

auto SmartKineticSolver::operator=(SmartKineticSolver other) -> SmartKineticSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto SmartKineticSolver::console() -> void
{
    pimpl->console();
}

auto SmartKineticSolver::summarize() -> void
{
    pimpl->summarize();
}
auto SmartKineticSolver::refresh() -> void
{
    pimpl->refresh();
}

auto SmartKineticSolver::setOptions(const KineticOptions& options) -> void
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
auto SmartKineticSolver::initialize(ChemicalState& state, double tstart, VectorConstRef be) -> void
{
    pimpl->initialize(state, tstart, be);
}

/*
auto SmartKineticSolver::step(ChemicalState& state, double t) -> double
{
    return pimpl->step(state, t);
}
*/

auto SmartKineticSolver::stepSmart(ChemicalState& state, double t, double tfinal) -> double
{
    return pimpl->stepSmart(state, t, tfinal);
}
/*
auto SmartKineticSolver::solveSmart(ChemicalState& state, double t, double dt) -> double
{
    return pimpl->solveSmart(state, t, dt);
}
*/
/*
auto SmartKineticSolver::stepKinetics(ChemicalState& state, double t, double dt, VectorConstRef be) -> double
{
    return pimpl->stepKinetics(state, t, dt,be);
}
*/

auto SmartKineticSolver::solve(ChemicalState& state, double t, double dt) -> void
{
    pimpl->solve(state, t, dt);
}

auto SmartKineticSolver::solveSmart(ChemicalState& state, double& t, double dt, Index step, Index cell, KineticResult& result) -> void
{
    pimpl->solveSmart(state, t, dt, step, cell, result);
}
auto SmartKineticSolver::solveSmartPostprocess(ChemicalState& state, double& t, double dt, Index step, Index cell, KineticResult& result) -> void
{
    pimpl->solveSmartPostprocess(state, t, dt, step, cell, result);
}

auto SmartKineticSolver::solveEquilibriumKinetics(ChemicalState& state, double t, double dt, VectorConstRef b, Index step, Index cell,  KineticResult& result) -> void
{
    pimpl->solveEquilibriumKinetics(state, t, dt, b, step, cell, result);
}
auto SmartKineticSolver::solveKineticsEquilibrium(ChemicalState& state, double t, double dt, VectorConstRef b,  Index step, Index icell, KineticResult& result) -> void
{
    pimpl->solveKineticsEquilibrium(state, t, dt, b, step, icell, result);
}
auto SmartKineticSolver::solvePath(ChemicalState& state, double t, double tfinal, std::string units, KineticResult& result) -> void
{
    pimpl->solvePath(state, t, tfinal, units, result);
}
/*
auto SmartKineticSolver::solveCustom(ChemicalState& state, double t, double tfinal, std::string units, VectorConstRef be, KineticPathResult& result) -> void
{
    pimpl->solveCustom(state, t, tfinal, units, be, result);
}
 */
/*
auto SmartKineticSolver::solveCustom(ChemicalState& state, double t, double tfinal, std::string units, KineticPathResult& result) -> void
{
    pimpl->solveKineticsEquilibrium(state, t, tfinal, units, result);
}
*/

auto SmartKineticSolver::output() -> ChemicalOutput
{
    pimpl->output = ChemicalOutput(pimpl->reactions);
    return pimpl->output;
}

auto SmartKineticSolver::plot() -> ChemicalPlot
{
    pimpl->plots.push_back(ChemicalPlot(pimpl->reactions));
    return pimpl->plots.back();
}

auto SmartKineticSolver::plots(unsigned num) -> std::vector<ChemicalPlot>
{
    for(unsigned i = 0; i < num; ++i) plot();
    return pimpl->plots;
}
auto  SmartKineticSolver::showTree(const Index& step) -> void {
    pimpl->showTree(step);
}


} // namespace Reaktoro
