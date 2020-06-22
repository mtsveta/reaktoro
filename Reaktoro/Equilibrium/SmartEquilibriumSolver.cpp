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

#include "SmartEquilibriumSolver.hpp"

// C++ includes
#include <algorithm>
#include <deque>
#include <numeric>
#include <tuple>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Profiling.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumResult.hpp>
#include <Reaktoro/ODML/ClusterConnectivity.hpp>
#include <Reaktoro/ODML/PriorityQueue.hpp>
#include <Reaktoro/Optimization/Canonicalizer.hpp>

namespace Reaktoro {
namespace detail {

/// Return the hash number of a vector.
/// https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector
template<typename Vec>
auto hash(Vec& vec) -> std::size_t
{
    using T = decltype(vec[0]);
    const std::hash<T> hasher;
    std::size_t seed = vec.size();
    for(const auto& i : vec)
        seed ^= hasher(i) + 0x9e3779b9 + (seed << 6u) + (seed >> 2u);
    return seed;
}

} // namespace detail


struct SmartEquilibriumSolver::Impl
{
    /// The chemical system instance
    ChemicalSystem system;

    /// The partition of the chemical system
    Partition partition;

    /// The canonicalizer used to determine primary and secondary species
    Canonicalizer canonicalizer;

    /// The options for the smart equilibrium calculation
    SmartEquilibriumOptions options;

    /// The result of the last smart equilibrium calculation.
    SmartEquilibriumResult result;

    /// The solver for the equilibrium calculations
    EquilibriumSolver solver;

    /// The chemical properties of the chemical system
    ChemicalProperties properties;

    /// The chemical potentials at the calculated equilibrium state
    ChemicalVector u;

    /// The amounts of the species in the chemical system
    Vector n;

    /// The amounts of the equilibrium species
    Vector ne;

    /// The amounts of the elements in the equilibrium partition
    Vector be;

    /// The storage for matrix du/db = du/dn * dn/db
    Matrix dudb;

    /// The storage for vector u(iprimary)
    Vector up;

    /// The storage for matrix Mbe = inv(u(iprimary)) * du(iprimary)/db.
    Matrix Mbe;

    /// The record of the knowledge database containing input, output, and derivatives data.
    struct Record
    {
        /// The temperature of the equilibrium state (in units of K).
        double T;

        /// The pressure of the equilibrium state (in units of Pa).
        double P;

        /// The amounts of elements in the equilibrium state (in units of mol).
        Vector be;

        /// The calculated equilibrium state at `T`, `P`, `be`.
        ChemicalState state;

        /// The chemical properties at the calculated equilibrium state.
        ChemicalProperties properties;

        /// The sensitivity derivatives at the calculated equilibrium state.
        EquilibriumSensitivity sensitivity;

        /// The matrix used to compute relative change of chemical potentials due to change in `be`.
        Matrix Mbe;
    };

    /// The cluster storing learned input-output data with same classification.
    struct Cluster
    {
        /// The indices of the primary species for this cluster.
        VectorXi iprimary;

        /// The hash of the indices of the primary species for this cluster.
        std::size_t label = 0;

        /// The records stored in this cluster with learning data.
        std::deque<Record> records;

        /// The priority queue for the records based on their usage count.
        PriorityQueue priority;
    };

    /// The database containing the learned input-output data.
    struct Database
    {
        /// The clusters containing the learned input-output data points.
        std::deque<Cluster> clusters;

        /// The connectivity matrix of the clusters to determine how we move from one to another when searching.
        ClusterConnectivity connectivity;

        /// The priority queue for the clusters based on their usage counts.
        PriorityQueue priority;
    };

    /// The database with learned input-output data points.
    Database database;

    /// Construct an SmartEquilibriumSolver::Impl instance with given chemical system.
    explicit Impl(const ChemicalSystem& system)
    : Impl(Partition(system))
    {
    }

    //----------------------------------------------------------------------------------
    // Old algorithms variables:
    //----------------------------------------------------------------------------------

    std::deque<Index> priority;

    std::deque<Index> ranking;

    /// The indices of the equilibrium species
    Indices ies;
    /// The indices of the equilibrium elements
    Indices iee;

    struct TreeNode
    {
        Vector be;
        ChemicalState state;
        ChemicalProperties properties;
        EquilibriumSensitivity sensitivity;
        Matrix Mb;
        VectorXi imajor;
    };
    /// The tree used to save the calculated equilibrium states and respective sensitivities
    std::deque<TreeNode> tree;

    /// Construct an SmartEquilibriumSolver::Impl instance with given partition of the chemical system.
    explicit Impl(const Partition& partition)
    : system(partition.system()), partition(partition),
      properties(partition.system()), solver(partition)
    {
        // Initialize the canonicalizer with the formula matrix Ae of the equilibrium species
        canonicalizer.compute(partition.formulaMatrixEquilibriumPartition());

        // Initialize indices of the equilibrium and kinetic species
        ies = partition.indicesEquilibriumSpecies();
        iee = partition.indicesEquilibriumElements();

    }

    /// Set the options for the equilibrium calculation.
    auto setOptions(const SmartEquilibriumOptions& options_) -> void
    {
        options = options_;

        // Tweak the options for the Gibbs energy minimization during learning operations.
        options.learning.hessian = GibbsHessian::Exact; // ensure the use of an exact Hessian of the Gibbs energy function
        options.learning.optimum.tolerance = 1e-10; // ensure the use of a stricter residual tolerance for the Gibbs energy minimization

        solver.setOptions(options.learning);
    }

    /// Learn how to perform a full equilibrium calculation (with tracking)
    auto learn_priority_based_acceptance_potential(ChemicalState& state, double T, double P, VectorConstRef be) -> void
    {
        //---------------------------------------------------------------------
        // GIBBS ENERGY MINIMIZATION CALCULATION DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(EQUILIBRIUM_STEP);

        // Calculate the equilibrium state using conventional Gibbs energy minimization approach
        // and store the result of the Gibbs energy minimization calculation performed during learning
        result.learning.gibbs_energy_minimization = solver.solve(state, T, P, be);

        result.timing.learn_gibbs_energy_minimization = toc(EQUILIBRIUM_STEP);

        //---------------------------------------------------------------------
        // CHEMICAL PROPERTIES UPDATE STEP DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(CHEMICAL_PROPERTIES_STEP);

        // Update the chemical properties of the system
        properties = solver.properties();

        result.timing.learn_chemical_properties = toc(CHEMICAL_PROPERTIES_STEP);

        //---------------------------------------------------------------------
        // ERROR CONTROL MATRICES ASSEMBLING STEP DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(ERROR_CONTROL_MATRICES);

        // The amounts of the species at the calculated equilibrium state
        n = state.speciesAmounts();

        // The amounts of the equilibrium species amounts at the calculated equilibrium state
        ne = n(ies);

        const auto& x = properties.moleFractions();
        Vector xe = x.val(ies);

        // Define the canonicalizer based on formular matrix A and priority weights set by the species amounts
        const auto& A = system.formulaMatrix();
        const auto& Ae = partition.formulaMatrixEquilibriumPartition();

        Canonicalizer canonicalizer(Ae);
        canonicalizer.updateWithPriorityWeights(ne);

        // Set tolerances for species' elements and fractions
        const auto eps_n = options.amount_fraction_cutoff * sum(ne);
        const auto eps_x = options.mole_fraction_cutoff;

        // Get indices of major species 'imajor'
        auto imajorminor = canonicalizer.Q();
        const auto nummajor = std::partition(imajorminor.begin(), imajorminor.end(),
                                             [&](Index i) { return ne[i] >= eps_n && xe[i] >= eps_x; })
                              - imajorminor.begin();
        const auto imajor = imajorminor.head(nummajor);

        // Fetch chemical potentials and their derivatives
        const auto u = properties.chemicalPotentials();
        const auto& dndb = solver.sensitivity().dndb;
        const auto& dudn = u.ddn;
        const Matrix dudb = dudn * dndb;

        // Calculate matrix Mb for the error control based on the potentials
        const Vector um = u.val(imajor);
        const auto dumdb = rows(dudb, imajor);
        const auto dumdbe = dudb(imajor, iee);
        const Matrix Mb = diag(inv(um)) * dumdb;
        const Matrix Mbe = diag(inv(um)) * dumdbe;

        result.timing.learn_error_control_matrices = toc(ERROR_CONTROL_MATRICES);

        //---------------------------------------------------------------------
        // STORAGE STEP DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(STORAGE_STEP);

        // Store the computed solution into the knowledge tree
        tree.push_back({be, state, properties, solver.sensitivity(), Mbe, imajor});

        // Update the priority queue
        // ----------------------------------------------
        // Add new element in the priority queue
        priority.push_back(priority.size());
        // Set its rank to zero
        ranking.push_back(0);

        result.timing.learn_storage = toc(STORAGE_STEP);
    }

    /// Learn how to perform a full equilibrium calculation (with tracking)
    auto learn_nnsearch(ChemicalState& state, double T, double P, VectorConstRef be) -> void
    {
        //---------------------------------------------------------------------
        // GIBBS ENERGY MINIMIZATION CALCULATION DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(EQUILIBRIUM_STEP);

        // Calculate the equilibrium state using conventional Gibbs energy minimization approach
        // and store the result of the Gibbs energy minimization calculation performed during learning
        result.learning.gibbs_energy_minimization = solver.solve(state, T, P, be);

        result.timing.learn_gibbs_energy_minimization = toc(EQUILIBRIUM_STEP);

        //---------------------------------------------------------------------
        // CHEMICAL PROPERTIES UPDATE STEP DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(CHEMICAL_PROPERTIES_STEP);

        // Update the chemical properties of the system
        properties = solver.properties();

        result.timing.learn_chemical_properties = toc(CHEMICAL_PROPERTIES_STEP);

        //---------------------------------------------------------------------
        // STORAGE STEP DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(STORAGE_STEP);

        // Store the computed solution into the knowledge tree
        tree.push_back({be, state, properties, solver.sensitivity(), Matrix::Zero(10, 10), VectorXi::Zero(10)});

        result.timing.learn_storage = toc(STORAGE_STEP);

    }

    /// Learn how to perform a full equilibrium calculation (with tracking)
    auto learn(ChemicalState& state, double T, double P, VectorConstRef be) -> void
    {
        //---------------------------------------------------------------------
        // GIBBS ENERGY MINIMIZATION CALCULATION DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(EQUILIBRIUM_STEP);

        // Perform a full chemical equilibrium calculation
        result.learning.gibbs_energy_minimization = solver.solve(state, T, P, be);

        result.timing.learn_gibbs_energy_minimization = toc(EQUILIBRIUM_STEP);

        //---------------------------------------------------------------------
        // CHEMICAL PROPERTIES UPDATE STEP DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(CHEMICAL_PROPERTIES_STEP);

        properties = solver.properties();

        result.timing.learn_chemical_properties = toc(CHEMICAL_PROPERTIES_STEP);

        //---------------------------------------------------------------------
        // SENSITIVITY MATRIX COMPUTATION STEP DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(SENSITIVITY_STEP);

        const auto& sensitivity = solver.sensitivity();

        result.timing.learn_sensitivity_matrix = toc(SENSITIVITY_STEP);

        //---------------------------------------------------------------------
        // ERROR CONTROL MATRICES ASSEMBLING STEP DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(ERROR_CONTROL_MATRICES);

        // The indices of the equilibrium species and elements
        const auto& ies = partition.indicesEquilibriumSpecies();
        const auto& iee = partition.indicesEquilibriumElements();

        // The number of equilibrium species and elements
        const auto& Ne = partition.numEquilibriumSpecies();
        const auto& Ee = partition.numEquilibriumElements();

        // The amounts of the species at the calculated equilibrium state
        n = state.speciesAmounts();

        // The amounts of the equilibrium species amounts at the calculated equilibrium state
        ne = n(ies);

        // Update the canonical form of formula matrix Ae so that we can identify primary species
        canonicalizer.updateWithPriorityWeights(ne);

        // The order of the equilibrium species as (primary, secondary)
        const auto& iorder = canonicalizer.Q();

        // Assemble the vector of indices of equilibrium species as (primary, secondary)
        VectorXi ips(ies.size());
        for(auto i = 0; i < ips.size(); ++i)  // TODO: type Indices should be alias to VectorXi, to avoid such kind of codes
            ips[i] = ies[iorder[i]];

        // The number of primary species among the equilibrium species (Np <= Ne)
        const auto& Np = canonicalizer.numBasicVariables();

        // Store the indices of primary and secondary species in state
        state.equilibrium().setIndicesEquilibriumSpecies(ips, Np);

        // The indices of the primary species at the calculated equilibrium state
        VectorXiConstRef iprimary = ips.head(Np);

        // The chemical potentials at the calculated equilibrium state
        u = properties.chemicalPotentials();

        // Auxiliary references to the derivatives dn/db and du/dn
        const auto& dndb = solver.sensitivity().dndb;
        const auto& dudn = u.ddn;

        // Compute the matrix du/db = du/dn * dn/db
        dudb = dudn * dndb;

        // The vector u(iprimary) with chemical potentials of primary species
        up.noalias() = u.val(iprimary);

        // The matrix du(iprimary)/dbe with derivatives of chemical potentials (primary species only)
        const auto dupdbe = dudb(iprimary, iee);

        // Compute matrix Mbe = 1/up * dup/db
        Mbe.noalias() = diag(inv(up)) * dupdbe;

        result.timing.learn_error_control_matrices = toc(ERROR_CONTROL_MATRICES);

        //---------------------------------------------------------------------
        // STORAGE STEP DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(STORAGE_STEP);

        // Generate the hash number for the indices of primary species in the state
        const auto label = detail::hash(iprimary);

        // Find the index of the cluster that has same primary species
        auto iter = std::find_if(database.clusters.begin(), database.clusters.end(),
            [&](const Cluster& cluster) { return cluster.label == label; });

        // If cluster is found, store the new record in it, otherwise, create a new cluster
        if(iter < database.clusters.end())
        {
            auto& cluster = *iter;
            cluster.records.push_back({T, P, be, state, properties, sensitivity, Mbe});
            cluster.priority.extend();
        }
        else
        {
            // Create a new cluster
            Cluster cluster;
            cluster.iprimary = iprimary;
            cluster.label = label;
            cluster.records.push_back({T, P, be, state, properties, sensitivity, Mbe});
            cluster.priority.extend();

            // Append the new cluster in the database
            database.clusters.push_back(cluster);
            database.connectivity.extend();
            database.priority.extend();
        }

        result.timing.learn_storage = toc(STORAGE_STEP);
    }

    /// Estimate the equilibrium state using sensitivity derivatives (profiling the expences)
    auto estimate(ChemicalState& state, double T, double P, VectorConstRef be) -> void
    {
        // Set the estimate status to false at the beginning
        result.estimate.accepted = false;

        // Skip estimation if no cluster exists yet
        if(database.clusters.empty())
            return;

        // Auxiliary relative and absolute tolerance parameters
        const auto reltol = options.reltol;

        // The threshold used to determine elements with insignificant amounts
        const auto eps_b = options.amount_fraction_cutoff * sum(be);

        // The current set of primary species in the chemical state
        const auto& iprimary = state.equilibrium().indicesPrimarySpecies();

        Vector dbe;

        // The function that checks if a record in the database pass the error test.
        // It returns (`success`, `error`, `iprimaryspecies`), where
        // * `success` is true if error test succeeds, false otherwise.
        // * `error` is the first error value violating the tolerance
        // * `iprimaryspecies` is the index of the primary species that fails the error test
        auto pass_error_test = [&be, &dbe, &reltol, &eps_b](const Record& record) -> std::tuple<bool, double, Index>
        {
            using std::abs;
            using std::max;
            const auto& state0 = record.state;
            const auto& be0 = record.be;
            const auto& Mbe0 = record.Mbe;
            const auto& isue0 = state0.equilibrium().indicesStrictlyUnstableElements();

            dbe.noalias() = be - be0;

            // Check if state0 has strictly unstable elements (i.e. elements with zero amounts)
            // which cannot be used for Taylor estimation if positive amounts for those elements are given.
            if((dbe(isue0).array() > eps_b).any())
            {
                assert((be0(isue0).array() < eps_b).any()); // ensure this condition is never broken (effective during debug only)
                return { false, 9999, -1 };
            }

            double error = 0.0;
            const auto size = Mbe0.rows();
            for(auto i = 1; i <= size; ++i) {
                error = max(error, abs(Mbe0.row(size - i) * dbe)); // start checking primary species with least amount first
                if(error >= reltol)
                    return { false, error, size - i };
            }

            return { true, error, -1 };
        };

        // Generate the hash number for the indices of primary species in the state
        const auto label = detail::hash(iprimary);

        // The function that identifies the starting cluster index
        auto index_starting_cluster = [&]() -> Index
        {
            // If no primary species, then return number of clusters to trigger use of total usage counts of clusters
            if(iprimary.size() == 0)
                return database.clusters.size();

            // Find the index of the cluster with the same set of primary species (search those with highest count first)
            for(auto icluster : database.priority.order())
                if(database.clusters[icluster].label == label)
                    return icluster;

            // In no cluster with the same set of primary species if found, then return number of clusters
            return database.clusters.size();
        };

        // The index of the starting cluster
        const auto icluster = index_starting_cluster();

        // The ordering of the clusters to look for (starting with icluster)
        const auto& clusters_ordering = database.connectivity.order(icluster);

        //---------------------------------------------------------------------
        // SEARCH STEP DURING THE ESTIMATE PROCESS
        //---------------------------------------------------------------------
        tic(SEARCH_STEP);

        // Iterate over all clusters (starting with icluster)
        for(auto jcluster : clusters_ordering)
        {
            // Fetch records from the cluster and the order they have to be processed in
            const auto& records = database.clusters[jcluster].records;
            const auto& records_ordering = database.clusters[jcluster].priority.order();

            // Iterate over all records in current cluster (using the  order based on the priorities)
            for(auto irecord : records_ordering)
            {
                const auto& record = records[irecord];

                //---------------------------------------------------------------------
                // ERROR CONTROL STEP DURING THE ESTIMATE PROCESS
                //---------------------------------------------------------------------
                tic(ERROR_CONTROL_STEP);

                // Check if the current record passes the error test
                const auto [success, error, iprimaryspecies] = pass_error_test(record);

                result.timing.estimate_error_control += toc(ERROR_CONTROL_STEP);

                if(success)
                {
                    //---------------------------------------------------------------------
                    // TAYLOR PREDICTION STEP DURING THE ESTIMATE PROCESS
                    //---------------------------------------------------------------------
                    tic(TAYLOR_STEP);

                    const auto& ies = partition.indicesEquilibriumSpecies();
                    const auto& iee = partition.indicesEquilibriumElements();

                    // Fetch reference values
                    const auto& be0 = record.be;
                    const auto& n0 = record.state.speciesAmounts();
                    const auto& P0 = record.state.pressure();
                    const auto& T0 = record.state.temperature();
                    const auto& y0 = record.state.equilibrium().elementChemicalPotentials();
                    const auto& z0 = record.state.equilibrium().speciesStabilities();
                    const auto& ips0 = record.state.equilibrium().indicesEquilibriumSpecies();
                    const auto& Np0 = record.state.equilibrium().numPrimarySpecies();
                    const auto& dndb0 = record.sensitivity.dndb;
                    const auto& dndP0 = record.sensitivity.dndP;
                    const auto& dndT0 = record.sensitivity.dndT;

                    // Fetch reference values restricted to equilibrium species only
                    const auto& dnedbe0 = dndb0(ies, iee);
                    const auto& dnedbPe0 = dndP0(ies);
                    const auto& dnedbTe0 = dndT0(ies);
                    const auto& ne0 = n0(ies);

                    // Perform Taylor extrapolation
                    //ne.noalias() = ne0 + dnedbe0 * (be - be0);
                    ne.noalias() = ne0 + dnedbe0 * (be - be0) + dnedbPe0 * (P - P0) + dnedbTe0 * (T - T0);

                    // Check if all projected species amounts are positive
                    const double ne_min = min(ne);
                    const double ne_sum = sum(ne);
                    const auto eps_n = options.amount_fraction_cutoff * ne_sum;
                    if(ne_min <= -eps_n)
                        continue;

                    result.timing.estimate_search = toc(SEARCH_STEP);

                    // After the search is finished successfully
                    //---------------------------------------------------------------------

                    // Update the amounts of elements for the equilibrium species
                    n(ies) = ne;

                    // Update the chemical state result with estimated amounts
                    state = record.state; // ATTENTION: If this changes one day, make sure indices of equilibrium primary/secondary species, and indices of strictly unstable species/elements are also transfered from reference state to new state
                    state.setSpeciesAmounts(n);

                    // Update the chemical properties of the system
                    properties = record.properties;  // TODO: We need to estimate properties = properties0 + variation : THIS IS A TEMPORARY SOLUTION!!!

                    result.timing.estimate_taylor = toc(TAYLOR_STEP);

                    //---------------------------------------------------------------------
                    // DATABASE PRIORITY UPDATE STEP DURING THE ESTIMATE PROCESS
                    //---------------------------------------------------------------------
                    tic(PRIORITY_UPDATE_STEP);

                    // Increment priority of the current record (irecord) in the current cluster (jcluster)
                    database.clusters[jcluster].priority.increment(irecord);

                    // Increment priority of the current cluster (jcluster) with respect to starting cluster (icluster)
                    database.connectivity.increment(icluster, jcluster);

                    // Increment priority of the current cluster (jcluster)
                    database.priority.increment(jcluster);

                    result.timing.estimate_database_priority_update = toc(PRIORITY_UPDATE_STEP);

                    result.estimate.accepted = true;

                    return;
                }
            }
        }

        result.estimate.accepted = false;

    }

    /// Estimate the equilibrium state using sensitivity derivatives (profiling the expences)
    auto estimate_priority_based_acceptance_potential(ChemicalState& state, double T, double P, VectorConstRef be) -> void
    {
        result.estimate.accepted = false;

        // Skip estimation if no previous full computation has been done
        if(tree.empty())
            return;

        // Relative and absolute tolerance parameters
        const auto reltol = options.reltol;

        //---------------------------------------------------------------------
        // SEARCH STEP DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(SEARCH_STEP);

        Vector rs;
        // Amounts of species
        Vector n, ns;
        // Variations of the potentials and species
        Vector du, dus, dns;
        // Variations of the elements
        Vector dbe;
        // Species' fractions
        Vector x;

        // Check if an entry in the database pass the error test.
        // It returns (`success`, `error`, `ispecies`), where
        //   - `success` is true if error test succeeds, false otherwise.
        //   - `error` is the first error violating the tolerance
        //   - `ispecies` is the index of the species that fails the error test
        auto pass_error_test = [&](const auto& node) -> std::tuple<bool, double, Index>
        {
            using std::abs;
            using std::max;
            const auto& be0 = node.be;
            const auto& Mb0 = node.Mb;
            const auto& imajor = node.imajor;

            dbe.noalias() = be - be0;

            double error = 0.0;
            for(auto i = 0; i < imajor.size(); ++i) {
                error = max(error, abs(Mb0.row(i) * dbe));
                if(error >= reltol)
                    return { false, error, imajor[i] };
            }

            return { true, error, n.size() };
        };

        auto inode_prev = priority.begin();
        for(auto inode=priority.begin(); inode!=priority.end(); ++inode)
        {
            const auto& node = tree[*inode];
            const auto& imajor = node.imajor;

            //---------------------------------------------------------------------
            // ERROR CONTROL STEP DURING THE LEARNING PROCESS
            //---------------------------------------------------------------------
            tic(ERROR_CONTROL_STEP);

            const auto [success, error, ispecies] = pass_error_test(node);

            if(success)
            {
                result.timing.estimate_error_control += toc(ERROR_CONTROL_STEP);

                //---------------------------------------------------------------------
                // TAYLOR EXTRAPOLATION STEP DURING THE LEARNING PROCESS
                //---------------------------------------------------------------------
                tic(TAYLOR_STEP);

                const auto& be0 = node.be;
                const auto& P0 = node.state.pressure();
                const auto& T0 = node.state.temperature();
                const auto& n0 = node.state.speciesAmounts();
                const auto& dndb0 = node.sensitivity.dndb;
                const auto& dndP0 = node.sensitivity.dndP;
                const auto& dndT0 = node.sensitivity.dndT;

                // Fetch reference values restricted to equilibrium species only
                const auto& ne0 = n0(ies);
                const auto& dnedbe0 = dndb0(ies, iee);
                const auto& dnedbPe0 = dndP0(ies);
                const auto& dnedbTe0 = dndT0(ies);

                // Perform Taylor extrapolation
                ne.noalias() = ne0 + dnedbe0 * (be - be0) + dnedbPe0 * (P - P0) + dnedbTe0 * (T - T0);

                result.timing.estimate_taylor = toc(TAYLOR_STEP);

                // Check if all projected species amounts are positive
                const double ne_min = min(ne);
                const double ne_sum = sum(ne);
                const auto eps_n = options.amount_fraction_cutoff * ne_sum;

                if(ne_min <= -eps_n)
                    continue;

                result.timing.estimate_search = toc(SEARCH_STEP);

                //-----------------------------------------------------------------------
                // DATABASE PRIORITY AND RANKING UPDATE STEP DURING THE ESTIMATE PROCESS
                //-----------------------------------------------------------------------
                tic(PRIORITY_UPDATE_STEP);

                // Increase ranking of of the used node
                ranking[*inode] += 1;

                // Make sure the indicies in the priority are ordered such that:
                // rank[priority[i-1]] > rank[priority[i]]
                auto comp = [&](Index l, Index r) { return ranking[l] > ranking[r]; };
                if( !((inode == priority.begin()) || (ranking[*inode_prev] >= ranking[*inode])) ) {
                    std::stable_sort(priority.begin(), inode + 1, comp);
                }

                result.timing.estimate_database_priority_update = toc(PRIORITY_UPDATE_STEP);

                // Update the amounts of elements for the equilibrium species
                //state = node.state; // this line was removed because it was destroying kinetics simulations
                state.setSpeciesAmounts(ne, ies);

                // Update the chemical properties of the system
                properties = node.properties;  // FIXME: We actually want to estimate properties = properties0 + variation : THIS IS A TEMPORARY SOLUTION!!!
                result.estimate.accepted = true;
                return;
            }
            else {
                inode_prev = inode;
                continue;
            }
        }

        result.estimate.accepted = false;

    }

    /// Estimate the equilibrium state using sensitivity derivatives,
    /// where search is performed by the NN algorithm and acceptance cirteria is using n and ln(a)
    auto estimate_nnsearch_acceptance_based_on_lna(ChemicalState& state, double T, double P, VectorConstRef be) -> void
    {
        // Skip estimation if no previous full computation has been done
        if(tree.empty())
            return;

        //---------------------------------------------------------------------
        // SEARCH STEP DURING THE ESTIMATE PROCESS
        //---------------------------------------------------------------------
        tic(SEARCH_STEP);

        // Comparison function based on the Euclidean distance
        auto distancefn = [&](const TreeNode& a, const TreeNode& b)
        {
            const auto& be_a = a.be;
            const auto& be_b = b.be;
            return (be_a - be).squaredNorm() < (be_b - be).squaredNorm();  // TODO: We need to extend this later with T and P contributions too (Allan, 26.06.2019)
        };

        // Find the entry with minimum "input" distance
        auto it = std::min_element(tree.begin(), tree.end(), distancefn);

        result.timing.estimate_search = toc(SEARCH_STEP);

        //---------------------------------------------------------------------
        // TAYLOR STEP DURING THE ESTIMATE PROCESS
        //---------------------------------------------------------------------
        tic(TAYLOR_STEP);

        // Get all the data stored in the reference element
        const auto& be0 = it->be;
        const auto& state0 = it->state;
        const auto& properties0 = it->properties;
        const auto& sensitivity0 = it->sensitivity;
        const auto& n0 = state0.speciesAmounts();
        const auto& ne0 = n0(ies);

        // Get the sensitivity derivatives dln(a) / dn (assuming all species are equilibrium species)
        MatrixConstRef dlna_dn = properties0.lnActivities().ddn;
        VectorConstRef lna0 = properties0.lnActivities().val;

        // Get the sensitivity derivatives w.r.t. the equilibrium species dln(ae) / dne
        MatrixConstRef dlnae_dne = dlna_dn(ies, ies);
        VectorConstRef lnae0 = lna0(ies);

        /// Auxiliary vectors delta(n) and delta(lna) in estimate function version v0
        Vector dne, dlnae;
        dne.noalias() = sensitivity0.dndb(ies, iee) * (be - be0); // delta(ne) = dne/db * (b - b0)
        ne.noalias() = ne0 + dne;                                                   // n = n0 + delta(n)
        dlnae.noalias() = dlnae_dne * dne;                                          // delta(ln(a)) = d(lna)/dn * delta(n)

        result.timing.estimate_taylor = toc(TAYLOR_STEP);

        //---------------------------------------------------------------------
        // ERROR CONTROL STEP DURING THE ESTIMATE PROCESS
        //---------------------------------------------------------------------
        tic(ERROR_CONTROL_STEP);

        // Fetch mole fractions
        Vector xe = properties0.moleFractions().val(ies);

        // Perform the check for the negative amounts
        const bool amount_check = ne.minCoeff() > options.cutoff;

        // Check in the loop mole fractions, negative amount check, and variations of the ln(a)
        bool variation_check = true;
        for(Index i = 0; i < ne.size(); ++i)
        {
            // If the fraction is too small, skip the variational check
            if(xe[i] < options.mole_fraction_cutoff)
                continue;

            // Absolute tolerance parameters
            const auto abstol = 1e-2;
            // Perform the variational check
            if(std::abs(dlnae[i]) > abstol + options.reltol * std::abs(lnae0[i])) {
                variation_check = false;    // variation test failed
                result.estimate.failed_with_species = system.species(ies[i]).name();
                result.estimate.failed_with_amount = ne[i];
                result.estimate.failed_with_chemical_potential = lna0[i] + dlnae[i]; // TODO: change to the chemical potential

                break;
            }
        }

        result.timing.estimate_error_control = toc(ERROR_CONTROL_STEP);

        // Check if smart estimation failed with respect to variation of chemical potentials or amounts
        if(!variation_check || !amount_check){
            return;
        }

        // Set the output chemical state to the approximate amounts of species
        // Assign small values to all the amount in the interval [cutoff, 0] (instead of mirroring above)
        for(unsigned int i = 0; i < ne.size(); ++i) if(ne[i] < 0) ne[i] = options.learning.epsilon;

        // Update equilibrium species
        state.setSpeciesAmounts(ne, ies);

        // Set the estimate accepted status to true
        result.estimate.accepted = true;

    }

    /// Estimate the equilibrium state using sensitivity derivatives (profiling the expences) v.1
    auto estimate_nnsearch_acceptance_based_on_residual(ChemicalState& state, double T, double P, VectorConstRef be) -> void
    {
        result.estimate.accepted = false;

        // Skip estimation if no previous full computation has been done
        if(tree.empty())
            return;

        //---------------------------------------------------------------------
        // SEARCH STEP DURING THE ESTIMATE PROCESS
        //---------------------------------------------------------------------
        tic(SEARCH_STEP);

        // Comparison function based on the Euclidean distance
        auto distancefn = [&](const TreeNode& a, const TreeNode& b)
        {
            const auto& be_a = a.be;
            const auto& be_b = b.be;
            return (be_a - be).squaredNorm() < (be_b - be).squaredNorm();  // TODO: We need to extend this later with T and P contributions too (Allan, 26.06.2019)
        };

        // Find the entry with minimum "input" distance
        auto it = std::min_element(tree.begin(), tree.end(), distancefn);

        result.timing.estimate_search = toc(SEARCH_STEP);

        //---------------------------------------------------------------------
        // TAYLOR STEP DURING THE ESTIMATE PROCESS
        //---------------------------------------------------------------------
        tic(TAYLOR_STEP);

        // Get all the data stored in the reference element
        const auto& be0 = it->be;
        const auto& state0 = it->state;
        const auto& properties0 = it->properties;
        const auto& sensitivity0 = it->sensitivity;
        const auto& T0 = state0.temperature();
        const auto& n0 = state0.speciesAmounts();
        const auto& y0 = state0.equilibrium().elementChemicalPotentials();
        const auto& z0 = state0.equilibrium().speciesStabilities();
        const auto u0 = properties0.chemicalPotentials();
        const auto x0 = properties0.moleFractions();

        // Amounts of species related to equilibrium and kinetics species
        const auto& ne0 = n0(ies);

        // Constant to scale the residual of the equilibrium species
        const auto RT = universalGasConstant*T0;

        /// Auxiliary vectors restricted to the equilibrium species
        Vector ye, ze, xe, ue, re;
        Vector dne;

        dne.noalias() = sensitivity0.dndb(ies, iee) * (be - be0);

        ye.noalias() = y0(iee);
        ne.noalias() = ne0 + dne;
        ze.noalias() = z0(ies);

        result.timing.estimate_taylor = toc(TAYLOR_STEP);

        //---------------------------------------------------------------------
        // ERROR CONTROL STEP DURING THE ESTIMATE PROCESS
        //---------------------------------------------------------------------
        tic(ERROR_CONTROL_STEP);

        // Negative cutoff control
        bool neg_amount_check = ne.minCoeff() > options.cutoff;

        // If cutoff test didn't pass, estimation has failded
        if(!neg_amount_check){
            result.timing.estimate_error_control = toc(ERROR_CONTROL_STEP);
            return;
        }

        // Calculate u(bar) = u(ref) + dudT(ref)*dT + dudP(ref)*dP + dudn(ref)*dn
        ue = u0.val(ies) + u0.ddn(ies, ies) * dne;

        // Calculate x(bar) = x(ref) + dxdn(ref)*dn + dxdP(ref)*dP + dxdn(ref)*dn
        xe = x0.val(ies) + x0.ddn(ies, ies) * dne;

        MatrixConstRef Ae = partition.formulaMatrixEquilibriumPartition();

        // Calculate the equilibrium residuals of the equilibrium species
        re = abs(ue - tr(Ae)*ye - ze)/RT;  // TODO: We should actually collect the entries in u and z corresponding to equilibrium species

        // Eliminate species with mole fractions below cutoff from the residual analysis
        for(auto i = 0; i < ne.size(); ++i)
            if(xe[i] < options.mole_fraction_cutoff)
                re[i] = 0.0; // set their residuals to zero

        // Estimate the residual error of the trial Taylor approximation for the equilibrium species
        Index ispecies;
        const double error = re.maxCoeff(&ispecies);
        bool equilibrium_variation_check = error <= options.reltol;

        result.timing.estimate_error_control = toc(ERROR_CONTROL_STEP);

         // Check if smart estimation failed with respect to variation of chemical potentials or amounts
        if(!equilibrium_variation_check)
            return;

        // Update three properties of the chemical state (for better initial guess of the GEM problem)
        state.setSpeciesAmounts(ne, ies);


        // Set the estimate accepted status to true
        result.estimate.accepted = true;

        // Update the chemical properties of the system
        properties = properties0;  // FIXME: We actually want to estimate properties = properties0 + variation : THIS IS A TEMPORARY SOLUTION!!!
    }

    /// Solve the equilibrium problem with given initial state
    auto solve(ChemicalState& state) -> SmartEquilibriumResult
    {
        const auto& iee = partition.indicesEquilibriumElements();
        const auto& ies = partition.indicesEquilibriumSpecies();
        const auto T = state.temperature();
        const auto P = state.pressure();
        be = state.elementAmountsInSpecies(ies)(iee);
        return solve(state, T, P, be);
    }

    /// Solve the equilibrium problem with given problem definition
    auto solve(ChemicalState& state, const EquilibriumProblem& problem) -> SmartEquilibriumResult
    {
        const auto T = problem.temperature();
        const auto P = problem.pressure();
        const auto& iee = partition.indicesEquilibriumElements();
        be = problem.elementAmounts()(iee);
        return solve(state, T, P, be);
    }

    /// Solve the equilibrium problem with given state, and input temperature, pressure, and element amounts
    auto solve(ChemicalState& state, double T, double P, VectorConstRef be) -> SmartEquilibriumResult
    {
        tic(SOLVE_STEP);

        // Absolutely ensure an exact Hessian of the Gibbs energy function is used in the calculations
        setOptions(options);

        // Reset the result of the last smart equilibrium calculation
        result = {};

        //-----------------------------------------------------------------------------------------------------
        // ESTIMATE EQUILIBRIUM SPECIES IN CHEMICAL STATE DURING THE SMART-EQUILIBRIUM SOLVE STEP
        //-----------------------------------------------------------------------------------------------------
        tic(ESTIMATE_STEP);

        // Perform a smart estimate of the chemical state
        //estimate(state, T, P, be);
        //estimate_nnsearch_acceptance_based_on_lna(state, T, P, be);
        //estimate_nnsearch_acceptance_based_on_residual(state, T, P, be);
        estimate_priority_based_acceptance_potential(state, T, P, be);
        //result.estimate.accepted = false;

        result.timing.estimate=toc(ESTIMATE_STEP);

        //-----------------------------------------------------------------------------------------------------
        // LEARN/FULLY CALCULATE EQUILIBRIUM SPECIES IN CHEMICAL STATE DURING THE SMART-EQUILIBRIUM SOLVE STEP
        //-----------------------------------------------------------------------------------------------------
        tic(LEARN_STEP);

        // Perform a learning step if the smart prediction is not sactisfatory
        if(!result.estimate.accepted){
            //learn(state, T, P, be);
            //learn_nnsearch(state, T, P, be);
            learn_priority_based_acceptance_potential(state, T, P, be);
        }

        result.timing.learn = toc(LEARN_STEP);

        result.timing.solve = toc(SOLVE_STEP);

        return result;
    }

    auto solve(ChemicalState& state, double T, double P, VectorConstRef be, Index istep, Index icell) -> SmartEquilibriumResult
    {
        tic(SOLVE_STEP);

        // Absolutely ensure an exact Hessian of the Gibbs energy function is used in the calculations
        setOptions(options);

        // Reset the result of the last smart equilibrium calculation
        result = {};

        //-----------------------------------------------------------------------------------------------------
        // ESTIMATE EQUILIBRIUM SPECIES IN CHEMICAL STATE DURING THE SMART-EQUILIBRIUM SOLVE STEP
        //-----------------------------------------------------------------------------------------------------
        tic(ESTIMATE_STEP);

        // Perform a smart estimate of the chemical state
        estimate(state, T, P, be);
        //estimate_nnsearch_acceptance_based_on_lna(state, T, P, be);
        //estimate_nnsearch_acceptance_based_on_potential(state, T, P, be);
        //estimate_priority_based_acceptance_potential(state, T, P, be);

        result.timing.estimate=toc(ESTIMATE_STEP);

        //-----------------------------------------------------------------------------------------------------
        // LEARN/FULLY CALCULATE EQUILIBRIUM SPECIES IN CHEMICAL STATE DURING THE SMART-EQUILIBRIUM SOLVE STEP
        //-----------------------------------------------------------------------------------------------------
        tic(LEARN_STEP);

        // Perform a learning step if the smart prediction is not sactisfatory
        if(!result.estimate.accepted){
            //learn_nnsearch(state, T, P, be);
            //learn_priority_based_acceptance_potential(state, T, P, be);
            learn(state, T, P, be);
        }
        //std::cout << "accepted? " << result.estimate.accepted << std::endl;
        //getchar();

        result.timing.learn = toc(LEARN_STEP);

        result.timing.solve = toc(SOLVE_STEP);

        /*
        // Print
        if(istep == 0 || istep == 19 || istep == 2399)
        {
            const auto Ae = partition.formulaMatrixEquilibriumPartition();
            const auto& ies = partition.indicesEquilibriumSpecies();
            const auto& n = state.speciesAmounts();
            const auto& ne = n(ies);

            Vector res = abs(Ae*ne - be);
            if(icell == 0) std::cout << "res on istep " << istep << "\n" << tr(res) << std::endl;
            else std::cout << tr(res)  << std::endl;
        }
        */
        return result;
    }

    auto outputClusterInfo() const -> void
    {
        //std::cout << "***********************************************************************************" << std::endl;
        //std::cout << "Clusters ordered by order of their creation" << std::endl;
        //std::cout << "***********************************************************************************" << std::endl;

        Index i = 0;
        for(auto cluster : database.clusters)
        {
            std::cout << "CLUSTER #" << i << std::endl;
            std::cout << "  RANK OF CLUSTER: " << database.priority.priorities()[i] << std::endl;
            std::cout << "  PRIMARY SPECIES: ";
            for(auto j : database.clusters[i].iprimary)
                std::cout << system.species(j).name() << " ";
            std::cout << std::endl;
            std::cout << "  NUMBER OF RECORDS: " << database.clusters[i].records.size() << std::endl;
            //std::cout << "  RANK OF RECORDS: ";
            //for(auto j : database.clusters[i].priority.order())
            //    std::cout << database.clusters[i].priority.priorities()[j] << " ";
            std::cout << std::endl;
            std::cout << std::endl;
            i++;
        }
        /*
        for(auto i : database.priority.order())
        {
            std::cout << "CLUSTER #" << i << std::endl;
            std::cout << "  RANK OF CLUSTER: " << database.priority.priorities()[i] << std::endl;
            std::cout << "  PRIMARY SPECIES: ";
            for(auto j : database.clusters[i].iprimary)
                std::cout << system.species(j).name() << " ";
            std::cout << std::endl;
            std::cout << "  NUMBER OF RECORDS: " << database.clusters[i].records.size() << std::endl;
            std::cout << "  RANK OF RECORDS: ";
            for(auto j : database.clusters[i].priority.order())
                std::cout << database.clusters[i].priority.priorities()[j] << " ";
            std::cout << std::endl;
            std::cout << "  NEXT CLUSTER: " << database.connectivity.order(i)[1] << std::endl;
            std::cout << "  NEXT CLUSTER PRIMARY SPECIES: ";
            for(auto j : database.clusters[database.connectivity.order(i)[1]].iprimary)
                std::cout << system.species(j).name() << " ";
            std::cout << std::endl;
            std::cout << std::endl;
        }
        */
    }
};

SmartEquilibriumSolver::SmartEquilibriumSolver()
{
    RuntimeError("Cannot proceed with EquilibriumSolver().",
                 "EquilibriumSolver() constructor is deprecated. "
                 "Use constructors EquilibriumSolver(const ChemicalSystem&) or EquilibriumSolver(const Partition&) instead.");
}
SmartEquilibriumSolver::SmartEquilibriumSolver(const ChemicalSystem& system)
: pimpl(new Impl(Partition(system)))
{}

SmartEquilibriumSolver::SmartEquilibriumSolver(const Partition& partition)
: pimpl(new Impl(partition))
{}

SmartEquilibriumSolver::SmartEquilibriumSolver(const SmartEquilibriumSolver& other)
: pimpl(new Impl(*other.pimpl))
{}

auto SmartEquilibriumSolver::operator=(SmartEquilibriumSolver other) -> SmartEquilibriumSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

SmartEquilibriumSolver::~SmartEquilibriumSolver() = default;

auto SmartEquilibriumSolver::setPartition(const Partition& partition) -> void
{
    RuntimeError("Cannot proceed with EquilibriumSolver::setPartition.",
                 "EquilibriumSolver::setPartition is deprecated. "
                 "Use constructor EquilibriumSolver(const Partition&) instead.");
}

auto SmartEquilibriumSolver::setOptions(const SmartEquilibriumOptions& options) -> void
{
    pimpl->setOptions(options);
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, double T, double P, VectorConstRef be) -> SmartEquilibriumResult
{
    return pimpl->solve(state, T, P, be);
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, const EquilibriumProblem& problem) -> SmartEquilibriumResult
{
    return pimpl->solve(state, problem);
}

auto SmartEquilibriumSolver::properties() const -> const ChemicalProperties&
{
    return pimpl->properties;
}

auto SmartEquilibriumSolver::result() const -> const SmartEquilibriumResult&
{
    return pimpl->result;
}

auto SmartEquilibriumSolver::outputClusterInfo() const -> void
{
    return pimpl->outputClusterInfo();
}

} // namespace Reaktoro

