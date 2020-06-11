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
#include <tuple>

// Eigen includes
#include <Reaktoro/deps/eigen3/Eigen/LU>

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
#include <Reaktoro/Optimization/Canonicalizer.hpp>

namespace Reaktoro {

struct SmartEquilibriumSolver::Impl
{
    /// The chemical system instance
    ChemicalSystem system;

    /// The partition of the chemical system
    Partition partition;

    /// The options for the smart equilibrium calculation
    SmartEquilibriumOptions options;

    /// The result of the last smart equilibrium calculation.
    SmartEquilibriumResult result;

    /// The solver for the equilibrium calculations
    EquilibriumSolver solver;

    /// The chemical properties of the chemical system
    ChemicalProperties properties;

    /// The amounts of the elements in the equilibrium partition
    Vector be;

    /// The solution of the equilibrium problem
    Vector n, y, z, x, u, r;

    /// Auxiliary vectors
    Vector dn, dy, dz;

    std::deque<Index> priority;

    std::deque<Index> ranking;

    //----------------------------------------------------------------------------------
    // TODO: move these vars into their corresponding fuctions
    /// The vector of amounts of equilibrium species
    Vector ne, ze, xe, ue, re;

    /// Auxiliary vectors
    Vector dne, dye, dze;

    /// The indices of the equilibrium and kinetic species
    Indices ies, iks;

    /// Auxiliary vectors delta(n) and delta(lna) in estimate function version v0
    Vector delta_ne, delta_lnae;
    //----------------------------------------------------------------------------------

    Vector a_dn;
    Vector a_dy;
    Vector a_dz;
    Vector a_n0;
    Vector a_y0;
    Vector a_z0;
    Vector a_n;
    Vector a_y;
    Vector a_z;
    Vector b_dn;
    Vector b_dy;
    Vector b_dz;
    Vector b_n0;
    Vector b_y0;
    Vector b_z0;
    Vector b_n;
    Vector b_y;
    Vector b_z;
    Vector a_u;
    Vector b_u;
    Vector a_x;
    Vector b_x;
    Vector a_r;
    Vector b_r;

    /// A class used to store the node of tree for smart equilibrium calculations.
    struct TreeNode
    {
        Vector be;
        Vector bebar;
        ChemicalState state;
        ChemicalProperties properties;
        EquilibriumSensitivity sensitivity;
        Matrix Mb;
        VectorXi imajor;
    };

    /// The tree used to save the calculated equilibrium states and respective sensitivities
    std::deque<TreeNode> tree;

    /// Construct a default SmartEquilibriumSolver::Impl instance.
    Impl()
    {}

    /// Construct an SmartEquilibriumSolver::Impl instance with given partition.
    Impl(const Partition& partition_)
    : partition(partition_), system(partition_.system()), solver(partition_)
    {
        setPartition(partition_);
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

    /// Set the partition of the chemical system.
    auto setPartition(const Partition& partition_) -> void
    {
        partition = partition_;

        // Initialize indices of the equilibrium and kinetic species
        ies = partition.indicesEquilibriumSpecies();
        iks = partition.indicesKineticSpecies();
        /*
        std::cout << "ies: ";
        for(auto i : ies)
            std::cout << i << ", ";
        std::cout << std::endl;

        std::cout << "iks: ";
        for(auto i : iks)
            std::cout << i << ", ";
        std::cout << std::endl;
        */
        solver.setPartition(partition_);
    }


    /// Learn how to perform a full equilibrium calculation (with tracking)
    auto learn(ChemicalState& state, double T, double P, VectorConstRef be) -> void
    {
        // Calculate the equilibrium state using conventional Gibbs energy minimization approach
        timeit( solver.solve(state, T, P, be),
                result.timing.learn_gibbs_energy_minimization= );

        // Store the result of the Gibbs energy minimization calculation performed during learning
        result.learning.gibbs_energy_minimization = solver.result();

        tic(0);

        // Adding into prioroity from the back ()
        priority.push_back(priority.size()); // adding to the priority deque according to the order of training in ODML algorithm
        ranking.push_back(0); // adding zero-rank of just added reference point

        // Update the chemical properties of the system
        properties = solver.properties();

        const auto RT = universalGasConstant * T;

        const auto& n = state.speciesAmounts();

        const auto u = properties.chemicalPotentials();
        const auto x = properties.moleFractions();

        const auto& dndb = solver.sensitivity().dndb;
        const auto& dudn = u.ddn;

        const Matrix dudb = dudn * dndb;

        const double nsum = sum(n);

        const auto& A = system.formulaMatrix();

        Canonicalizer canonicalizer(A);
        canonicalizer.updateWithPriorityWeights(n);

        const auto eps_n = options.amount_fraction_cutoff * nsum;
        const auto eps_x = options.mole_fraction_cutoff;

        const auto S = canonicalizer.S();

        auto imajorminor = canonicalizer.Q();
        const auto nummajor = std::partition(imajorminor.begin(), imajorminor.end(),
            [&](Index i) { return n[i] >= eps_n && x.val[i] >= eps_x; }) - imajorminor.begin();

        const auto imajor = imajorminor.head(nummajor);

        // const auto iprimary = canonicalizer.indicesBasicVariables();
        // const auto isecondary = canonicalizer.indicesNonBasicVariables();

        // Indices imajor;
        // imajor.reserve(isecondary.size());
        // for(auto i = 0; i < isecondary.size(); ++i)
        //     if(n[isecondary[i]] >= eps_n && x.val[isecondary[i]] >= eps_x)
        //         imajor.push_back(i);

        // VectorXi isecondary_major = isecondary(imajor);

        const Vector um = u.val(imajor);
        const auto dumdb = rows(dudb, imajor);

        const Matrix Mb = diag(inv(um)) * dumdb;

        // Store the computed solution into the knowledge tree
        tree.push_back({be, be/sum(be), state, properties, solver.sensitivity(), Mb, imajor});

        toc(0, result.timing.learning_storage);
    }

    /// Estimate the equilibrium state using sensitivity derivatives (profiling the expences)
    auto estimate(ChemicalState& state, double T, double P, VectorConstRef be) -> void
    {
        result.estimate.accepted = false;

        // Skip estimation if no previous full computation has been done
        if(tree.empty())
            return;

        // Relative and absolute tolerance parameters
        const auto reltol = options.reltol;
        const auto abstol = options.abstol;

        //---------------------------------------------------------------------------------------
        // Step 1: Search for the reference element (closest to the new state input conditions)
        //---------------------------------------------------------------------------------------
        tic(0);

        Vector rs;
        // Amounts of speceis
        Vector n, ns;
        // Variations of the potentials and species
        Vector du, dus, dns;
        Vector x;
        //Vector bebar = be/sum(be);
        Vector dbe;

        double nmin, nsum, uerror;
        Index inmin, iuerror;

        double nerror;
        Index inerror;


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

            const auto [success, error, ispecies] = pass_error_test(node);

            if(success)
            {
                const auto& be0 = node.be;
                const auto& n0 = node.state.speciesAmounts();
                const auto& dndb0 = node.sensitivity.dndb;

                n = n0 + dndb0 * (be - be0);

                nmin = n.minCoeff(&inmin);
                ntot = sum(n);

                // if(nmin/ntot < -options.amount_fraction_cutoff)
                if(nmin < -1.0e-5)
                    continue;

                toc(0, result.timing.estimate_search);

                ranking[*inode] += 1;

                auto comp = [&](Index l, Index r) { return ranking[l] > ranking[r]; };
                //std::sort(priority.begin(), priority.end(), comp);
                //std::stable_sort(priority.begin(), inode + 1, comp);
                if ( !((inode == priority.begin()) || (ranking[*inode_prev] >= ranking[*inode])) ) {
                    std::stable_sort(priority.begin(), inode + 1, comp);
                }

                state.setSpeciesAmounts(n);
                // state.setElementDualPotentials(y);
                // state.setSpeciesDualPotentials(z);

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
        return;
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
        const auto abstol = options.abstol;

        //---------------------------------------------------------------------------------------
        // Step 1: Search for the reference element (closest to the new state input conditions)
        //---------------------------------------------------------------------------------------
        tic(0);

        Vector rs;
        // Amounts of speceis
        Vector n, ns;
        // Variations of the potentials and species
        Vector du, dus, dns;
        Vector x;
        //Vector bebar = be/sum(be);
        Vector dbe;

        double nmin, nsum, uerror;
        Index inmin, iuerror;

        double nerror;
        Index inerror;


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

            const auto [success, error, ispecies] = pass_error_test(node);

            if(success)
            {
                const auto& be0 = node.be;
                const auto& n0 = node.state.speciesAmounts();
                const auto& y0 = node.state.elementDualPotentials();
                const auto& z0 = node.state.speciesDualPotentials();
                const auto& dndb0 = node.sensitivity.dndb;

                n = n0 + dndb0 * (be - be0);

                nmin = min(n(imajor));
                nsum = sum(n);

                const auto eps_n = options.amount_fraction_cutoff * nsum;

                if(nmin < -eps_n)
                    continue;

                toc(0, result.timing.estimate_search);

                ranking[*inode] += 1;

                auto comp = [&](Index l, Index r) { return ranking[l] > ranking[r]; };
                //std::sort(priority.begin(), priority.end(), comp);
                //std::stable_sort(priority.begin(), inode + 1, comp);
                if ( !((inode == priority.begin()) || (ranking[*inode_prev] >= ranking[*inode])) ) {
                    std::stable_sort(priority.begin(), inode + 1, comp);
                }

                state.setSpeciesAmounts(n);
                state.setElementDualPotentials(y0);
                state.setSpeciesDualPotentials(z0);

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
        return;
    }
    /// Estimate the equilibrium state using sensitivity derivatives,
    /// where search is performed by the NN algorithm and acceptance cirteria is using n and ln(a)
    auto estimate_nnsearch_acceptance_based_on_lna(ChemicalState& state, double T, double P, VectorConstRef be) -> void
    {
        // Skip estimation if no previous full computation has been done
        if(tree.empty())
            return;
        // Relative and absolute tolerance parameters
        const auto reltol = options.reltol;
        const auto abstol = options.abstol;

        //---------------------------------------------------------------------------------------
        // Step 1: Search for the reference element (closest to the new state input conditions)
        //---------------------------------------------------------------------------------------
        tic(0);

        // Comparison function based on the Euclidean distance
        auto distancefn = [&](const TreeNode& a, const TreeNode& b)
        {
            const auto& be_a = a.be;
            const auto& be_b = b.be;
            return (be_a - be).squaredNorm() < (be_b - be).squaredNorm();  // TODO: We need to extend this later with T and P contributions too (Allan, 26.06.2019)
        };

        // Find the entry with minimum "input" distance
        auto it = std::min_element(tree.begin(), tree.end(), distancefn);

        toc(0, result.timing.estimate_search);

        //----------------------------------------------------------------------------
        // Step 2: Calculate predicted state with a first-order Taylor approximation
        //----------------------------------------------------------------------------
        tic(1);

        // Get all the data stored in the reference element
        const auto& be0 = it->be;
        const auto& state0 = it->state;
        const auto& properties0 = it->properties;
        const auto& sensitivity0 = it->sensitivity;
        const auto& n0 = state0.speciesAmounts();
        const auto& N = ies.size() + iks.size();
        const auto& ne0 = n0(ies);
        const auto& nk0 = n0(iks);

        // Get the sensitivity derivatives dln(a) / dn (assuming all species are equilibrium species)
        MatrixConstRef dlna_dn = properties0.lnActivities().ddn;
        VectorConstRef lna0 = properties0.lnActivities().val;

        // Get the sensitivity derivatives w.r.t. the equilibrium species dln(ae) / dne
        MatrixConstRef dlnae_dne = dlna_dn(ies, ies);
        VectorConstRef lnae0 = lna0(ies);

        // Calculate perturbation of equilibrium species ne
        delta_ne.noalias() = sensitivity0.dndb * (be - be0); // delta(n) = dn/db * (b - b0)
        ne.noalias() = ne0 + delta_ne;                         // n = n0 + delta(n)
        delta_lnae.noalias() = dlnae_dne * delta_ne;             // delta(ln(a)) = d(lna)/dn * delta(n)

        toc(1, result.timing.estimate_mat_vec_mul);

        //----------------------------------------------
        // Step 3: Checking the acceptance criterion
        //----------------------------------------------
        tic(2);

        const auto& x = properties0.moleFractions().val;
        const auto& xe = x(ies);

        // Perform the check for the negative amounts
        const bool amount_check = ne.minCoeff() > -1e-5;

        // Check in the loop mole fractions, negative amount check, and variations of the ln(a)
        //bool amount_check = true;
        bool variation_check = true;
        const double fraction_tol = abstol * 1e-2;

        //const double fraction_tol = abstol;
        for(Index i = 0; i < ne.size(); ++i)
        {
            // If the fraction is too small, skip the variational check
            if(xe(i, 0) < fraction_tol)
                continue;

            // Perform the variational check
            if(std::abs(delta_lnae[i]) > abstol + reltol * std::abs(lnae0[i])) {
                variation_check = false;    // variation test failed
                result.estimate.failed_with_species = system.species(i).name();
                result.estimate.failed_with_amount = ne[i];
                result.estimate.failed_with_chemical_potential = lna0[i] + delta_lnae[i]; // TODO: change to the chemical potential

                break;
            }
        }

        toc(2, result.timing.estimate_acceptance);

        // Check if smart estimation failed with respect to variation of chemical potentials or amounts
        if(!variation_check || !amount_check){
            return;
        }

        // Set the output chemical state to the approximate amounts of species
        // Assign small values to all the amount in the interval [cutoff, 0] (instead of mirroring above)
        //for(unsigned int i = 0; i < ne.size(); ++i) if(ne[i] >= cutoff && ne[i] < 0) ne[i] = options.learning.epsilon;
        for(unsigned int i = 0; i < ne.size(); ++i) if(ne[i] < 0) ne[i] = options.learning.epsilon;

        // Update equilibrium species
        state.setSpeciesAmounts(ne, ies);

        // Set the estimate accepted status to true
        result.estimate.accepted = true;

        // Remember index of reference element used
        result.estimate.reference_state_index = std::distance(tree.begin(), it);

    }

    /// Estimate the equilibrium state using sensitivity derivatives (profiling the expences) v.1
    auto estimate_nnsearch_acceptance_based_on_residual(ChemicalState& state, double T, double P, VectorConstRef be) -> void
    {
        result.estimate.accepted = false;

        // Skip estimation if no previous full computation has been done
        if(tree.empty())
            return;

        MatrixConstRef Ae = partition.formulaMatrixEquilibriumPartition();

        //---------------------------------------------------------------------------------------
        // Step 1: Search for the reference element (closest to the new state input conditions)
        //---------------------------------------------------------------------------------------
        tic(0);

        // Comparison function based on the Euclidean distance (scaled)
        auto distancefn_scaled = [&](const TreeNode& a, const TreeNode& b)
        {
            Vector be_a = a.be/sum(a.be);
            Vector be_b = b.be/sum(b.be);
            Vector be_x = be/sum(be);

            return (be_a - be_x).squaredNorm() < (be_b - be_x).squaredNorm();  // TODO: We need to extend this later with T and P contributions too (Allan, 26.06.2019)
        };
        // Comparison function based on the Euclidean distance
        auto distancefn = [&](const TreeNode& a, const TreeNode& b)
        {
            const auto& be_a = a.be;
            const auto& be_b = b.be;
            return (be_a - be).squaredNorm() < (be_b - be).squaredNorm();  // TODO: We need to extend this later with T and P contributions too (Allan, 26.06.2019)
        };

        // Find the entry with minimum "input" distance
        auto it = std::min_element(tree.begin(), tree.end(), distancefn);

        toc(0, result.timing.estimate_search);

        //----------------------------------------------------------------------------
        // Step 2: Calculate predicted state with a first-order Taylor approximation
        //----------------------------------------------------------------------------
        tic(1);

        // Get all the data stored in the reference element
        const auto& be0 = it->be;
        const auto& state0 = it->state;
        const auto& properties0 = it->properties;
        const auto& sensitivity0 = it->sensitivity;
        const auto& T0 = state0.temperature();
        const auto& P0 = state0.pressure();
        const auto& n0 = state0.speciesAmounts();
        const auto& y0 = state0.elementDualPotentials();
        const auto& z0 = state0.speciesDualPotentials();
        const auto u0 = properties0.chemicalPotentials();
        const auto x0 = properties0.moleFractions();

        // Amounts of species related to equilibrium and kinetics species
        const auto& ne0 = n0(ies);
        const auto& ze0 = z0(ies);
        const auto& ue0 = u0.val(ies);
        const auto& xe0 = x0.val(ies);

        // Constant to scale the residual of the equilibrium species
        const auto RT = universalGasConstant*T0;

        // Calculate perturbation of n
        dn.noalias() = sensitivity0.dndb * (be - be0);
        dy.noalias() = sensitivity0.dydb * (be - be0);
        dz.noalias() = sensitivity0.dzdb * (be - be0);

        dne = dn;
        dze = dz(ies);

        y.noalias() = y0;
        z.noalias() = z0;
        //y.noalias() = y0 + dy * RT; // This approximation increases the number of learnings bu not improves the accuracy
        //z.noalias() = z0 + dz * RT; //

        ne.noalias() = ne0 + dne;
        ze.noalias() = ze0;
        //ze.noalias() = ze0 + dze * RT;

        // Negative cutoff control
        bool neg_amount_check = ne.minCoeff() > -1e-5;

        // If cutoff test didn't pass, estimation has failded
        if(neg_amount_check == false)
            return;

        for(auto i = 0; i < ne.size(); ++i)
            if(ne[i] <= 0.0)
                ne[i] = 1e-12;

        // Recompute the change in n
        dne = ne - ne0;

        // Calculate u(bar) = u(ref) + dudT(ref)*dT + dudP(ref)*dP + dudn(ref)*dn
        ue = ue0 + u0.ddn(ies, ies) *dne;

        // Calculate x(bar) = x(ref) + dxdn(ref)*dn
        xe = xe0 + x0.ddn(ies, ies) * dne;

        const auto& Aey = tr(Ae)*y;
        const auto& Aeye = Aey(ies);

        // Calculate the equilibrium residuals of the equilibrium species
        re = abs(ue - Aey(ies) - ze)/RT;  // TODO: We should actually collect the entries in u and z corresponding to equilibrium species

        // Eliminate species with mole fractions below cutoff from the residual analysis
        for(auto i = 0; i < ne.size(); ++i)
            if(xe[i] < options.mole_fraction_cutoff)
                re[i] = 0.0; // set their residuals to zero

        // Estimate the residual error of the trial Taylor approximation
        Index ispecies;
        const double error = re.maxCoeff(&ispecies);
        bool is_error_acceptable = error <= options.tol;

        toc(1, result.timing.estimate_mat_vec_mul);

        //----------------------------------------------
        // Step 3: Checking the acceptance criterion
        //----------------------------------------------
        tic(2);

        // Update three properties of the chemical state (for better initial guess of the GEM problem)
        state.setSpeciesAmounts(ne, ies);
        z(ies) = ze;
        state.setElementDualPotentials(y);
        state.setSpeciesDualPotentials(z);

        toc(2, result.timing.estimate_acceptance);

        // Check if smart estimation failed with respect to variation of chemical potentials or amounts
        if(is_error_acceptable == false)
            return;

        // Set the estimate accepted status to true
        result.estimate.accepted = true;

        // Update the chemical properties of the system
        properties = properties0;  // FIXME: We actually want to estimate properties = properties0 + variation : THIS IS A TEMPORARY SOLUTION!!!
    }

    /// Estimate the equilibrium state using sensitivity derivatives (profiling the expences) v.1
    auto estimate_nnsearch_acceptance_based_on_residual_no_quartz(ChemicalState& state, double T, double P, VectorConstRef be) -> void
    {
        result.estimate.accepted = false;

        // Skip estimation if no previous full computation has been done
        if(tree.empty())
            return;

        MatrixConstRef Ae = partition.formulaMatrixEquilibriumPartition();

        //---------------------------------------------------------------------------------------
        // Step 1: Search for the reference element (closest to the new state input conditions)
        //---------------------------------------------------------------------------------------
        tic(0);

        // Comparison function based on the Euclidean distance (scaled)
        auto distancefn_scaled = [&](const TreeNode& a, const TreeNode& b)
        {
            Vector be_a = a.be/sum(a.be);
            Vector be_b = b.be/sum(b.be);
            Vector be_x = be/sum(be);

            return (be_a - be_x).squaredNorm() < (be_b - be_x).squaredNorm();  // TODO: We need to extend this later with T and P contributions too (Allan, 26.06.2019)
        };
        // Comparison function based on the Euclidean distance
        auto distancefn = [&](const TreeNode& a, const TreeNode& b)
        {
            const auto& be_a = a.be;
            const auto& be_b = b.be;
            return (be_a - be).squaredNorm() < (be_b - be).squaredNorm();  // TODO: We need to extend this later with T and P contributions too (Allan, 26.06.2019)
        };

        // Find the entry with minimum "input" distance
        auto it = std::min_element(tree.begin(), tree.end(), distancefn);

        toc(0, result.timing.estimate_search);

        //----------------------------------------------------------------------------
        // Step 2: Calculate predicted state with a first-order Taylor approximation
        //----------------------------------------------------------------------------
        tic(1);

        // Get all the data stored in the reference element
        const auto& be0 = it->be;
        const auto& state0 = it->state;
        const auto& properties0 = it->properties;
        const auto& sensitivity0 = it->sensitivity;
        const auto& T0 = state0.temperature();
        const auto& P0 = state0.pressure();
        const auto& n0 = state0.speciesAmounts();
        const auto& y0 = state0.elementDualPotentials();
        const auto& z0 = state0.speciesDualPotentials();
        const auto u0 = properties0.chemicalPotentials();
        const auto x0 = properties0.moleFractions();

        // Amounts of species related to equilibrium and kinetics species
        const auto& ne0 = n0(ies);
        const auto& ze0 = z0(ies);
        const auto& ue0 = u0.val(ies);
        const auto& xe0 = x0.val(ies);

        // Constant to scale the residual of the equilibrium species
        const auto RT = universalGasConstant*T0;

        // Calculate perturbation of n
        dn.noalias() = sensitivity0.dndb * (be - be0);
        dy.noalias() = sensitivity0.dydb * (be - be0);
        dz.noalias() = sensitivity0.dzdb * (be - be0);

        dne = dn;
        dze = dz(ies);

        y.noalias() = y0;
        z.noalias() = z0;

        ne.noalias() = ne0 + dne;
        ze.noalias() = ze0;

        // Negative cutoff control
        bool neg_amount_check = ne.minCoeff() > -1e-5;

         // If cutoff test didn't pass, estimation has failded
        if(neg_amount_check == false)
            return;

        for(auto i = 0; i < ne.size(); ++i)
            if(ne[i] <= 0.0)
                ne[i] = 1e-12;

        // Recompute the change in n309.593
        dne = ne - ne0;

        // Calculate u(bar) = u(ref) + dudT(ref)*dT + dudP(ref)*dP + dudn(ref)*dn
        u = u0.val + u0.ddn * dn;
        ue = ue0 + u0.ddn(ies, ies) * dne;

        // Calculate x(bar) = x(ref) + dxdn(ref)*dn
        xe = xe0 + x0.ddn(ies, ies) * dne;

        const auto& Aey = tr(Ae)*y;
        const auto& Aeye = Aey(ies);

        // Calculate the equilibrium residuals of the equilibrium species
        re = abs(ue - Aey(ies) - ze)/RT;  // TODO: We should actually collect the entries in u and z corresponding to equilibrium species
        // TODO: fix this exclusion of quartz in smarter
        Vector re_no_quartz(re.size()-1);
        for(auto i = 0; i < ne.size()-1; ++i)
            re_no_quartz[i] = re[i];

        // Eliminate species with mole fractions below cutoff from the residual analysis
        for(auto i = 0; i < ne.size(); ++i)
            if(xe[i] < options.mole_fraction_cutoff)
                re[i] = 0.0; // set their residuals to zero

        // Estimate the residual error of the trial Taylor approximation
        Index ispecies;
        const double error = re.maxCoeff(&ispecies);
        const double error_no_quartz = re_no_quartz.maxCoeff(&ispecies);

        bool is_error_acceptable = error_no_quartz <= options.tol;

        toc(1, result.timing.estimate_mat_vec_mul);

        //----------------------------------------------
        // Step 3: Checking the acceptance criterion
        //----------------------------------------------
        tic(2);

        // Update three properties of the chemical state (for better initial guess of the GEM problem)
        state.setSpeciesAmounts(ne, ies);
        z(ies) = ze;
        state.setElementDualPotentials(y);
        state.setSpeciesDualPotentials(z);

        toc(2, result.timing.estimate_acceptance);

        // Check if smart estimation failed with respect to variation of chemical potentials or amounts
        if(is_error_acceptable == false)
            return;

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
        setPartition(problem.partition());
        const auto T = problem.temperature();
        const auto P = problem.pressure();
        const auto& iee = partition.indicesEquilibriumElements(); // This statement needs to be after setPartition
        be = problem.elementAmounts()(iee);
        return solve(state, T, P, be);
    }

    auto solve(ChemicalState& state, double T, double P, VectorConstRef be) -> SmartEquilibriumResult
    {
        tic(0);

        // Absolutely ensure an exact Hessian of the Gibbs energy function is used in the calculations
        setOptions(options);

        // Reset the result of the last smart equilibrium calculation
        result = {};

        // Perform a smart estimate of the chemical state
        //timeit( estimate_nnsearch_acceptance_based_on_lna(state, T, P, be),
        //        result.timing.estimate= );
        //timeit( estimate_nnsearch_acceptance_based_on_residual(state, T, P, be),
        //         result.timing.estimate= );
        timeit(estimate_nnsearch_acceptance_based_on_residual_no_quartz(state, T, P, be),
              result.timing.estimate= );

        // Perform a learning step if the smart prediction is not sactisfatory
        if(!result.estimate.accepted)
            timeit( learn(state, T, P, be), result.timing.learn= );

        toc(0, result.timing.solve);

        return result;
    }

    auto solve(ChemicalState& state, double T, double P, VectorConstRef be, Index istep, Index icell) -> SmartEquilibriumResult
    {
        tic(SOLVE_STEP);

        // Absolutely ensure an exact Hessian of the Gibbs energy function is used in the calculations
        setOptions(options);

        // Reset the result of the last smart equilibrium calculation
        result = {};

        // Perform a smart estimate of the chemical state
        timeit(estimate_nnsearch_acceptance_based_on_residual(state, T, P, be),
                result.timing.estimate= );

        // Perform a learning step if the smart prediction is not sactisfatory
        if(!result.estimate.accepted)
            timeit( learn(state, T, P, be), result.timing.learn= );

        toc(SOLVE_STEP, result.timing.solve);

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
        return result;
    }

};

SmartEquilibriumSolver::SmartEquilibriumSolver()
: pimpl(new Impl())
{}

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

SmartEquilibriumSolver::~SmartEquilibriumSolver()
{}

auto SmartEquilibriumSolver::setOptions(const SmartEquilibriumOptions& options) -> void
{
    pimpl->setOptions(options);
}

auto SmartEquilibriumSolver::setPartition(const Partition& partition) -> void
{
    pimpl->setPartition(partition);
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, double T, double P, VectorConstRef be) -> SmartEquilibriumResult
{
    return pimpl->solve(state, T, P, be);
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, double T, double P, VectorConstRef be, Index istep, Index icell) -> SmartEquilibriumResult
{
    return pimpl->solve(state, T, P, be, istep, icell);
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

} // namespace Reaktoro

