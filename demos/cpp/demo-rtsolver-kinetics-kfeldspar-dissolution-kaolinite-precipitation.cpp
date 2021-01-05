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

// C++ includes
#include <sstream>      // for using stringstream
#include <iomanip>      // for setprecition

#if defined _WIN32      // for creating a new folder
#include <windows.h>
#ifdef __MINGW32__
#include <sys/stat.h>
#endif
#else
#include <sys/stat.h>
#endif

// Reaktoro includes
#include <Reaktoro/Reaktoro.hpp>

using namespace Reaktoro;

struct Params
{
    // Discretization params
    int ncells = 0; // the number of cells in the spacial discretization
    int nsteps = 0; // the number of steps in the reactive transport simulation
    double xl = 0.0; // the x-coordinates of the left boundaries
    double xr = 0.0; // the x-coordinates of the right boundaries
    double dx = 0.0; // the space step (in units of m)
    double dt = 0.0; // the time step (in units of s)

    // Physical params
    double D = 0.0; // the diffusion coefficient (in units of m2/s)
    double v = 0.0; // the Darcy velocity (in units of m/s)
    double T = 0.0; // the temperature (in units of degC)
    double P = 0.0; // the pressure (in units of bar)

    // Kinetic and equilibrium solvers' parameters
    bool use_smart_equilibrium_solver = false;
    bool use_smart_kinetics_solver = false;

    double smart_equilibrium_reltol = 0.0;
    double smart_equilibrium_abstol = 0.0;
    double smart_equilibrium_cutoff = 0.0;
    double smart_equilibrium_tol = 0.0;

    double smart_kinetics_reltol = 0.0;
    double smart_kinetics_abstol = 0.0;
    double smart_kinetics_tol = 0.0;

    bool output_results = true;

    GibbsHessian hessian = GibbsHessian::Exact;

    std::string activity_model = "";

};


struct RTKineticsResults
{
     // Conventional kinetic and conventional equilibrium schemes' times
    // *********************************************************************************//

    /// Total CPU time (in s) required by conventional kinetic and equilibrium schemes.
    double conv_kin_conv_eq_total = 0.0;

    /// Total CPU time (in s) required by conventional kinetic and equilibrium schemes
    /// excluding the costs for the chemical properties evaluation.
    double conv_kin_conv_eq_total_ideal_properties = 0.0;

     // Smart kinetic and conventional equilibrium schemes' times
    // *********************************************************************************//

    /// Total CPU time (in s) required by smart kinetic and conventional equilibrium schemes.
    double smart_kin_conv_eq_total = 0.0;

    /// Total CPU time (in s) required by smart kinetic and conventional equilibrium schemes
    /// excluding the costs for the search of the closest reference states.
    double smart_kin_conv_eq_total_ideal_search = 0.0;

    /// Total CPU time (in s) required by smart kinetic and conventional equilibrium schemes
    /// excluding the costs for the search and storage of the closest reference states.
    double smart_kin_conv_eq_total_ideal_search_store = 0.0;

    /// Total CPU time (in s) required by smart kinetic and conventional equilibrium schemes
    /// excluding the costs for the search and storage of the closest reference states
    /// and the chemical properties evaluation.
    double smart_kin_conv_eq_total_ideal_search_store_properties = 0.0;

     // Conventional kinetic and smart equilibrium schemes' times
    // *********************************************************************************//

    double conv_kin_smart_eq_total = 0.0;

    /// Total CPU time (in s) required by conventional kinetic and smart equilibrium schemes
    /// excluding the costs for the search of the closest reference states.
    double conv_kin_smart_eq_total_ideal_search = 0.0;

    /// Total CPU time (in s) required by conventional kinetic and smart equilibrium schemes
    /// excluding the costs for the search and store of the closest reference states.
    double conv_kin_smart_eq_total_ideal_search_store = 0.0;

    /// Total CPU time (in s) required by conventional kinetic and smart equilibrium schemes
    /// excluding the costs for the search and store of the closest reference states
    /// and the chemical properties evaluation.
    double conv_kin_smart_eq_total_ideal_search_store_properties = 0.0;


    // Smart kinetic and smart equilibrium schemes' times
    // *********************************************************************************//

    double smart_kin_smart_eq_total = 0.0;

    /// Total CPU time (in s) required by conventional kinetic and smart equilibrium schemes
    /// excluding the costs for the search of the closest reference states.
    double smart_kin_smart_eq_total_ideal_search = 0.0;

    /// Total CPU time (in s) required by conventional kinetic and smart equilibrium schemes
    /// excluding the costs for the search and store of the closest reference states.
    double smart_kin_smart_eq_total_ideal_search_store = 0.0;

    /// Total CPU time (in s) required by conventional kinetic and smart equilibrium schemes
    /// excluding the costs for the search and store of the closest reference states
    /// and the chemical properties evaluation.
    double smart_kin_smart_eq_total_ideal_search_store_properties = 0.0;

    // Total times of reactive transport steps
    // *********************************************************************************//

    /// Total time taken to perform all time steps using conventional kinetics and conventional equilibrium algorithm
    double time_reactive_transport_conv_kin_conv_eq = 0.0;

    /// Total time taken to perform all time steps using smart kinetics and conventional equilibrium algorithm
    double time_reactive_transport_smart_kin_conv_eq = 0.0;

    /// Total time taken to perform all time steps using conventional kinetics and smart equilibrium algorithm
    double time_reactive_transport_conv_kin_smart_eq = 0.0;

    /// Total time taken to perform all time steps using smart kinetics and smart equilibrium algorithm
    double time_reactive_transport_smart_kin_smart_eq = 0.0;

    // Total times of reactive transport steps
    // *********************************************************************************//

    /// Accumulated timing information of all equilibrium calculations.
    EquilibriumTiming equilibrium_timing;

    /// The accumulated timing information of all smart equilibrium calculations.
    SmartEquilibriumTiming smart_equilibrium_timing;

    /// Accumulated timing information of all kinetic calculations.
    KineticTiming kinetic_timing;

    /// The accumulated timing information of all smart kinetic calculations.
    SmartKineticTiming smart_kinetic_timing;

};

/// Forward declaration
auto mkdir(const std::string& folder) -> bool;
auto outputConsole(const Params& params) -> void;
auto makeResultsFolder(const Params& params) -> std::string;
auto runReactiveTransport(const Params& params, RTKineticsResults& results) -> void;

using MineralCatalystFunction = std::function<ChemicalScalar(const ChemicalProperties&)>;
auto mineralCatalystFunctionActivity(const MineralCatalyst& catalyst, const ChemicalSystem& system) -> MineralCatalystFunction
{
    const double power = catalyst.power;
    const std::string species = catalyst.species;
    const Index ispecies = system.indexSpeciesWithError(species);

    MineralCatalystFunction fn = [=](const ChemicalProperties& properties) mutable
    {
        const ChemicalVector& ln_a = properties.lnActivities();
        ChemicalScalar ai = exp(ln_a[ispecies]);
        ChemicalScalar res = pow(ai, power);
        return res;
    };

    return fn;
}
auto mineralCatalystFunctionPartialPressure(const MineralCatalyst& catalyst, const ChemicalSystem& system) -> MineralCatalystFunction
{
    const auto gas         = catalyst.species;                           // the species of the catalyst
    const auto power       = catalyst.power;                             // the power of the catalyst
    const auto idx_phase   = system.indexPhase("Gaseous");               // the index of the gaseous phase
    const auto gases       = names(system.phase(idx_phase).species());   // the names of the gaseous species
    const auto igases      = system.indicesSpecies(gases);               // the indices of the gaseous species
    const auto ifirst      = system.indexFirstSpeciesInPhase(idx_phase); // the index of the first gaseous species
    const auto igas        = index(gas, gases);                          // the index of the gaseous species
    const auto num_gases   = gases.size();                               // the number of gases

    ChemicalScalar res;

    MineralCatalystFunction fn = [=](const ChemicalProperties& properties) mutable
    {
        // The pressure and composition of the system
        const auto P = properties.pressure();
        const auto n = properties.composition();

        // The molar composition of the gaseous species
        const auto ng = rows(n, ifirst, num_gases);

        // The total number of moles in the gaseous phase
        const auto ngsum = sum(ng);

        // The mole fraction of the gas
        const auto xi = ng[igas]/ngsum;

        // The pressure in units of bar
        const auto Pbar = convertPascalToBar(P);

        // Evaluate the mineral catalyst function
        res = pow(xi * Pbar, power);

        return res;
    };

    return fn;
}
auto mineralCatalystFunction(const MineralCatalyst& catalyst, const ChemicalSystem& system) -> MineralCatalystFunction
{
    if (catalyst.quantity == "a" || catalyst.quantity == "activity")
        return mineralCatalystFunctionActivity(catalyst, system);
    else
        return mineralCatalystFunctionPartialPressure(catalyst, system);
}
auto mineralMechanismFunction(const MineralMechanism& mechanism, const Reaction& reaction, const ChemicalSystem& system) -> ReactionRateFunction
{
    // The number of chemical species in the system
    const unsigned num_species = system.numSpecies();

    // The universal gas constant (in units of kJ/(mol*K))
    const double R = 8.3144621e-3;

    // Create the mineral catalyst functions
    std::vector<MineralCatalystFunction> catalysts;
    for(const MineralCatalyst& catalyst : mechanism.catalysts)
        catalysts.push_back(mineralCatalystFunction(catalyst, system));

    // Auxiliary variables
    ChemicalScalar aux, f, g;

    // Define the mineral mechanism function
    ReactionRateFunction fn = [=](const ChemicalProperties& properties) mutable
    {
        // The temperature and pressure of the system
        const Temperature T = properties.temperature();

        // The result of this function evaluation
        ChemicalScalar res(num_species);

        // Calculate the saturation index of the mineral
        const auto lnK = reaction.lnEquilibriumConstant(properties);
        const auto lnQ = reaction.lnReactionQuotient(properties);
        const auto lnOmega = lnQ - lnK;

        // Calculate the rate constant for the current mechanism
        const auto kappa = mechanism.kappa * exp(-mechanism.Ea/R * (1.0/T - 1.0/298.15));

        // Calculate the saturation index
        const auto Omega = exp(lnOmega);

        // Calculate the p and q powers of the saturation index Omega
        const auto pOmega = pow(Omega, mechanism.p);
        const auto qOmega = pow(1 - pOmega, mechanism.q);

        // Calculate the function f
        f = kappa * qOmega;

        // Calculate the function g
        g = ChemicalScalar(num_species, 1.0);

        for(const MineralCatalystFunction& catalyst : catalysts)
            g *= catalyst(properties);

        // Calculate the resulting mechanism function
        res = f * g;

        return res;
    };

    return fn;
}

int main()
{
    Time start = time();

    // Step 1: Initialise auxiliary time-related constants
    int second = 1;
    int minute = 60;
    int hour = 60 * minute;
    int day = 24 * hour;
    int week = 7 * day;
    // int month = 30 * day;
    // int year = 365 * day;

    // Step 2: Define parameters for the reactive transport simulation
    Params params = {};

    // Define discretization parameters
    params.xl = 0.0; // the x-coordinates of the left boundaries
    /*
    params.xr = 0.1; // the x-coordinates of the right boundaries
    params.ncells = 10; // the number of cells in the spacial discretization
    */
    ///*
    params.xr = 50.0; // the x-coordinates of the right boundaries
    params.ncells = 50; // the number of cells in the spacial discretization
    //*/
    params.nsteps = 10; // the number of steps in the reactive transport simulation
    params.dx = (params.xr - params.xl) / params.ncells; // the time step (in units of s)
    params.dt = 0.01 * second; // the time step (in units of s)

    // Define physical and chemical parameters
    params.D = 5.0e-3;     // the diffusion coefficient (in units of m2/s)
    params.v = 10.0;       // the Darcy velocity (in units of m/s)
    params.T = 25.0;       // the temperature (in units of degC)
    params.P = 1.0;          // the pressure (in units of bar)

    // Define parameters of the equilibrium solvers
    params.smart_equilibrium_reltol = 1e-1;
    params.smart_equilibrium_abstol = 1e-12;

    // Define parameters of the kinetics solvers
    params.smart_kinetics_reltol = 1e-1;
    params.smart_kinetics_abstol = 1e-2;
    params.smart_kinetics_tol = 5e-1; // for the residual check

    params.activity_model = "hkf-full";

    std::cout << "total time                          : " << elapsed(start) << std::endl;

    // Output
    outputConsole(params);

    // RTKineticsResults
    RTKineticsResults results;

    // Execute reactive transport with different solvers
    params.use_smart_kinetics_solver = false; params.use_smart_equilibrium_solver = false; runReactiveTransport(params, results);

    results.conv_kin_conv_eq_total = results.kinetic_timing.solve;
    results.conv_kin_conv_eq_total_ideal_properties = results.kinetic_timing.solve - results.kinetic_timing.integrate_chemical_properties;

    std::cout << "total time                          : " << elapsed(start) << std::endl;

    return 0;
}
auto runReactiveTransport(const Params& params, RTKineticsResults& results) -> void
{

    double molar_mass_kfeldspar = 278.3315; // g / mol
    double molar_mass_kaolinite = 258.071; // g / mol // or 258.1604

    // Step **: Create the results folder
    auto folder = makeResultsFolder(params);

    // Step **: Define chemical equilibrium solver options
    EquilibriumOptions equilibrium_options;

    // Step **: Define smart chemical equilibrium solver options
    SmartEquilibriumOptions smart_equilibrium_options;
    smart_equilibrium_options.reltol = params.smart_equilibrium_reltol;

    // Step **: Define chemical kinetic solver options
    KineticOptions kinetic_options;
    kinetic_options.equilibrium = equilibrium_options;
    kinetic_options.use_smart_equilibrium_solver = params.use_smart_equilibrium_solver;

    // Step **: Define smart chemical kinetic solver options
    SmartKineticOptions smart_kinetic_options;
    smart_kinetic_options.reltol = params.smart_kinetics_reltol;
    smart_kinetic_options.abstol = params.smart_kinetics_abstol;
    smart_kinetic_options.tol = params.smart_kinetics_tol;
    smart_kinetic_options.learning = kinetic_options;
    smart_kinetic_options.learning.equilibrium = equilibrium_options;

    // Step **: Construct the chemical system with its phases and species (using ChemicalEditor)
    ChemicalEditor editor;
    std::string elements = "H Cl S O Ca Na K Mg C Si Al";

    if(params.activity_model == "hkf-full")
        editor.addAqueousPhaseWithElements(elements);
    else if(params.activity_model == "pitzer-full")
        editor.addAqueousPhaseWithElements(elements)
            .setChemicalModelPitzerHMW()
            .setActivityModelDrummondCO2();
    else if(params.activity_model == "dh-full")
        editor.addAqueousPhaseWithElements(elements)
            .setChemicalModelDebyeHuckel()
            .setActivityModelDrummondCO2();

    std::vector<std::string> minerals = {"K-Feldspar", "Kaolinite"};
    //std::vector<std::string> minerals = {"K-Feldspar", "Kaolinite", "Gibbsite", "Quartz", "Muscovite", "Halloysite"};
    for(const auto& mineral : minerals)
        editor.addMineralPhase(mineral);

    // Step **: Create the ChemicalSystem object using the configured editor
    ChemicalSystem system(editor);

    // ------------------------------------------------------------------------------------------
    // Definition of K-feldspar(K(AlSi3)O8) reaction
    // ------------------------------------------------------------------------------------------

    double ssa_m2mol_kfeldspar = 10.0; // defined in PHREEQC, m2/mol
    double ssa_m2kg_kfeldspar = ssa_m2mol_kfeldspar / molar_mass_kfeldspar / 1e-3;
    double mass_kfeldspar = 10 * 1e-3; // 10 g
    double sa_m2_kfeldspar = ssa_m2kg_kfeldspar * mass_kfeldspar;
    double vsa_m2m3_kfeldspar = 1000.0;

    std::string eq_str_kfeldspar = "K-Feldspar + 4*H+ = Al+++ + K+ + 3*SiO2(aq) + 2*H2O(l)";
//    ReactionEquation eq_kfeldspar(eq_str_kfeldspar);
//    Reaction kfeldspar_reaction(eq_kfeldspar, system);
//    kfeldspar_reaction.setName("K-Feldspar dissolution");
//    kfeldspar_reaction.setRate(kfeldspar_rate_func);
    MineralReaction min_react_kfeldspar("K-Feldspar");
    min_react_kfeldspar.setEquation(eq_str_kfeldspar)
            .addMechanism("logk = -11.5 mol/(m2*s); Ea = 38.2 kJ/mol; p = 1.0; q = 1.0")
            .setSurfaceArea(sa_m2_kfeldspar, "m2");

    // ------------------------------------------------------------------------------------------
    // Definition of Kaolinite(Al2Si2O5(OH)4) reaction
    // ------------------------------------------------------------------------------------------
    double ssa_m2mol_kaolinite = 1.0; // defined in PHREEQC, m2/mol
    double ssa_m2kg_kaolinite = ssa_m2mol_kaolinite / molar_mass_kaolinite / 1e-3;
    double mass_kaolinite = 10 * 1e-3; // 10 g
    double sa_m2_kaolinite = ssa_m2kg_kaolinite * mass_kaolinite;
    double vsa_m2m3_kaolinite = 1000.0;

    std::string eq_str_kaolinite = "2*Al+++ + 2*SiO2(aq) + 5*H2O(l) = Kaolinite + 6*H+";
    //std::string kaolinite_eq_str = "Kaolinite + 6*H+ = 2*Al+++ + 2*SiO2(aq) + 5*H2O(l)";
    //ReactionEquation kaolinite_equation(kaolinite_eq_str);
    MineralReaction min_react_kaolinite("Kaolinite");
    min_react_kaolinite.setEquation(eq_str_kaolinite)
            .addMechanism("logk = -13.2 mol/(m2*s); Ea = 71.0 kJ/mol; p = 1/9; q = 2.0")
            .setSurfaceArea(sa_m2_kaolinite, "m2");

    editor.addMineralReaction(min_react_kfeldspar);
    editor.addMineralReaction(min_react_kaolinite);

    // Step **: Create the ReactionSystem instances
    ReactionSystem reactions(editor);
    //ReactionSystem reacts(system, {react_kfeldspar, react_kaolinite});

    //std::cout << "system = \n" << system << std:: endl;
    std::cout << "reactions number = " << reactions.numReactions() << std:: endl;

    // Step **: Create the ReactionSystem instances
    Partition partition(system);
    std::vector<std::string> kin_species = {"Kaolinite", "K-Feldspar"};
    partition.setKineticSpecies(kin_species);

    std::cout << "num kinetic species = " << partition.numKineticSpecies() << std:: endl;

    // Step **: Define the initial condition (IC) of the reactive transport modeling problem
    EquilibriumProblem problem_ic(system);
    problem_ic.setTemperature(params.T, "celsius");
    problem_ic.setPressure(params.P, "bar");
    problem_ic.add("H2O",   1.0, "kg");
    problem_ic.add("K", 1e-3, "mol");
    problem_ic.add("NaCl", 1e-3, "mol");
    problem_ic.add("SiO2", 1e-3, "mol");
    problem_ic.add("CO2", 1e-4, "mol");

    // Step **: Calculate the equilibrium states for the IC and BC
    ChemicalState state_ic = equilibrate(problem_ic);
    state_ic.setSpeciesAmount("K-Feldspar", 10.0, "mol");
    state_ic.setSpeciesAmount("Kaolinite", 10.0, "mol");

    // Step **: Define the boundary condition (BC, rainfall)
    EquilibriumInverseProblem problem_bc(system);
    problem_bc.setTemperature(params.T, "celsius");
    problem_bc.setPressure(params.P, "bar");
    problem_bc.add("H2O",   1.0, "kg");
    problem_bc.add("K",  10e-6, "mol");
    problem_bc.add("Al", 10-8, "mol");
    problem_bc.add("SiO2", 10e-6, "mol");
    problem_bc.add("Cl",   std::pow(10, -4.98), "mol");
    problem_bc.pH(5.0);

    ChemicalState state_bc = equilibrate(problem_bc);

    // Step **: Scale the boundary condition state
    state_bc.scaleVolume(1.0, "m3");

    // Step **: Scale the volumes of the phases in the initial condition
    state_ic.scaleVolume(1.0, "m3");
    state_ic.scalePhaseVolume("Aqueous", 0.1, "m3");    // 10% if the 1.0m3
    state_ic.scalePhaseVolume("K-Feldspar", 0.882, "m3");   // 0.882 = 0.98 * 0.9 (0.9 is 90% of 1.0m3, 0.98 is 98% K-Feldspar of the rock)
    state_ic.scalePhaseVolume("Kaolinite", 0.018, "m3");  // 0.882 = 0.02 * 0.9 (0.9 is 90% of 1.0m3, 0.02 is 2% Kaolinite of the rock)

    std::cout << "state_ic = \n" << state_ic << std:: endl;
    getchar();
    std::cout << "state_bc = \n" << state_bc << std:: endl;
    getchar();

    // Reaction function (with provided surface area )
    ReactionRateFunction kfeldspar_rate_func = [&min_react_kfeldspar, &system](const ChemicalProperties& properties) -> ChemicalScalar {

        // The number of chemical species in the system
        const unsigned num_species = system.numSpecies();

        // The mineral reaction rate using specified surface area
        ChemicalScalar res(num_species);

        // The universal gas constant (in units of kJ/(mol*K))
        const double R = 8.3144621e-3;

        // The temperature and pressure of the system
        const Temperature T = properties.temperature();

        // Create a Reaction instance
        Reaction reaction(min_react_kfeldspar.equation(), system);

        for(const MineralMechanism& mechanism : min_react_kfeldspar.mechanisms()){
            // Create the mineral catalyst functions
            std::vector<MineralCatalystFunction> catalysts;
            for(const MineralCatalyst& catalyst : mechanism.catalysts)
                catalysts.push_back(mineralCatalystFunction(catalyst, system));\

            // Auxiliary variables
            ChemicalScalar aux, f, g;

            // Calculate the saturation index of the mineral
            const auto lnK = reaction.lnEquilibriumConstant(properties);
            const auto lnQ = reaction.lnReactionQuotient(properties);
            const auto lnOmega = lnQ - lnK;

            // Calculate the rate constant for the current mechanism (Arrhenius equation)
            const auto kappa = mechanism.kappa * exp(-mechanism.Ea/R * (1.0/T - 1.0/298.15));

            // Calculate the saturation index
            const auto Omega = exp(lnOmega);

            // Calculate the p and q powers of the saturation index Omega
            const auto pOmega = pow(Omega, mechanism.p);
            const auto qOmega = pow(1 - pOmega, mechanism.q);

            std::cout << "lnK = " << lnK << std:: endl;
            std::cout << "lnQ = " << lnQ << std:: endl;
            std::cout << "lnOmega = " << lnOmega << std:: endl;
            std::cout << "Omega = " << Omega << std:: endl;
            std::cout << "p = " << mechanism.p << std:: endl;
            std::cout << "q = " << mechanism.q << std:: endl;
            std::cout << "pOmega = " << pOmega << std:: endl;
            std::cout << "qOmega = " << qOmega << std:: endl;
            std::cout << "kappa = " << kappa << std:: endl;


            // Calculate the function f
            f = kappa * qOmega;

            // Calculate the function g
            g = ChemicalScalar(num_species, 1.0);

            // Calculate the resulting mechanism function
            res = f * g;

            // Multiply the mechanism contributions by the surface area of the mineral
            res *= min_react_kfeldspar.surfaceArea();
            std::cout << "kappa = " << kappa << std:: endl;


        }
        std::cout << "res = " << res << std:: endl;

        return res;

//        double logk = -11.5;
//        double k = std::pow(std::exp(1), logk);
//        double Ea = 71.0;
//        double m = 1.0;
//        double n = 1.0;
//        // Rate = A * k * (Omega^m - 1)^n
//        // k is the growth (dissolution) rate constant
//        // m, n > 0
//        // Omega^m is saturation ration
//        // Rate > 0, precipitation (solution is supersaturated w.r.t. mineral)
//        // Rate < 0, dissolution (solution is undersaturated w.r.t. mineral)
//        // A is the reactive surface area (not constant), A = ? * specific surface area
//        // A ~ total surface area

    };

    ReactionRateFunction kaolinite_rate_func = [&min_react_kaolinite, &system](const ChemicalProperties& properties) -> ChemicalScalar {

        // The number of chemical species in the system
        const unsigned num_species = system.numSpecies();

        // The mineral reaction rate using specified surface area
        ChemicalScalar res(num_species);

        // The universal gas constant (in units of kJ/(mol*K))
        const double R = 8.3144621e-3;

        // The temperature and pressure of the system
        const Temperature T = properties.temperature();

        // Create a Reaction instance
        Reaction reaction(min_react_kaolinite.equation(), system);

        for(const MineralMechanism& mechanism : min_react_kaolinite.mechanisms()){
            // Create the mineral catalyst functions
            std::vector<MineralCatalystFunction> catalysts;
            for(const MineralCatalyst& catalyst : mechanism.catalysts)
                catalysts.push_back(mineralCatalystFunction(catalyst, system));\

            // Auxiliary variables
            ChemicalScalar aux, f, g;

            // Calculate the saturation index of the mineral
            const auto lnK = reaction.lnEquilibriumConstant(properties);
            const auto lnQ = reaction.lnReactionQuotient(properties);
            const auto lnOmega = lnQ - lnK;

            // Calculate the rate constant for the current mechanism
            const auto kappa = mechanism.kappa * exp(-mechanism.Ea/R * (1.0/T - 1.0/298.15));

            // Calculate the saturation index
            const auto Omega = exp(lnOmega);

            // Calculate the p and q powers of the saturation index Omega
            const auto pOmega = pow(Omega, mechanism.p);
            const auto qOmega = pow(1 - pOmega, mechanism.q);

            // Calculate the function f
            f = kappa * qOmega;

            // Calculate the function g
            g = ChemicalScalar(num_species, 1.0);

            // Calculate the resulting mechanism function
            res = f * g;

            // Multiply the mechanism contributions by the surface area of the mineral
            res *= min_react_kaolinite.surfaceArea();

        }
        return res;
//        ChemicalScalar val;
//        double A = 0.0; // reactive surface area of the mineral [m2 / m3_total_rock]
//        double k_kaolinite = 0.0; // rate constant [mol / m2 / s]
//        double Omega = 0.0; // function of saturation state, f(Omega)
//        double logk = -13.2;
//        double k = std::pow(std::exp(1), -logk);
//
//        double Ea = 71.0;
//        double m = 1/9;
//        double n = 2.0;
        // Omega = IAP / K_eq
        // IAP ion activity product (written for the dissolution reaction)
        // K_eq equilibrium constant

        // Rate = A * k * (Omega^m - 1)^n
        // k is the growth (dissolution) rate constant
        // m, n > 0
        // Omega^m is saturation ration
        //
        // Rate > 0, precipitation (solution is supersaturated w.r.t. mineral)
        // Rate < 0, dissolution (solution is undersaturated w.r.t. mineral)
        //  A * k_kaolinite * std::pow((std::pow(Omega, m) - 1), n)
    };

    ChemicalProperties props = state_ic.properties();
    auto evaluate_pH = ChemicalProperty::pH(system);
    auto evaluate_alkalinity = ChemicalProperty::alkalinity(system);
    auto pH = evaluate_pH(props);

    std::cout << "pH  = " << pH.val << std::endl;

    ChemicalScalar kfeldspar_r = kfeldspar_rate_func(props);
    std::cout << "kfeldspar_r  = " << kfeldspar_r.val << std::endl;
    std::cout << "kfeldspar_r_ddn  = " << kfeldspar_r.ddn << std::endl;
    getchar();

    ChemicalScalar kaolinite_r = kaolinite_rate_func(props);
    std::cout << "kaolinite_r  = " << kaolinite_r.val << std::endl;
    std::cout << "kaolinite_r_ddn  = " << kaolinite_r.ddn << std::endl;
    getchar();

    // Step **: Create the mesh for the column
    Mesh mesh(params.ncells, params.xl, params.xr);

    // Step **: Create a chemical field object with every cell having state given by state_ic
    ChemicalField field(mesh.numCells(), state_ic);

    // Step **: Define the options for the reactive transport solver
    ReactiveTransportOptions reactive_transport_options;
    reactive_transport_options.use_smart_equilibrium_solver = params.use_smart_equilibrium_solver;
    reactive_transport_options.use_smart_kinetic_solver = params.use_smart_kinetics_solver;
    reactive_transport_options.equilibrium = equilibrium_options;
    reactive_transport_options.smart_equilibrium = smart_equilibrium_options;
    reactive_transport_options.kinetics = kinetic_options; // TODO: think about better structure of the kinetic and equilibrium options
    reactive_transport_options.smart_kinetics = smart_kinetic_options;

    // Step **: Define the reactive transport modeling
    ReactiveTransportSolver rtsolver(reactions, partition);
    rtsolver.setOptions(reactive_transport_options);
    rtsolver.setMesh(mesh);
    rtsolver.setVelocity(params.v);
    rtsolver.setDiffusionCoeff(params.D);
    rtsolver.setBoundaryState(state_bc);
    rtsolver.setTimeStep(params.dt);
    rtsolver.initialize();

    // Step **: Define the quantities that should be output for every cell, every time step
    ChemicalOutput output(rtsolver.output());
    output.add("pH");
    output.add("speciesMolality(H+)");
    output.add("speciesMolality(Al+++)");
    output.add("speciesMolality(K+)");
    output.add("speciesMolality(SiO2(aq))");
    output.add("speciesMolality(CO2(aq))");
    output.add("phaseVolume(K-Feldspar)");
    output.add("phaseVolume(Kaolinite)");
    output.add("speciesAmount(K-Feldspar)");
    output.add("speciesAmount(Kaolinite)");
    output.filename(folder + "/" + "test.txt");

    // Step **: Create RTProfiler to track the timing and results of reactive transport
    ReactiveTransportProfiler profiler;

    // Step **: Set initial time and counter of steps in time
    double t = 0.0;
    int step = 0;

    tic(TRANSPORT);

    // Reactive transport simulations in the cycle
    while (step < params.nsteps)
    {
        // Print some progress
        std::cout << "Step " << step << " of " << params.nsteps << std::endl;

        // Perform one reactive transport time step (with profiling of some parts of the transport simulations)
        rtsolver.stepKinetics(field);

        // Update the profiler after every call to step method
        profiler.update(rtsolver.result());

        // Increment time step and number of time steps
        t += params.dt;

        //getchar();
        step += 1;
    }

    if(params.use_smart_kinetics_solver && params.use_smart_equilibrium_solver) results.time_reactive_transport_smart_kin_smart_eq = toc(TRANSPORT);
    if(!params.use_smart_kinetics_solver && !params.use_smart_equilibrium_solver) results.time_reactive_transport_conv_kin_conv_eq = toc(TRANSPORT);
    if(!params.use_smart_kinetics_solver && params.use_smart_equilibrium_solver) results.time_reactive_transport_conv_kin_smart_eq = toc(TRANSPORT);
    if(params.use_smart_kinetics_solver && !params.use_smart_equilibrium_solver) results.time_reactive_transport_smart_kin_conv_eq = toc(TRANSPORT);

    // Step **: Collect the analytics related to reactive transport performance
    auto analysis = profiler.analysis();
    auto rt_results = profiler.results();

    // Step **: Generate json output file with collected profiling data
    if(params.use_smart_kinetics_solver && params.use_smart_equilibrium_solver) JsonOutput(folder + "/" + "analysis-smart-kin-smart-eq.json") << analysis;
    if(!params.use_smart_kinetics_solver && !params.use_smart_equilibrium_solver) JsonOutput(folder + "/" + "analysis-conventional-kin-conventional-eq.json") << analysis;
    if(!params.use_smart_kinetics_solver && params.use_smart_equilibrium_solver) JsonOutput(folder + "/" + "analysis-conventional-kin-smart-eq.json") << analysis;
    if(params.use_smart_kinetics_solver && !params.use_smart_equilibrium_solver) JsonOutput(folder + "/" + "analysis-smart-kin-conventional-eq.json") << analysis;

    // Step **: Save equilibrium timing to compare the speedup of smart equilibrium solver versus conventional one
    if(params.use_smart_kinetics_solver && params.use_smart_equilibrium_solver) results.smart_kin_smart_eq_total = 0.0;
    if(!params.use_smart_kinetics_solver && !params.use_smart_equilibrium_solver) results.conv_kin_conv_eq_total = 0.0;
    if(!params.use_smart_kinetics_solver && params.use_smart_equilibrium_solver) results.conv_kin_smart_eq_total = 0.0;
    if(params.use_smart_kinetics_solver && !params.use_smart_equilibrium_solver) results.smart_kin_conv_eq_total = 0.0;

    //if(params.use_smart_equilibrium_solver) results.smart_equilibrium_timing = analysis.smart_equilibrium.timing;
    //else results.equilibrium_timing = analysis.equilibrium.timing;
}

/// Make directory for Windows and Linux
auto mkdir(const std::string& folder) -> bool
{
#if defined _WIN32
    // Replace slash by backslash
    std::transform(begin(folder), end(folder), begin(folder),
                   [](char ch) { return ch == '/' ? '\\' : ch; });
    return 0 != CreateDirectory(folder.c_str(), NULL);
#else
    // Create the directory with Read + Write + Execute rights for user, group, and others
    return ::mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
#endif
}

/// Create results file with parameters of the test
auto makeResultsFolder(const Params& params) -> std::string
{
    struct stat status = {0};               // structure to get the file status

    std::ostringstream eqreltol_stream, eqabstol_stream, dt_stream, kinreltol_stream, kinabstol_stream;
    dt_stream << params.dt;
    eqreltol_stream << std::scientific << std::setprecision(1) << params.smart_equilibrium_reltol;
    eqabstol_stream << std::scientific << std::setprecision(1) << params.smart_equilibrium_abstol;
    kinreltol_stream << std::scientific << std::setprecision(1) << params.smart_kinetics_reltol;
    kinabstol_stream << std::scientific << std::setprecision(1) << params.smart_kinetics_abstol;
    std::string test_tag = "-dt-" + dt_stream.str() +
                           "-ncells-" + std::to_string(params.ncells) +
                           "-nsteps-" + std::to_string(params.nsteps) +
                           "-activity-model-" + params.activity_model +
                           "-conv-kin-conv-eq";
    std::string folder = "../plotting-results/rt" + test_tag;
    //std::string folder = "../rt-sa-5000-no-cvode" + test_tag;
    //std::string folder = "../rt-sa-5000" + test_tag;
    if (stat(folder.c_str(), &status) == -1) mkdir(folder);

    std::cout << "\nsolver                         : "
              << (params.use_smart_kinetics_solver ? "smart_kin & " : "conv_kin & ")
              << (params.use_smart_equilibrium_solver ? "smart_eq" : "conv_eq") << std::endl;

    return folder;
}

auto outputConsole(const Params& params) -> void {

    // Log the parameters in the console
    std::cout << "dt      : " << params.dt << std::endl;
    std::cout << "ncells  : " << params.ncells << std::endl;
    std::cout << "nsteps  : " << params.nsteps << std::endl;
    std::cout << "D       : " << params.D << std::endl;
    std::cout << "v       : " << params.v << std::endl;
    std::cout << "CFD     : " << params.v * params.dt / params.dx << std::endl;
    std::cout << "T       : " << params.T << std::endl;
    std::cout << "P       : " << params.P << std::endl;
    std::cout << "activity model : " << params.activity_model << std::endl;
}

