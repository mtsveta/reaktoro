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
    Index ncells = 0; // the number of cells in the spacial discretization
    Index nsteps = 0; // the number of steps in the reactive transport simulation
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

    double smart_kinetics_reltol = 0.0;
    double smart_kinetics_abstol = 0.0;
    double smart_kinetics_tol = 0.0;

    bool output_results = true;

    GibbsHessian hessian = GibbsHessian::Exact;

    std::string activity_model = "";

    std::string smart_method = "";

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

    /// Total CPU time (in s) required for equilibrium in the conventional kinetic using equilibrium schemes
    double conv_kin_conv_eq_total_equilibration = 0.0;

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

    // Rate of the smart kinetic estimation w.r.t to the total chemical kinetics calculation
    double smart_kin_conv_eq_acceptance_rate = 0.0;

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


    /// Total CPU time (in s) required for smart equilibrium in the conventional kinetic using smart equilibrium schemes
    double conv_kin_smart_eq_total_smart_equilibration = 0.0;

    // Rate of the smart equilibrium estimation w.r.t to the total chemical kinetics calculation
    double conv_kin_smart_eq_equilibrium_acceptance_rate = 0.0;

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

    // Rate of the smart kinetics estimation w.r.t to the total chemical kinetics calculation
    double smart_kin_smart_eq_acceptance_rate = 0.0;

    // Rate of the smart equilibrium estimation w.r.t to the total chemical kinetics calculation
    double  smart_kin_smart_eq_equilibrium_acceptance_rate = 0.0;

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
    EquilibriumTiming equilibrium_timing = {};

    /// The accumulated timing information of all smart equilibrium calculations.
    SmartEquilibriumTiming smart_equilibrium_timing = {};

    /// Accumulated timing information of all kinetic calculations.
    KineticTiming kinetic_timing = {};

    /// The accumulated timing information of all smart kinetic calculations.
    SmartKineticTiming smart_kinetic_timing = {};

};

/// Forward declaration
auto mkdir(const std::string& folder) -> bool;
auto outputConsole(const Params& params) -> void;
auto makeResultsFolder(const Params& params) -> std::string;
auto runReactiveTransport(const Params& params, RTKineticsResults& results) -> void;
auto outputStatisticsSmartKineticsConventionalEquilibrium(const Params& params, RTKineticsResults& results) -> void;
auto outputStatisticsSmartKineticsSmartEquilibrium(const Params& params, RTKineticsResults& results) -> void;
auto outputStatisticsConventionalKineticsConventionalEquilibrium(const Params& params, RTKineticsResults& results) -> void;
auto outputStatisticsConventionalKineticsSmartEquilibrium(const Params& params, RTKineticsResults& results) -> void;

int main()
{
    Time start = time();

    // Step 1: Initialise auxiliary time-related constants
    // int second = 1;
    int minute = 60;
    int hour = 60 * minute;
    int day = 24 * hour;

    // Step 2: Define parameters for the reactive transport simulation
    Params params = {};

    // Define discretization parameters
    params.xl = 0.0;                                        // the x-coordinates of the left boundaries
    params.xr = 100.0;                                      // the x-coordinates of the right boundaries
    params.ncells = 100;                                    // the number of cells in the spacial discretization
    params.nsteps = 1000;                                   // the number of steps in the reactive transport simulation
    params.dx = (params.xr - params.xl) / params.ncells;    // the time step (in units of s)
    params.dt = 0.1*day;                                    // the time step (in units of s)

    // Define physical and chemical parameters
    params.D = 0.0;     // the diffusion coefficient (in units of m2/s)
    params.v = 3e-5;    // the Darcy velocity (in units of m/s)
    params.T = 25.0;    // the temperature (in units of degC)
    params.P = 1.01325; // the pressure (in units of bar)


//    // Run clustering algorithm
//    params.smart_method = "kin-clustering-eq-clustering";
//    params.smart_equilibrium_reltol = 1e-3;
//    params.smart_kinetics_tol = 1e-3;
//    params.smart_kinetics_reltol = 1e-1;
//    params.smart_kinetics_abstol = 1e-4;

    // Run priority-based queue algorithm
    params.smart_method = "kin-priority-eq-priority";
    params.smart_equilibrium_reltol = 1e-2;
    params.smart_kinetics_tol = 1e-2;
    params.smart_kinetics_reltol = 1e-1;
    params.smart_kinetics_abstol = 1e-5;
    params.output_results = true;

    //
    // // Run nn-search algorithm
    // params.smart_method = "kin-nnsearch-eq-nnsearch";
    // params.smart_equilibrium_reltol = 1e-1;
    // params.smart_kinetics_reltol = 1e-1;
    // params.smart_kinetics_abstol = 1e-4;

    params.output_results = true;

    params.hessian = GibbsHessian::Exact;
    //params.hessian = GibbsHessian::Approximation;

    params.activity_model = "dk-full";
    //params.activity_model = "pitzer-full";

    // Output
    outputConsole(params);

    // RTKineticsResults
    RTKineticsResults results;

    // ------------------------------------------------------------------------------------------------------------- //
    // CONVENTIONAL kinetics & SMART equilibrium
    // ------------------------------------------------------------------------------------------------------------- //
    //params.use_smart_kinetics_solver = false; params.use_smart_equilibrium_solver = true; runReactiveTransport(params, results);

    // ------------------------------------------------------------------------------------------------------------- //
    // SMART kinetics & CONVENTIONAL equilibrium
    // ------------------------------------------------------------------------------------------------------------- //
    //params.use_smart_kinetics_solver = true; params.use_smart_equilibrium_solver = false; runReactiveTransport(params, results);

    // ------------------------------------------------------------------------------------------------------------- //
    /// CONVENTIONAL kinetics & CONVENTIONAL equilibrium
    // ------------------------------------------------------------------------------------------------------------- //
    params.use_smart_kinetics_solver = false; params.use_smart_equilibrium_solver = false; runReactiveTransport(params, results);

    // ------------------------------------------------------------------------------------------------------------- //
    // SMART kinetics & SMART equilibrium
    // ------------------------------------------------------------------------------------------------------------- //
    params.use_smart_kinetics_solver = true;  params.use_smart_equilibrium_solver = true;  runReactiveTransport(params, results);

    // **************************************************************************************************************///
    // SPEED-UP analysis
    // **************************************************************************************************************///
    if(results.smart_kin_conv_eq_total != 0){
        std::cout << "*********************************************************************" << std::endl;
        std::cout << "*********************************************************************" << std::endl;
        std::cout << "speed up of using smart kinetics solver                     : "
                  << results.conv_kin_conv_eq_total / results.smart_kin_conv_eq_total << std::endl;
        std::cout << "speed up ... (with ideal search)                            : "
                  << results.conv_kin_conv_eq_total / results.smart_kin_conv_eq_total_ideal_search << std::endl;
        std::cout << "speed up ... (with ideal search & store)                    : "
                  << results.conv_kin_conv_eq_total / results.smart_kin_conv_eq_total_ideal_search_store << std::endl;
        std::cout << "speed up ... (with ideal search & store & properties eval.) : "
                  << results.conv_kin_conv_eq_total / results.smart_kin_conv_eq_total_ideal_search_store_properties << std::endl;

        std::cout << "time RT conv.kin.& conv.eq.  : " << results.time_reactive_transport_conv_kin_conv_eq << std::endl;
        std::cout << "time RT smart.kin.& conv.eq. : " << results.time_reactive_transport_smart_kin_conv_eq << std::endl;
        std::cout << "speedup                      : " << results.time_reactive_transport_conv_kin_conv_eq
                                                          / results.time_reactive_transport_smart_kin_conv_eq << std::endl;
    }
    if(results.conv_kin_smart_eq_total != 0) {
        std::cout << "*********************************************************************" << std::endl;
        std::cout << "*********************************************************************" << std::endl;
        std::cout << "speed up of using smart equilibrium solver                  : "
                  << results.conv_kin_conv_eq_total / results.conv_kin_smart_eq_total << std::endl;
        std::cout << "speed up ... (with ideal search)                            : "
                  << results.conv_kin_conv_eq_total / results.conv_kin_smart_eq_total_ideal_search << std::endl;
        std::cout << "speed up ... (with ideal search & store)                    : "
                  << results.conv_kin_conv_eq_total / results.conv_kin_smart_eq_total_ideal_search_store << std::endl;
        std::cout << "speed up ... (with ideal search & store & properties eval.) : "
                  << results.conv_kin_conv_eq_total / results.conv_kin_smart_eq_total_ideal_search_store_properties
                  << std::endl;

        std::cout << "speed up in equilibration    : "
                  << results.conv_kin_conv_eq_total_equilibration /
                     results.conv_kin_smart_eq_total_smart_equilibration << std::endl;
        std::cout << "time RT conv.kin.& conv.eq.  : " << results.time_reactive_transport_conv_kin_conv_eq << std::endl;
        std::cout << "time RT conv.kin.& smart.eq. : " << results.time_reactive_transport_conv_kin_smart_eq
                  << std::endl;
        std::cout << "speedup                      : " << results.time_reactive_transport_conv_kin_conv_eq
                                                          / results.time_reactive_transport_conv_kin_smart_eq << std::endl;
    }
    if(results.smart_kin_smart_eq_total != 0){
        std::cout << "*********************************************************************" << std::endl;
        std::cout << "*********************************************************************" << std::endl;
        std::cout << "speed up of using smart kinetic solver (smart equilibrium)  : "
                  << results.conv_kin_conv_eq_total / results.smart_kin_smart_eq_total << std::endl;
        std::cout << "speed up ... (with ideal search)                            : "
                  << results.conv_kin_conv_eq_total / results.smart_kin_smart_eq_total_ideal_search << std::endl;
        std::cout << "speed up ... (with ideal search & store)                    : "
                  << results.conv_kin_conv_eq_total / results.smart_kin_smart_eq_total_ideal_search_store << std::endl;
        std::cout << "speed up ... (with ideal search & store & properties eval.) : "
                  << results.conv_kin_conv_eq_total / results.smart_kin_smart_eq_total_ideal_search_store_properties << std::endl;

        std::cout << "time RT conv.kin.& conv.eq.  : " << results.time_reactive_transport_conv_kin_conv_eq << std::endl;
        std::cout << "time RT smart.kin.& smart.eq. : " << results.time_reactive_transport_smart_kin_smart_eq << std::endl;
        std::cout << "speedup                      : " << results.time_reactive_transport_conv_kin_conv_eq
                                                          / results.time_reactive_transport_smart_kin_smart_eq << std::endl;
    }
    std::cout << "total time                          : " << elapsed(start) << std::endl;

    return 0;
}
auto runReactiveTransport(const Params& params, RTKineticsResults& results) -> void
{
    // Step **: Create the results folder
    auto folder = makeResultsFolder(params);

    // Step **: Define chemical equilibrium solver options
    EquilibriumOptions equilibrium_options;
    equilibrium_options.hessian = params.hessian;

    // Step **: Define smart chemical equilibrium solver options
    SmartEquilibriumOptions smart_equilibrium_options;
    smart_equilibrium_options.reltol = params.smart_equilibrium_reltol;
    smart_equilibrium_options.learning.hessian = params.hessian;
    smart_equilibrium_options.smart_method = params.smart_method;

    // Step **: Define chemical kinetic solver options
    KineticOptions kinetic_options;
    kinetic_options.equilibrium = equilibrium_options;
    kinetic_options.smart_equilibrium = smart_equilibrium_options;
    kinetic_options.use_smart_equilibrium_solver = params.use_smart_equilibrium_solver;

    // Step **: Define smart chemical kinetic solver options
    SmartKineticOptions smart_kinetic_options;
    smart_kinetic_options.reltol = params.smart_kinetics_reltol;
    smart_kinetic_options.abstol = params.smart_kinetics_abstol;
    smart_kinetic_options.tol = params.smart_kinetics_tol;
    smart_kinetic_options.learning = kinetic_options;
    smart_kinetic_options.learning.equilibrium = equilibrium_options;
    smart_kinetic_options.smart_equilibrium = smart_equilibrium_options;
    smart_kinetic_options.use_smart_equilibrium_solver = params.use_smart_equilibrium_solver;
    smart_kinetic_options.smart_method = params.smart_method;

    // Step **: Construct the chemical system with its phases and species (using ChemicalEditor)
    Database database("supcrt07.xml");

    DebyeHuckelParams dhModel{};
    dhModel.setPHREEQC();

    // Step **: Construct the chemical system with its phases and species (using ChemicalEditor)
    ChemicalEditor editor(database);

    StringList selected_species = {"Ca(HCO3)+", "CO3--", "CaCO3(aq)", "Ca++", "CaSO4(aq)", "CaOH+", "Cl-",
                                   "FeCl++", "FeCl2(aq)", "FeCl+", "Fe++", "FeOH+",  "FeOH++", "Fe+++",
                                   "H2(aq)", "HSO4-", "H2S(aq)", "HS-", "H2O(l)",  "H+", "OH-", "HCO3-",
                                   "K+", "KSO4-",
                                   "Mg++", "MgSO4(aq)", "MgCO3(aq)", "MgOH+", "Mg(HCO3)+",
                                   "Na+", "NaSO4-",
                                   "O2(aq)",
                                   "S5--", "S4--", "S3--", "S2--", "SO4--"};
    StringList selected_elements = "C Ca Cl Fe H K Mg Na O S";

    if(params.activity_model == "hkf-full"){
        // HKF full system
        editor.addAqueousPhaseWithElements(selected_elements);
    }
    else if(params.activity_model == "hkf-selected-species"){
        // HKF selected species
        editor.addAqueousPhase(selected_species);
    }
    else if(params.activity_model == "pitzer-full"){
        // Pitzer full system
        editor.addAqueousPhaseWithElements(selected_elements)
                .setChemicalModelPitzerHMW()
                .setActivityModelDrummondCO2();
    }
    else if(params.activity_model == "pitzer-selected-species"){
        // Pitzer selected species
        editor.addAqueousPhase(selected_species)
                .setChemicalModelPitzerHMW()
                .setActivityModelDrummondCO2();
    }
    else if(params.activity_model == "dk-full"){
        // Debye-Huckel full system
        editor.addAqueousPhaseWithElements(selected_elements)
                .setChemicalModelDebyeHuckel()
                .setActivityModelDrummondCO2();
    }
    else if(params.activity_model == "dk-selected-species"){
        // Debye-Huckel selected species
        editor.addAqueousPhase(selected_species)
                .setChemicalModelDebyeHuckel()
                .setActivityModelDrummondCO2();
    }
    editor.addMineralPhase("Siderite");
    editor.addMineralPhase("Pyrite");
    editor.addMineralPhase("Hematite");

    // Step **: Create the ChemicalSystem object using the configured editor
    ChemicalSystem system(editor);
    //std::cout << "system = \n" << system << std:: endl;

    // Hematite: Fe2O3
    // Fe2O3 + 4*H+ = 2*H2O(l) + 2*Fe++ + 0.5*O2(aq)
    // .setEquation("Hematite + 4*H+ = 2*H2O(l) + 2*Fe++ + 0.5*O2(aq)")
    // Fe2O3 + 6*H+ = 3*H2O(l) + 2*Fe+++
    // Fe2O3 + 6H+ = 2Fe+3 + 3H2O
    //Fe2O3 + 6 H+ = 2 Fe+3 + 3 H2O
    /*
    MineralReaction reaction = editor.addMineralReaction("Hematite")
            .setEquation("Hematite + 4*H+ = 2*H2O(l) + 2*Fe++ + 0.5*O2(aq)")
            .addMechanism("logk = -14.60 mol/(m2*s); Ea = 66.2 kJ/mol")
            .addMechanism("logk = -9.39 mol/(m2*s); Ea = 66.2 kJ/mol; a[H+] = 1.0")
            .setSpecificSurfaceArea(1000, "cm2/g");
    */
    /*
    MineralReaction reaction = editor.addMineralReaction("Hematite")
            .setEquation("Hematite + 6*H+ = 3*H2O(l) + 2*Fe+++")
            .addMechanism("logk = -14.60 mol/(m2*s); Ea = 66.2 kJ/mol")
            .addMechanism("logk = -9.39 mol/(m2*s); Ea = 66.2 kJ/mol; a[H+] = 1.0")
            .setSpecificSurfaceArea(1000, "cm2/g");
    */
    //setEquation("Hematite + 4*H+ = 2*H2O(l) + 2*Fe++ + 0.5*O2(aq)")
    // Fe2O3 + 6 H3O+ = 2 Fe+3 + 9 H2O - does not work because H3O+ species is not present in the chemical system
//    MineralReaction reaction = editor.addMineralReaction("Hematite")
//            .setEquation("Hematite + 6*H+ = 3*H2O(l) + 2*Fe+++") // works but hematite precipitates first
//            .addMechanism("logk = -14.60 mol/(m2*s); Ea = 66.2 kJ/mol")
//            .addMechanism("logk = -9.39 mol/(m2*s); Ea = 66.2 kJ/mol; a[H+] = 1.0")
//            .setSpecificSurfaceArea(1000, "cm2/g");

    // Create reaction for the Hematite (Fe2O3) mineral
    //std::string eq_str_hematite = "Hematite + 6*H+ = 3*H2O(l) + 2*Fe+++";
    std::string eq_str_hematite = "Hematite + 6*H+ = 3*H2O(l) + 2*Fe+++";
    MineralReaction min_reaction_hematite = editor.addMineralReaction("Hematite")
            .setEquation(eq_str_hematite)
            .setSpecificSurfaceArea(1.0, "m2/g");
    Reaction reaction_hematite = createReaction(min_reaction_hematite, system);
    reaction_hematite.setName("Hematite reaction");
    ReactionRateFunction rate_func_hematite_shell = [&min_reaction_hematite, &reaction_hematite, &system](const ChemicalProperties& properties) -> ChemicalScalar {

        // The number of chemical species in the system
        const unsigned num_species = system.numSpecies();

        // The mineral reaction rate using specified surface area
        ChemicalScalar res(num_species, 0.0);

        // The universal gas constant (in units of kJ/(mol*K))
        const double R = 8.3144621e-3;

        // The temperature and pressure of the system
        const Temperature T = properties.temperature();

        // Create a Reaction instance
        Reaction reaction(min_reaction_hematite.equation(), system);

        // Auxiliary variables
        ChemicalScalar f(num_species, 1.0);

        // Calculate the saturation index of the mineral
        const auto lnK = reaction_hematite.lnEquilibriumConstant(properties);
        const auto lnQ = reaction_hematite.lnReactionQuotient(properties);
        const auto lnOmega = lnQ - lnK;

        // Calculate the saturation index
        const auto Omega = exp(lnOmega).val;

        // The composition of the chemical system
        const auto n = properties.composition().val;
        const auto n0 = reaction_hematite.initialAmounts();

        // The index of the mineral
        const Index imineral = system.indexSpeciesWithError(min_reaction_hematite.mineral());

        // The number of moles of the mineral
        auto nm = n[imineral];
        auto nm0 = n0[imineral];


        VectorConstRef lna = properties.lnActivities().val;
        const Index i_h = system.indexSpeciesWithError("H+");
        double activity_h = std::exp(lna(i_h));
        const Index i_hs = system.indexSpeciesWithError("HS-");
        double activity_hs = std::exp(lna(i_hs));

        // The molar mass of the mineral (in units of kg/mol)
        const double molar_mass = system.species(imineral).molarMass();

        /*
         * Hematite
        # kinetic data extracted from 04pal/kha 04pou/kro
        # warning dissolution only
        # Confidence level: 4
        -start
        1 SRmin = SR("Hematite")
        10 moles = 0
        20 If (m <= 0) and (SRmin < 1) Then GoTo 250
        30 S = PARM(1) # Default value
        40 Mm = 159.7 # molar mass in g/mol
        50 If (SRmin > 1) Then GoTo 250
        ########## start dissolution bloc ##########
        60 knu = 2.51E-15 * exp((-66200 / 8.314) * ((1 / TK) - (1 / 298.15)))
        70 k1 = 0.000000000407 * exp((-66200 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("H3O+") ^ 1.0)
        80 k2 = 3.50E-9 * exp((-40000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("HS-") ^ 0.5)
        90 k = knu + k1 + k2
        100 rate = S * m * Mm * ((m/m0)^(2/3)) * k * (1 - SRmin) # by default
        110 moles = rate * Time
        120 REM Do not dissolve more than what is available
        130 IF (moles > M) THEN moles = M
        ########## end dissolution bloc ##########
        250 Save moles
        -end
         */
        // If (m <= 0) and (SRmin < 1) Then GoTo 250
        if(nm <= 0 && Omega < 1) // the is no way to precipitate further
            res = ChemicalScalar(num_species, 0.0);
        // S = PARM(1) # Default value
        const auto ssa = min_reaction_hematite.specificSurfaceArea();

        if(Omega < 1) // dissolution kinetics
        {
            /*
             * ########## start dissolution bloc ##########
            60 knu = 2.51E-15 * exp((-66200 / 8.314) * ((1 / TK) - (1 / 298.15)))
            70 k1 = 0.000000000407 * exp((-66200 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("H3O+") ^ 1.0)
            80 k2 = 3.50E-9 * exp((-40000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("HS-") ^ 0.5)
            90 k = knu + k1 + k2
            100 rate = S * m * Mm * ((m/m0)^(2/3)) * k * (1 - SRmin) # by default
            110 moles = rate * Time
            120 REM Do not dissolve more than what is available
            130 IF (moles > M) THEN moles = M
            ########## end dissolution bloc ##########
             */

            // knu = 2.51E-15 * exp((-66200 / 8.314) * ((1 / TK) - (1 / 298.15)))
            const auto kappa_1 = 2.51e-15 * exp(- 66.2 / R * (1.0/T - 1.0/298.15));

            // k1 = 0.000000000407 * exp((-66200 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("H3O+") ^ 1.0)
            const auto kappa_2 = 0.000000000407 * exp( -66.2 / R * (1.0/T - 1.0/298.15)) * std::pow(activity_h, 1.0);

            // Sulfide catalyzer (sulfide promotion)
            // k2 = 3.50E-9 * exp((-40000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("HS-") ^ 0.5)
            const auto kappa_3 = 3.50e-9 * exp(- 40 / R * (1.0/T - 1.0/298.15)) * std::pow(activity_hs, 0.5);

            const auto kappa = kappa_1 + kappa_2 + kappa_3;
            //std::cout << "kappa = " << kappa << std:: endl;

            // Calculate the resulting mechanism function
            // rate = S * m * Mm * ((m/m0)^(2/3)) * k * (1 - SRmin) # by default
            res += f * ssa * nm * molar_mass * pow(nm/nm0, 2/3) * kappa * (1 - Omega);

            // Do not dissolve more than what is available
            // IF (moles > M) THEN moles = M
//            double total_moles = nm.val; // current amount of mols of available minerals
//            if (res > nm.val) res += nm.val;

        }

        return res;

    };
    reaction_hematite.setRate(rate_func_hematite_shell);
    // Pyrite: FeS2
    // 4*FeS2 + 4*H2O(l) =  H+ + SO4-- + 4*Fe++ + 7*HS-
    // FeS2 + 2H+ + 2e- = Fe+2 + 2HS-
    //
    /*
    editor.addMineralReaction("Pyrite")
            .setEquation("Pyrite + H2O(l) = 0.25*H+ + 0.25*SO4-- + Fe++ + 1.75*HS-")
            .addMechanism("logk = -4.55 mol/(m2*s); Ea = 56.9 kJ/mol; a[H+] = 0.500")
            .setSpecificSurfaceArea(10, "cm2/g");
    */
    /*
    editor.addMineralReaction("Pyrite")
            .setEquation("Pyrite + 2*H+ + 2*e- = Fe++ + 2*HS-")
            .addMechanism("logk = -4.55 mol/(m2*s); Ea = 56.9 kJ/mol; a[H+] = 0.500")
            .setSpecificSurfaceArea(10, "cm2/g");
    */
    // Step **: Create the ReactionSystem instances
    ReactionSystem reactions(system, {reaction_hematite});

    //std::cout << "system = \n" << system << std:: endl;
    // getchar();
    //std::cout << "reactions number = " << reactions.numReactions() << std:: endl;

    // Step **: Create the ReactionSystem instances
    Partition partition(system);
    partition.setKineticSpecies(std::vector<std::string>{"Hematite"});
    //partition.setKineticSpecies(std::vector<std::string>{"Pyrite"});
    //partition.setKineticSpecies(std::vector<std::string>{"Hematite", "Pyrite"});
    //partition.setInertSpecies({"Quartz"});

    // Step **: Define the initial condition (IC) of the reactive transport modeling problem
    EquilibriumInverseProblem problem_ic(partition);
    problem_ic.setTemperature(params.T, "celsius");
    problem_ic.setPressure(params.P, "bar");
    problem_ic.add("H2O", 58.0, "kg");
    problem_ic.add("Cl-", 1122.3e-3, "kg");
    problem_ic.add("Na+", 624.08e-3, "kg");
    problem_ic.add("SO4--", 157.18e-3, "kg");
    problem_ic.add("Mg++", 74.820e-3, "kg");
    problem_ic.add("Ca++", 23.838e-3, "kg");
    problem_ic.add("K+", 23.142e-3, "kg");
    problem_ic.add("HCO3-", 8.236e-3, "kg");
    problem_ic.add("O2(aq)", 58e-12, "kg");
    problem_ic.pH(8.951);
    problem_ic.pE(8.676);

    ChemicalState state_ic = equilibrate(problem_ic);

//    std::cout << "state_ic = \n" << state_ic << std:: endl;
//    getchar();

    // Step **: Define the boundary condition (BC)  of the reactive transport modeling problem
    EquilibriumInverseProblem problem_bc(partition);
    problem_bc.setTemperature(params.T, "celsius");
    problem_bc.setPressure(params.P, "bar");
    problem_bc.add("H2O", 58.0, "kg");
    problem_bc.add("Cl-", 1122.3e-3, "kg");
    problem_bc.add("Na+", 624.08e-3, "kg");
    problem_bc.add("SO4--", 157.18e-3, "kg");
    problem_bc.add("Mg++", 74.820e-3, "kg");
    problem_bc.add("Ca++", 23.838e-3, "kg");
    problem_bc.add("K+", 23.142e-3, "kg");
    problem_bc.add("HCO3-", 8.236e-3, "kg");
    problem_bc.add("O2(aq)", 58e-12, "kg");
    problem_bc.add("HS-", 0.0196504, "mol");
    problem_bc.add("H2S(aq)", 0.167794, "mol");
    problem_bc.pH(5.726);
    //problem_bc.pE(8.220);

    // Step **: Calculate the equilibrium states for the IC and BC
    ChemicalState state_bc = equilibrate(problem_bc);
//
//    if (params.use_smart_equilibrium_solver){
//        std::cout << "state_bc = \n" << state_bc << std:: endl;
//        getchar();
//    }

    // Step **: Scale the boundary condition state
    state_bc.scaleVolume(1.0, "m3");

    // Set initial Hematite amount
    state_ic.setSpeciesAmount("Hematite", 0.5, "mol"); // MM(Hematite) = 115.86 g / mol

    // Step **: Scale the volumes of the phases in the initial condition
    state_ic.scalePhaseVolume("Aqueous", 0.1, "m3");    // 10% if the 1.0m3
    state_ic.scaleVolume(1.0, "m3");

    // Set initial value of Hemitite
    reaction_hematite.setInitialAmounts(state_ic.speciesAmounts());

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
    reactive_transport_options.kinetics = kinetic_options;
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
    if(params.output_results)
    {
        ChemicalOutput output(rtsolver.output());
        output.add("pH");
        output.add("speciesAmount(H+)");
        output.add("speciesAmount(HS-)");
        output.add("speciesAmount(S2--)");
        output.add("speciesAmount(CO3--)");
        output.add("speciesAmount(HSO4-)");
        output.add("speciesAmount(H2S(aq))");
        output.add("speciesAmount(Fe++)");
        output.add("speciesAmount(Siderite)");
        output.add("speciesAmount(Pyrite)");
        output.add("speciesAmount(Hematite)");
        output.add("phaseAmount(Siderite)");
        output.add("phaseAmount(Pyrite)");
        output.add("phaseAmount(Hematite)");
        output.add("elementAmount(C)");
        output.add("elementAmount(Ca)");
        output.add("elementAmount(Cl)");
        output.add("elementAmount(Fe)");
        output.add("elementAmount(H)");
        output.add("elementAmount(K)");
        output.add("elementAmount(Mg)");
        output.add("elementAmount(Na)");
        output.add("elementAmount(O)");
        output.add("elementAmount(S)");
        output.add("elementAmount(Z)");
        output.add("speciesAmount(Fe+++)");
        output.filename(folder + "/" + "test.txt");
    }
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
        //if (!(step % 1))
        //std::cout << "Step " << step << " of " << params.nsteps << std::endl;

        // Perform one reactive transport time step (with profiling of some parts of the transport simulations)
        rtsolver.stepKinetics(field);

        // Update the profiler after every call to step method
        profiler.update(rtsolver.result());

        // Increment time step and number of time steps
        t += params.dt;

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
    if(params.use_smart_kinetics_solver && params.use_smart_equilibrium_solver)
        JsonOutput(folder + "/" + "analysis-smart-kin-smart-eq.json") << analysis;
    if(!params.use_smart_kinetics_solver && !params.use_smart_equilibrium_solver)
        JsonOutput(folder + "/" + "analysis-conventional-kin-conventional-eq.json") << analysis;
    if(!params.use_smart_kinetics_solver && params.use_smart_equilibrium_solver)
        JsonOutput(folder + "/" + "analysis-conventional-kin-smart-eq.json") << analysis;
    if(params.use_smart_kinetics_solver && !params.use_smart_equilibrium_solver)
        JsonOutput(folder + "/" + "analysis-smart-kin-conventional-eq.json") << analysis;

    // Step **: Save equilibrium timing to compare the speedup of smart equilibrium solver versus conventional one
    if(params.use_smart_kinetics_solver && params.use_smart_equilibrium_solver) {
        results.smart_kinetic_timing = analysis.smart_kinetics.timing;
        results.smart_equilibrium_timing = analysis.smart_equilibrium.timing;
        results.smart_kin_smart_eq_acceptance_rate = analysis.smart_kinetics.smart_kinetics_estimate_acceptance_rate;
        results.smart_kin_smart_eq_equilibrium_acceptance_rate
                = analysis.smart_equilibrium.smart_equilibrium_estimate_acceptance_rate;
        outputStatisticsSmartKineticsSmartEquilibrium(params, results);
    }
    if(!params.use_smart_kinetics_solver && !params.use_smart_equilibrium_solver) {
        results.kinetic_timing = analysis.kinetics.timing;
        results.equilibrium_timing = analysis.equilibrium.timing;
        outputStatisticsConventionalKineticsConventionalEquilibrium(params, results);
    }
    if(!params.use_smart_kinetics_solver && params.use_smart_equilibrium_solver) {
        results.smart_equilibrium_timing = analysis.smart_equilibrium.timing;
        results.kinetic_timing = analysis.kinetics.timing;
        results.conv_kin_smart_eq_equilibrium_acceptance_rate
                = analysis.smart_equilibrium.smart_equilibrium_estimate_acceptance_rate;
        outputStatisticsConventionalKineticsSmartEquilibrium(params, results);
    }
    if(params.use_smart_kinetics_solver && !params.use_smart_equilibrium_solver) {
        results.smart_kinetic_timing = analysis.smart_kinetics.timing;
        results.equilibrium_timing = analysis.equilibrium.timing;
        results.smart_kin_conv_eq_acceptance_rate = analysis.smart_kinetics.smart_kinetics_estimate_acceptance_rate;
        outputStatisticsSmartKineticsConventionalEquilibrium(params, results);
    }
}

auto outputStatisticsSmartKineticsConventionalEquilibrium(const Params& params, RTKineticsResults& results) -> void
{
    results.smart_kin_conv_eq_total = results.smart_kinetic_timing.solve;
    results.smart_kin_conv_eq_total_ideal_search = results.smart_kinetic_timing.solve
                                                   - results.smart_kinetic_timing.estimate_search;
    results.smart_kin_conv_eq_total_ideal_search_store = results.smart_kinetic_timing.solve
                                                         - results.smart_kinetic_timing.estimate_search
                                                         - results.smart_kinetic_timing.learn_storage;
    results.smart_kin_conv_eq_total_ideal_search_store_properties = results.smart_kinetic_timing.solve
                                                                    - results.smart_kinetic_timing.learn_chemical_properties
                                                                    - results.smart_kinetic_timing.estimate_search
                                                                    - results.smart_kinetic_timing.learn_storage;

    std::cout << "-----------------------------------------------------" << std::endl;
    std::cout << " - solve                : " << results.smart_kinetic_timing.solve << std::endl;
    std::cout << "   - learn              : " << results.smart_kinetic_timing.learn << " (" << results.smart_kinetic_timing.learn / results.smart_kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "     - integrate             : " << results.smart_kinetic_timing.learn_integration << " (" << results.smart_kinetic_timing.learn_integration / results.smart_kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "       - chemical properties : " << results.smart_kinetic_timing.learn_chemical_properties << " (" << results.smart_kinetic_timing.learn_chemical_properties / results.smart_kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "       - integrate_reaction_rates      : " << results.smart_kinetic_timing.learn_reaction_rates << " (" << results.smart_kinetic_timing.learn_reaction_rates / results.smart_kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "       - equilibration       : " << results.smart_kinetic_timing.learn_equilibration << " (" << results.smart_kinetic_timing.learn_equilibration / results.smart_kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "     - store                 : " << results.smart_kinetic_timing.learn_storage << " (" << results.smart_kinetic_timing.learn_storage / results.smart_kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "   - estimate           : " << results.smart_kinetic_timing.estimate << " (" << results.smart_kinetic_timing.estimate / results.smart_kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "     - search                : " << results.smart_kinetic_timing.estimate_search << " (" << results.smart_kinetic_timing.estimate_search / results.smart_kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "     - acceptance            : " << results.smart_kinetic_timing.estimate_error_control << " (" << results.smart_kinetic_timing.estimate_error_control / results.smart_kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "   - equilibrate           : " << results.smart_kinetic_timing.equilibrate << " (" << results.smart_kinetic_timing.equilibrate / results.smart_kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "-----------------------------------------------------" << std::endl;
    std::cout << "-----------------------------------------------------" << std::endl;
    std::cout << " acceptance rate      : " << results.smart_kin_conv_eq_acceptance_rate << " / " << (1 - results.smart_kin_conv_eq_acceptance_rate) * params.ncells *params.nsteps << " learnings out of " << params.ncells *params.nsteps  << std::endl;
    std::cout << "-----------------------------------------------------" << std::endl;
    std::cout << "-----------------------------------------------------" << std::endl;
    std::cout << " - solve - search                 : " << results.smart_kin_conv_eq_total_ideal_search << std::endl;
    std::cout << " - solve - search - store         : " << results.smart_kin_conv_eq_total_ideal_search_store << std::endl;
    std::cout << " - solve - search - store - prop. : " << results.smart_kin_conv_eq_total_ideal_search_store_properties << std::endl;

}

auto outputStatisticsSmartKineticsSmartEquilibrium(const Params& params, RTKineticsResults& results) -> void
{
    results.smart_kin_smart_eq_total = results.smart_kinetic_timing.solve;
    results.smart_kin_smart_eq_total_ideal_search = results.smart_kinetic_timing.solve
                                                    - results.smart_kinetic_timing.estimate_search
                                                    - results.smart_equilibrium_timing.estimate_search;
    results.smart_kin_smart_eq_total_ideal_search_store = results.smart_kinetic_timing.solve
                                                          - results.smart_kinetic_timing.estimate_search
                                                          - results.smart_kinetic_timing.learn_storage
                                                          - results.smart_equilibrium_timing.estimate_search
                                                          - results.smart_equilibrium_timing.learn_storage;
    results.smart_kin_smart_eq_total_ideal_search_store_properties = results.smart_kinetic_timing.solve
                                                                     - results.smart_kinetic_timing.learn_chemical_properties
                                                                     - results.smart_kinetic_timing.estimate_search
                                                                     - results.smart_kinetic_timing.learn_storage
                                                                     - results.smart_equilibrium_timing.learn_chemical_properties
                                                                     - results.smart_equilibrium_timing.estimate_search
                                                                     - results.smart_equilibrium_timing.learn_storage;

    std::cout << "-----------------------------------------------------" << std::endl;
    std::cout << " - solve                : " << results.smart_kinetic_timing.solve << std::endl;
    std::cout << "   - learn              : " << results.smart_kinetic_timing.learn << " (" << results.smart_kinetic_timing.learn / results.smart_kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "     - integrate             : " << results.smart_kinetic_timing.learn_integration << " (" << results.smart_kinetic_timing.learn_integration / results.smart_kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "       - chemical properties      : " << results.smart_kinetic_timing.learn_chemical_properties << " (" << results.smart_kinetic_timing.learn_chemical_properties / results.smart_kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "       - integrate_reaction_rates : " << results.smart_kinetic_timing.learn_reaction_rates << " (" << results.smart_kinetic_timing.learn_reaction_rates / results.smart_kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "       - equilibration            : " << results.smart_kinetic_timing.learn_equilibration << " (" << results.smart_kinetic_timing.learn_equilibration / results.smart_kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "     - store                 : " << results.smart_kinetic_timing.learn_storage << " (" << results.smart_kinetic_timing.learn_storage / results.smart_kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "   - estimate           : " << results.smart_kinetic_timing.estimate << " (" << results.smart_kinetic_timing.estimate / results.smart_kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "     - search                : " << results.smart_kinetic_timing.estimate_search << " (" << results.smart_kinetic_timing.estimate_search / results.smart_kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "     - acceptance            : " << results.smart_kinetic_timing.estimate_error_control << " (" << results.smart_kinetic_timing.estimate_error_control / results.smart_kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "   - equilibrate        : " << results.smart_kinetic_timing.equilibrate << " (" << results.smart_kinetic_timing.equilibrate / results.smart_kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "     - learning              : " << results.smart_equilibrium_timing.learn << " (" << results.smart_equilibrium_timing.learn / results.smart_kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "       - store                    : " << results.smart_equilibrium_timing.learn_storage << " (" << results.smart_equilibrium_timing.learn_storage / results.smart_kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "     - estimation            : " << results.smart_equilibrium_timing.estimate << " (" << results.smart_equilibrium_timing.estimate / results.smart_kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "       - search                   : " << results.smart_equilibrium_timing.estimate_search << " (" << results.smart_equilibrium_timing.estimate_search / results.smart_kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "-----------------------------------------------------" << std::endl;
    std::cout << " smart kinetics acceptance rate      : " << results.smart_kin_smart_eq_acceptance_rate << " ( "
              << results.smart_kin_smart_eq_acceptance_rate * 100 << " % ) / "
              << (1 - results.smart_kin_smart_eq_acceptance_rate) * params.ncells *params.nsteps
              << " learnings out of " << params.ncells *params.nsteps
              << " ( " << (1 - results.smart_kin_smart_eq_acceptance_rate) * 100 << "% )" << std::endl;
    std::cout << " smart equilibrium acceptance rate   : " << results.smart_kin_smart_eq_equilibrium_acceptance_rate << " ( "
              << results.smart_kin_smart_eq_equilibrium_acceptance_rate * 100 << " % ) / "
              << (1 - results.smart_kin_smart_eq_equilibrium_acceptance_rate) * params.ncells *params.nsteps
              << " learnings states out of " << params.ncells *params.nsteps
              << " ( " << (1 - results.smart_kin_smart_eq_equilibrium_acceptance_rate) * 100 << "% )" << std::endl;
    std::cout << "-----------------------------------------------------" << std::endl;
    std::cout << " - solve - search                 : " << results.smart_kin_smart_eq_total_ideal_search << std::endl;
    std::cout << " - solve - search - store         : " << results.smart_kin_smart_eq_total_ideal_search_store << std::endl;
    std::cout << " - solve - search - store - prop. : " << results.smart_kin_smart_eq_total_ideal_search_store_properties << std::endl;
}

auto outputStatisticsConventionalKineticsConventionalEquilibrium(const Params& params, RTKineticsResults& results) -> void
{
    results.conv_kin_conv_eq_total = results.kinetic_timing.solve;
    results.conv_kin_conv_eq_total_ideal_properties = results.kinetic_timing.solve - results.kinetic_timing.integrate_chemical_properties;
    results.conv_kin_conv_eq_total_equilibration = results.kinetic_timing.integrate_equilibration;

    std::cout << "-----------------------------------------------------" << std::endl;
    std::cout << " - solve                   : " << results.kinetic_timing.solve << std::endl;
    std::cout << "   - initialize            : " << results.kinetic_timing.initialize << " (" << results.kinetic_timing.initialize / results.kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "   - integrate             : " << results.kinetic_timing.integrate << " (" << results.kinetic_timing.integrate / results.kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "     - chemical properties      : " << results.kinetic_timing.integrate_chemical_properties << " (" << results.kinetic_timing.integrate_chemical_properties / results.kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "     - integrate_reaction_rates : " << results.kinetic_timing.integrate_reaction_rates << " (" << results.kinetic_timing.integrate_reaction_rates / results.kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "     - equilibration            : " << results.kinetic_timing.integrate_equilibration << " (" << results.kinetic_timing.integrate_equilibration / results.kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "   - equilibrate           : " << results.kinetic_timing.equilibrate << " (" << results.kinetic_timing.equilibrate / results.kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "-----------------------------------------------------" << std::endl;
}

auto outputStatisticsConventionalKineticsSmartEquilibrium(const Params& params, RTKineticsResults& results) -> void
{
    results.conv_kin_smart_eq_total = results.kinetic_timing.solve;
    results.conv_kin_smart_eq_total_ideal_search = results.kinetic_timing.solve - results.smart_equilibrium_timing.estimate_search;
    results.conv_kin_smart_eq_total_ideal_search_store = results.kinetic_timing.solve
                                                         - results.smart_equilibrium_timing.estimate_search
                                                         - results.smart_equilibrium_timing.learn_storage;
    results.conv_kin_smart_eq_total_ideal_search_store_properties = results.kinetic_timing.solve
                                                                    - results.smart_equilibrium_timing.estimate_search
                                                                    - results.smart_equilibrium_timing.learn_storage
                                                                    - results.kinetic_timing.integrate_chemical_properties;
    results.conv_kin_smart_eq_total_smart_equilibration = results.kinetic_timing.equilibrate;

    std::cout << "-----------------------------------------------------" << std::endl;
    std::cout << " - solve                   : " << results.kinetic_timing.solve << std::endl;
    std::cout << "   - initialize            : " << results.kinetic_timing.initialize << " (" << results.kinetic_timing.initialize / results.kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "   - integrate             : " << results.kinetic_timing.integrate << " (" << results.kinetic_timing.integrate / results.kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "     - chemical properties : " << results.kinetic_timing.integrate_chemical_properties << " (" << results.kinetic_timing.integrate_chemical_properties / results.kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "     - integrate_reaction_rates      : " << results.kinetic_timing.integrate_reaction_rates << " (" << results.kinetic_timing.integrate_reaction_rates / results.kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "     - equilibration       : " << results.kinetic_timing.equilibrate << " (" << results.kinetic_timing.equilibrate / results.kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "   - equilibrate             : " << results.kinetic_timing.equilibrate << " (" << results.kinetic_timing.equilibrate / results.kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "     - learning          : " << results.smart_equilibrium_timing.learn << " (" << results.smart_equilibrium_timing.learn / results.kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "       - store           : " << results.smart_equilibrium_timing.learn_storage << " (" << results.smart_equilibrium_timing.learn_storage / results.kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "     - estimation        : " << results.smart_equilibrium_timing.estimate << " (" << results.smart_equilibrium_timing.estimate / results.kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "       - search          : " << results.smart_equilibrium_timing.estimate_search << " (" << results.smart_equilibrium_timing.estimate_search / results.kinetic_timing.solve * 100 << " %)" << std::endl;
    std::cout << "-----------------------------------------------------" << std::endl;
    std::cout << " smart equilibrium acceptance rate   : " << results.conv_kin_smart_eq_equilibrium_acceptance_rate << " ( "
              << results.conv_kin_smart_eq_equilibrium_acceptance_rate * 100 << " % ) / "
              << (1 - results.conv_kin_smart_eq_equilibrium_acceptance_rate) * params.ncells *params.nsteps
              << " learnings states out of " << params.ncells *params.nsteps
              << " ( " << (1 - results.conv_kin_smart_eq_equilibrium_acceptance_rate) * 100 << "% )" << std::endl;
    std::cout << "-----------------------------------------------------" << std::endl;
    std::cout << " - solve - search              : " << results.kinetic_timing.solve
                                                        - results.smart_equilibrium_timing.estimate_search << std::endl;
    std::cout << " - solve - search - store      : " << results.kinetic_timing.solve
                                                        - results.smart_equilibrium_timing.estimate_search
                                                        - results.smart_equilibrium_timing.learn_storage << std::endl;
    std::cout << " - solve - search - store - properties : " << results.kinetic_timing.solve
                                                                - results.kinetic_timing.integrate_chemical_properties
                                                                - results.smart_equilibrium_timing.learn_storage
                                                                - results.smart_equilibrium_timing.estimate_search << std::endl;
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

    std::ostringstream eqreltol_stream,
            eqabstol_stream,
            eqcutoff_stream,
            eqtol_stream,
            dt_stream,
            kinreltol_stream,
            kinabstol_stream,
            kintol_stream;
    dt_stream << params.dt;
    eqreltol_stream << std::scientific << std::setprecision(1) << params.smart_equilibrium_reltol;
    eqcutoff_stream << std::scientific << std::setprecision(1) << params.smart_equilibrium_cutoff;
    kinreltol_stream << std::scientific << std::setprecision(1) << params.smart_kinetics_reltol;
    kinabstol_stream << std::scientific << std::setprecision(1) << params.smart_kinetics_abstol;
    kintol_stream << std::scientific << std::setprecision(1) << params.smart_kinetics_tol;
    std::string test_tag = "-dt-" + dt_stream.str() +
                           "-ncells-" + std::to_string(params.ncells) +
                           "-nsteps-" + std::to_string(params.nsteps) +
                           "-" + params.activity_model +
                           (params.use_smart_kinetics_solver ? "-smart-kin" : "-conv-kin") +
                           (params.use_smart_equilibrium_solver ? "-smart-eq"  : "-conv-eq");      // name of the folder with results
    std::string smart_test_tag = "-" + params.smart_method +
                                 "-dt-" + dt_stream.str() +
                                 "-ncells-" + std::to_string(params.ncells) +
                                 "-nsteps-" + std::to_string(params.nsteps) +
                                 "-eqrel-" + eqreltol_stream.str() +
                                 "-kinrel-" + kinreltol_stream.str() +
                                 "-kinabs-" + kinabstol_stream.str() +
                                 "-kintol-" + kintol_stream.str() +
                                 "-" + params.activity_model +
                                 (params.use_smart_kinetics_solver ? "-smart-kin" : "-conv-kin") +
                                 (params.use_smart_equilibrium_solver ? "-smart-eq"  : "-conv-eq");      // name of the folder with results

    std::string tag = "../plotting-results/rt-scavenging-hematite-supcrt07"; // -kin-clustering-eq-clustering";
    //std::string tag = "../plotting-results/rt-scavenging-with-hematite-1000-pyrite"; // -kin-clustering-eq-clustering";
    //std::string tag = "../plotting-results/rt-scavenging-no-kinetics"; // -kin-clustering-eq-clustering";
    //std::string tag = "../plotting-results/rt-kin-priority-eq-clustering";
    //std::string tag = "../plotting-results/rt-kin-priority-eq-priority-primary-potentials";
    //std::string tag = "../plotting-results/rt-new-ode-solve-kin-priority-eq-priority-no-state-copying-no-A";
    //std::string tag = "../plotting-results/rt-new-ode-solve-kin-priority-eq-priority-success-false";
    //std::string tag = "../plotting-results/rt-new-ode-solve-kin-priority-eq-nn-residual";
    //std::string tag = "../plotting-results/rt-new-ode-solve-kin-nn-residual-eq-nn-residual";
    //std::string tag = "../plotting-results/rt-kin-nn-search-residual-eq-nn-search-lna";
    //std::string tag = "../plotting-results/rt-kin-nn-search-lna-eq-nn-search-lna";
    //std::string tag = "../plotting-results/rt-kin-nn-search-residual-eq-nn-search-residual";
    //std::string tag = "../plotting-results/rt-kin-nn-search-potentials";
    //std::string tag = "../plotting-results/rt-kin-nn-search-residual";
    //std::string tag = "../plotting-results/rt-kin-nn-search-lna";
    //std::string tag = "../plotting-results/rt-kin-nn-search-lna-f-J-update-of-rates-sens-props";
    //std::string tag = "../plotting-results/rt-kin-nn-search-lna-update-of-rates-sens-props";
    //std::string tag = "../plotting-results/rt-old-kinetics";
    std::string folder =
            (params.use_smart_kinetics_solver || params.use_smart_equilibrium_solver) ?
            tag + smart_test_tag :
            tag + test_tag;
    if (stat(folder.c_str(), &status) == -1) mkdir(folder);

    std::cout << "*********************************************************************" << std::endl;
    std::cout << "*********************************************************************" << std::endl;
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
    std::cout << "equilibrium reltol   : " << params.smart_equilibrium_reltol << std::endl;
    std::cout << "kinetics reltol      : " << params.smart_kinetics_reltol << std::endl;
    std::cout << "kinetics abstol      : " << params.smart_kinetics_abstol << std::endl;
    std::cout << "kinetics tol         : " << params.smart_kinetics_tol << std::endl;
    std::cout << "activity model       : " << params.activity_model << std::endl;
    std::cout << "smart method         : " << params.smart_method << std::endl;

}


/*
 * system =
====================================================================================================
Aqueous                  Siderite                 Pyrite                   Hematite
----------------------------------------------------------------------------------------------------
CO(aq)                   Siderite                 Pyrite                   Hematite
CO2(aq)
CO3--
Ca(HCO3)+
Ca++
CaCO3(aq)
CaCl+
CaCl2(aq)
CaOH+
CaSO4(aq)
Cl-
ClO-
ClO2-
ClO3-
ClO4-
Fe++
Fe+++
FeCl+
FeCl++
FeCl2(aq)
FeO(aq)
FeO+
FeO2-
FeOH+
FeOH++
H+
H2(aq)
H2O(l)
H2O2(aq)
H2S(aq)
H2S2O3(aq)
H2S2O4(aq)
HCO3-
HCl(aq)
HClO(aq)
HClO2(aq)
HFeO2(aq)
HFeO2-
HO2-
HS-
HS2O3-
HS2O4-
HSO3-
HSO4-
HSO5-
K+
KCl(aq)
KHSO4(aq)
KOH(aq)
KSO4-
Mg(HCO3)+
Mg++
MgCO3(aq)
MgCl+
MgOH+
MgSO4(aq)
Na+
NaCl(aq)
NaOH(aq)
NaSO4-
O2(aq)
OH-
S2--
S2O3--
S2O4--
S2O5--
S2O6--
S2O8--
S3--
S3O6--
S4--
S4O6--
S5--
S5O6--
SO2(aq)
SO3--
SO4--
====================================================================================================
Index                    Species                  Element                  Phase
----------------------------------------------------------------------------------------------------
0                        CO(aq)                   C                        Aqueous
1                        CO2(aq)                  Ca                       Siderite
2                        CO3--                    Cl                       Pyrite
3                        Ca(HCO3)+                Fe                       Hematite
4                        Ca++                     H
5                        CaCO3(aq)                K
6                        CaCl+                    Mg
7                        CaCl2(aq)                Na
8                        CaOH+                    O
9                        CaSO4(aq)                S
10                       Cl-                      Z
11                       ClO-
12                       ClO2-
13                       ClO3-
14                       ClO4-
15                       Fe++
16                       Fe+++
17                       FeCl+
18                       FeCl++
19                       FeCl2(aq)
20                       FeO(aq)
21                       FeO+
22                       FeO2-
23                       FeOH+
24                       FeOH++
25                       H+
26                       H2(aq)
27                       H2O(l)
28                       H2O2(aq)
29                       H2S(aq)
30                       H2S2O3(aq)
31                       H2S2O4(aq)
32                       HCO3-
33                       HCl(aq)
34                       HClO(aq)
35                       HClO2(aq)
36                       HFeO2(aq)
37                       HFeO2-
38                       HO2-
39                       HS-
40                       HS2O3-
41                       HS2O4-
42                       HSO3-
43                       HSO4-
44                       HSO5-
45                       K+
46                       KCl(aq)
47                       KHSO4(aq)
48                       KOH(aq)
49                       KSO4-
50                       Mg(HCO3)+
51                       Mg++
52                       MgCO3(aq)
53                       MgCl+
54                       MgOH+
55                       MgSO4(aq)
56                       Na+
57                       NaCl(aq)
58                       NaOH(aq)
59                       NaSO4-
60                       O2(aq)
61                       OH-
62                       S2--
63                       S2O3--
64                       S2O4--
65                       S2O5--
66                       S2O6--
67                       S2O8--
68                       S3--
69                       S3O6--
70                       S4--
71                       S4O6--
72                       S5--
73                       S5O6--
74                       SO2(aq)
75                       SO3--
76                       SO4--
77                       Siderite
78                       Pyrite
79                       Hematite
====================================================================================================
 */


/*
Hematite
Fe2O3 + 6 H3O+ = 2 Fe+3 + 9 H2O
log_k           -4.008
delta_h -30.845 kcal

Hematite
Fe2O3 + 6 H+ = 2 Fe+3 + 3 H2O
-log_k	-4.008
-delta_h -30.845 kcal
-Vm 30.39

Hematite            108
Fe2O3 + 6H+ = 2Fe+3 + 3H2O
log_k           -4.008
delta_h -30.845 kcal


Hematite
	Fe2O3 + 6 H+ = 2 Fe+3 + 3 H2O
	log_k		0.1086
	-delta_H	-129.415	kJ/mol
#	deltafH		-197.72		kcal/mol
	-analytic	-2.2015e2 -6.0290e-2 1.1812e4 8.0253e1 1.8438e2
#	Range		0-350
	-Vm		30.274
#	Extrapol	supcrt92
#	Ref		HDN+78

Hematite
        Fe2O3 + 6 H3O+ = 2 Fe+3 + 9 H2O
        log_k           -4.008
        delta_h -30.845 kcal

Hematite
        Fe2O3 +6.0000 H+  =  + 2.0000 Fe+++ + 3.0000 H2O
        log_k           0.1086
	-delta_H	-129.415	kJ/mol	# Calculated enthalpy of reaction	Hematite
#	Enthalpy of formation:	-197.72 kcal/mol
        -analytic -2.2015e+002 -6.0290e-002 1.1812e+004 8.0253e+001 1.8438e+002
#       -Range:  0-300

Hematite
        Fe2O3 + 6H+ = 2Fe+3 + 3H2O
        log_k   -4.008
        delta_h -30.845 kcal

Hematite
Fe2O3 + 6H+ = 2Fe+3 + 3H2O
log_k	-1.418
delta_h	-128.987	kJ


Hematite
Fe2O3 + 6 H+ = 2 Fe+3 + 3 H2O
-log_k	-4.008
-delta_h -30.845 kcal
-Vm 30.39


Hematite
Fe2O3     = 2.000Fe+3     - 6.000H+     + 3.000H2O
log_k    -1.020     #05GRI
delta_h -123.679    #kJ/mol
# Enthalpy of formation:           -831.811        #kJ/mol
-analytic -2.26876E+1 0E+0 6.46019E+3 0E+0 0E+0


Hematite            108
Fe2O3 + 6H+ = 2Fe+3 + 3H2O
log_k           -4.008
delta_h -30.845 kcal

 */