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
#include <Reaktoro/Optimization/Canonicalizer.hpp> // class needed for to fetch primary species

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
    int week = 7 * day;
    // int month = 30 * day;
    // int year = 365 * day;

    // Step 2: Define parameters for the reactive transport simulation
    Params params = {};

    // Define discretization parameters
    params.xl = 0.0; // the x-coordinates of the left boundaries
    params.xr = 1.0; // the x-coordinates of the right boundaries
    params.ncells = 100; // the number of cells in the spacial discretization
    params.nsteps = 2000; // the number of steps in the reactive transport simulation
    params.dx = (params.xr - params.xl) / params.ncells; // the time step (in units of s)
    params.dt = 30 * minute; // the time step (in units of s)

    // Define physical and chemical parameters
    params.D = 1.0e-9;     // the diffusion coefficient (in units of m2/s)
    params.v = 1.0 / week; // the Darcy velocity (in units of m/s)
    params.T = 160.0;                     // the temperature (in units of degC)
    params.P = 100;                      // the pressure (in units of bar)

//    // Run clustering algorithm
//    params.smart_method = "kin-clustering-eq-clustering";
//    params.smart_equilibrium_reltol = 1e-3;
//    params.smart_kinetics_tol = 5e-3;
//    params.smart_kinetics_reltol = 1e-1;
//    params.smart_kinetics_abstol = 1e-1;

    // Run clustering algorithm
    params.smart_method = "kin-clustering-extended-eq-clustering";
    params.smart_equilibrium_reltol = 1e-3;
    params.smart_kinetics_tol = 1e-3;
    params.smart_kinetics_reltol = 1e-1;
    params.smart_kinetics_abstol = 1e-1;

//// Run priority-based queue algorithm
//    params.smart_method = "kin-priority-eq-priority";
//    params.smart_equilibrium_reltol = 1e-3;
//    params.smart_kinetics_tol = 1e-4;
//    params.smart_kinetics_reltol = 1e-1;
//    params.smart_kinetics_abstol = 1e-5;

//// Run priority-based queue algorithm
//    params.smart_method = "kin-priority-potentials-eq-priority-potentials";
//    params.smart_equilibrium_reltol = 1e-3;
//    params.smart_kinetics_tol = 5e-3;
//    params.smart_kinetics_reltol = 1e-1;
//    params.smart_kinetics_abstol = 1e-5;
//    params.output_results = true;
//

//    // Run nn-search algorithm
//    params.smart_method = "kin-nnsearch-eq-nnsearch";
//    params.smart_equilibrium_reltol = 1e-1;
//    params.smart_kinetics_reltol = 1e-1;
//    params.smart_kinetics_abstol = 1e-5;
//    // TODO: smart_kinetics_abstol = 1e-5 is the most optimal tolerance, but there is the shift.

    params.output_results = true;

    params.hessian = GibbsHessian::Exact;
    //params.hessian = GibbsHessian::Approximation;

    params.activity_model = "hkf-full";
    //params.activity_model = "dk-full";
    //params.activity_model = "pitzer-full";

    //params.smart_method = "kin-clustering-eq-clustering";
    //params.smart_method = "kin-priority-eq-priority"; // kin-priority-eq-priority-major-potentials
    //params.smart_method = "kin-nnsearch-eq-nnsearch";

    // Output
    outputConsole(params);

    // RTKineticsResults
    RTKineticsResults results;

//    // **************************************************************************************************************///
//    //  CONVENTIONAL kinetics & SMART equilibrium
//    // **************************************************************************************************************///
//    // Execute reactive transport with different solvers
//    params.use_smart_kinetics_solver = false; params.use_smart_equilibrium_solver = true; runReactiveTransport(params, results);

    // **************************************************************************************************************///
    // SMART kinetics & CONVENTIONAL equilibrium
    // **************************************************************************************************************///
//    params.use_smart_kinetics_solver = true; params.use_smart_equilibrium_solver = false; runReactiveTransport(params, results);

//    // **************************************************************************************************************///
//    // SMART kinetics & SMART equilibrium
//    // **************************************************************************************************************///
//    params.use_smart_kinetics_solver = true;  params.use_smart_equilibrium_solver = true;  runReactiveTransport(params, results);

//    /// **************************************************************************************************************///
//    /// CONVENTIONAL kinetics & CONVENTIONAL equilibrium
//    /// **************************************************************************************************************///
    params.use_smart_kinetics_solver = false; params.use_smart_equilibrium_solver = false; runReactiveTransport(params, results);

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

        std::cout << "time RT smart.kin.& smart.eq. : " << results.time_reactive_transport_smart_kin_smart_eq << std::endl;
        std::cout << "time RT conv.kin.& conv.eq.  : " << results.time_reactive_transport_conv_kin_conv_eq << std::endl;
        std::cout << "speedup                      : " << results.time_reactive_transport_conv_kin_conv_eq
                                                          / results.time_reactive_transport_smart_kin_smart_eq << std::endl;
    }

    std::cout << "total time                          : " << elapsed(start) << std::endl;

    return 0;
}
auto runReactiveTransport(const Params& params, RTKineticsResults& results) -> void {
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
    ChemicalEditor editor;

    // Step **: Add aqueous phase, default chemical model (HKF extended Debye-Hückel model)
    if (params.activity_model == "hkf-full") {
        // HKF full system
        editor.addAqueousPhaseWithElements("H O Na Cl Ca Mg C");
    } else if (params.activity_model == "hkf-selected-species") {
        // HKF selected species
        editor.addAqueousPhase("H2O(l) H+ OH- Na+ Cl- Ca++ Mg++ HCO3- CO2(aq) CO3-- CaCl+ Ca(HCO3)+ MgCl+ Mg(HCO3)+");
    } else if (params.activity_model == "pitzer-full") {
        // Pitzer full system
        editor.addAqueousPhaseWithElements("H O Na Cl Ca Mg C")
                .setChemicalModelPitzerHMW()
                .setActivityModelDrummondCO2();
    } else if (params.activity_model == "pitzer-selected-species") {
        // Pitzer selected species
        editor.addAqueousPhase("H2O(l) H+ OH- Na+ Cl- Ca++ Mg++ HCO3- CO2(aq) CO3-- CaCl+ Ca(HCO3)+ MgCl+ Mg(HCO3)+")
                .setChemicalModelPitzerHMW()
                .setActivityModelDrummondCO2();
    } else if (params.activity_model == "dk-full") {
        // Debye-Huckel full system
        editor.addAqueousPhaseWithElements("H O Na Cl Ca Mg C")
                .setChemicalModelDebyeHuckel()
                .setActivityModelDrummondCO2();
    } else if (params.activity_model == "dk-selected-species") {
        // Debye-Huckel selected species
        editor.addAqueousPhase("H2O(l) H+ OH- Na+ Cl- Ca++ Mg++ HCO3- CO2(aq) CO3-- CaCl+ Ca(HCO3)+ MgCl+ Mg(HCO3)+")
                .setChemicalModelDebyeHuckel()
                .setActivityModelDrummondCO2();
    }
    // Step **: Add mineral phase
    editor.addMineralPhase("Calcite");
    editor.addMineralPhase("Dolomite");
    editor.addMineralPhase("Quartz");
    editor.addMineralPhase("Halite");

    // Step **: Create the ChemicalSystem object using the configured editor
    ChemicalSystem system(editor);

//    std::cout << "system  = " << system << std::endl;
//    getchar();

    std::string eq_str_calcite = "Calcite + H+  =  Ca++ + HCO3-";
    MineralReaction min_reaction_calcite = editor.addMineralReaction("Calcite")
            .setEquation(eq_str_calcite)
            .addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol")
            .addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0")
            .setSpecificSurfaceArea(5000, "cm2/g"); // 0.5 m2/g
    Reaction reaction_calcite = createReaction(min_reaction_calcite, system);
    reaction_calcite.setName("Calcite reaction");
    ReactionRateFunction rate_func_calcite_llnl = [&min_reaction_calcite, &reaction_calcite, &system](const ChemicalProperties& properties) -> ChemicalScalar {

        // The number of chemical species in the system
        const unsigned num_species = system.numSpecies();

        // The mineral reaction rate
        ChemicalScalar res(num_species, 0.0);
        // Auxiliary variable for calculating mineral reaction rate
        ChemicalScalar f(num_species, 1.0);

        // The universal gas constant (in units of kJ/(mol*K))
        const double R = 8.3144621e-3;

        // The temperature and pressure of the system
        const Temperature T = properties.temperature();

        // Calculate the saturation index of the mineral
        const auto lnK = reaction_calcite.lnEquilibriumConstant(properties);
        const auto lnQ = reaction_calcite.lnReactionQuotient(properties);
        const auto lnOmega = lnQ - lnK;

        // Calculate the saturation index
        double Omega = exp(lnOmega).val;
        //std::cout << "Omega (Calcite) = " << Omega << std::endl;

        // The composition of the chemical system
        const auto n = properties.composition().val;
        // The initial composition of the chemical system
        const auto n0 = reaction_calcite.initialAmounts();

        // Amount of elements
        const auto b = system.elementAmounts(n);

        // Calculate TOT("Ca")
        const Index i_ca = system.indexElementWithError("Ca");
        const auto tot_ca = b(i_ca);

        // Calculate TOT("C")
        const Index i_c = system.indexElementWithError("C");
        const auto tot_c = b(i_c);

        // The index of the mineral
        const Index imineral = system.indexSpeciesWithError(min_reaction_calcite.mineral());

        // The current and the initial number of moles of the mineral
        auto nm = n[imineral];
        auto nm0 = n0[imineral];

        VectorConstRef lna = properties.lnActivities().val;
        const Index i_h = system.indexSpeciesWithError("H+");
        double activity_h = std::exp(lna(i_h));

        // The molar mass of the mineral (in units of kg/mol)
        const double molar_mass = system.species(imineral).molarMass();

        // Specific surface area
        const auto ssa = min_reaction_calcite.specificSurfaceArea();

        // Specific surface area correction factor
        auto ssa_corr_factor = 1.0;

        // REM If SI too small, then No dissolution nor precipitation
        // IF (SatIndex > -1e-15) AND ( SatIndex < 1e-15) THEN GOTO 900
        if(lnOmega > -1e-15 && lnOmega < 1e-15)
            return res;
        // IF (SatRatio > 10000) THEN SatRatio = 10000
        if(Omega > 1e4)
            Omega = 1e4;
        // IF (M = 0) AND (SatIndex < PARM(2)) THEN GOTO 900
        // PARM(2): used when M0 = 0, (Nucleation) SatIndex threshold for precipitation
        double sat_index_threshold = 1e6;
        if(nm == 0 && lnOmega < sat_index_threshold)
            return res;
        // IF (M0 > 0) AND (M > 0) THEN t = (M / M0)^PARM(4)
        // PARM(4): exponent for M/M0 for surface area correction (=1 if no data)
        if(nm0 > 0 && nm > 0)
            ssa_corr_factor = pow(nm / nm0, 2 / 3);

        // k25a =	0.501187234
        // Eaa =	14400
        // na =	1
        // IF (k25a > 0) THEN KTa = k25a * EXP((-Eaa/R) * (1/TK - 1/298.15)) * ACT("H+")^na
        const auto kappa_acid = 0.501187234 * exp(- 14.4 / R * (1.0/T - 1.0/298.15)) * std::pow(activity_h, 1.0);
        // k25n =	1.54882E-06
        // Ean =	23500
        // IF (k25n>0) THEN KTn = k25n * EXP((-Ean/R) * (1/TK - 1/298.15))
        const auto kappa_neu = 1.54882e-06 * exp(- 23.5 / R * (1.0/T - 1.0/298.15));
        // k25b =	0
        // Eab =	0
        // nb =	0
        const auto kappa_base = 0.0;

        auto kappa_diss = kappa_neu + kappa_acid + kappa_base;
        auto kappa_pre = kappa_neu + kappa_acid + kappa_base;
        auto rate_constant_ratio = kappa_pre / kappa_diss; // PARM(5): when M0 > 0, rate constant ratio (K_precipitation / kappa_dissolution), (=1 if no data)

        /*
         * Calcite
            -start
            10 k25a =	0.501187234
            20 k25n =	1.54882E-06
            30 k25b =	0
            40 Eaa =	14400
            50 Ean =	23500
            60 Eab =	0
            70 na =	1
            80 nb =	0
            90  SatIndex = SI("Calcite")
            100 SatRatio = SR("Calcite")
            110 R = 8.314472 # (gas constant J/mol/Kelvin)
            120 REM If SI too small, then No dissolution nor precipitation
            130 IF (SatIndex > -1e-15) AND ( SatIndex < 1e-15) THEN GOTO 900
            135 IF (SatRatio > 10000) THEN SatRatio = 10000
            140 IF (M = 0) AND (SatIndex < PARM(2)) THEN GOTO 900
            160 t = 1    # surface area correction factor
            170 IF (M0 > 0) AND (M > 0) THEN t = (M / M0)^PARM(4)
            200 REM acid mechanism
            210   IF (k25a > 0) THEN KTa = k25a * EXP((-Eaa/R) * (1/TK - 1/298.15)) * ACT("H+")^na
            220 REM neutral mechanism
            230   IF (k25n>0) THEN KTn = k25n * EXP((-Ean/R) * (1/TK - 1/298.15))
            240 REM base mechanism
            250   KTb = 0     # Xu(2010)
            300 REM Dissolution
            310   moles = PARM(1) * M0 * t * (KTa+KTn+KTb) * (1 - SatRatio) * TIME
            320 IF (SatIndex > 0) THEN GOTO 400
            330   REM Do not dissolve more than what is available
            340   IF (moles > M) THEN moles = M
            350   GOTO 900
            400 REM Precipitation when M0 > 0
            410 IF (M0 = 0) THEN GOTO 500
            420   moles = moles * PARM(5)
            430   GOTO 600
            500 REM Precipitation when M0 = 0 (secondary mineral)
            520   surfArea = 1.0e-5  # default nucleation surface area
            530   IF (M > 0) THEN surfArea = M * PARM(3)
            540   moles = surfArea * (KTa+KTn+KTb)* PARM(5) * (1 - SatRatio) * TIME
            600 REM Do not precipitate more than the elements in solution
            610   maxMol = TOT("Ca")
            620   IF (maxMol > TOT("C")) THEN maxMol =TOT("C")
            690   IF (maxMol < -moles) THEN moles = -maxMol
            900 SAVE moles
            -end
        -end
         */

        // dissolution kinetics
        // moles = PARM(1) * M0 * t * (KTa+KTn+KTb) * (1 - SatRatio) * TIME
        res = f * ssa * nm0 * ssa_corr_factor * molar_mass * kappa_pre * (1 - Omega);

        if(lnOmega > 0) // lnOmega > 0 // precipitation kinetics
        {
            // REM Precipitation when M0 = 0 (secondary mineral)
            if(nm0 == 0)
            {
                // Default nucleation surface area
                auto nucl_sa = 1.0e-5;
                // IF (M > 0) THEN surfArea = M * PARM(3)
                if(nm) nucl_sa *= nm;

                // moles = surfArea * (KTa+KTn+KTb)* PARM(5) * (1 - SatRatio) * TIME
                res = f * nucl_sa * kappa_pre * rate_constant_ratio * (1 - Omega);
            }
            else // IF (M0 = 0) THEN GOTO 500
                // Set nucleation rate
                // moles = moles * PARM(5) // this taking into account that res was initialized by dissolution kinetic rate
                res *= f * rate_constant_ratio;

//            600 REM Do not precipitate more than the elements in solution
//            610   maxMol = TOT("Ca")
//            620   IF (maxMol > TOT("C")) THEN maxMol =TOT("C")
//            690   IF (maxMol < -moles) THEN moles = -maxMol
//            900 SAVE moles

            // TODO: implement if in the kinetic solver
//            // Implement upper bound in precipitation kinetics
//            auto max_mol = tot_ca;
//            if(max_mol > tot_c) max_mol = tot_c;
//            if(max_mol < -nm.val) return ChemicalScalar(num_species, 0.0);
        }
        else
        {
//            330   REM Do not dissolve more than what is available
//            340   IF (moles > M) THEN moles = M
//            350   GOTO 900
        }
        //std::cout << "res (Calcite) = " << res << std::endl;
        //getchar();
        return res;

    };
    reaction_calcite.setRate(rate_func_calcite_llnl);

    std::string eq_str_dolomite = "Dolomite  + 2*H+  =  Ca++ + Mg++ + 2*HCO3-";
    MineralReaction min_reaction_dolomite = editor.addMineralReaction("Dolomite")
            .setEquation(eq_str_dolomite)
            .addMechanism("logk = -7.53 mol/(m2*s); Ea = 52.2 kJ/mol")
            .addMechanism("logk = -3.19 mol/(m2*s); Ea = 36.1 kJ/mol; a[H+] = 0.5")
            .setSpecificSurfaceArea(5000, "cm2/g"); // 0.5 m2/g
//    Reaction reaction_dolomite = createReaction(min_reaction_dolomite, system);
//    reaction_dolomite.setName("Dolomite reaction");
//    ReactionRateFunction rate_func_dolomite_llnl = [&min_reaction_dolomite, &reaction_dolomite, &system](const ChemicalProperties& properties) -> ChemicalScalar {
//
//        // The number of chemical species in the system
//        const unsigned num_species = system.numSpecies();
//
//        // The mineral reaction rate
//        ChemicalScalar res(num_species, 0.0);
//        // Auxiliary variable for calculating mineral reaction rate
//        ChemicalScalar f(num_species, 1.0);
//
//        // The universal gas constant (in units of kJ/(mol*K))
//        const double R = 8.3144621e-3;
//
//        // The temperature and pressure of the system
//        const Temperature T = properties.temperature();
//
//        // Calculate the saturation index of the mineral
//        const auto lnK = reaction_dolomite.lnEquilibriumConstant(properties);
//        const auto lnQ = reaction_dolomite.lnReactionQuotient(properties);
//        const auto lnOmega = lnQ - lnK;
//
//        // Calculate the saturation index
//        auto Omega = exp(lnOmega).val;
//        std::cout << "Omega (Dolomite) = " << Omega << std::endl;
//
//        // The composition of the chemical system
//        const auto n = properties.composition();
//        // The initial composition of the chemical system
//        const auto n0 = reaction_dolomite.initialAmounts();
//
//        // Amount of elements
//        const auto b = system.elementAmounts(n.val);
//
//        // Calculate TOT("Mg")
//        const Index i_mg = system.indexElementWithError("Mg");
//        const auto tot_mg = b(i_mg);
//
//        // Calculate TOT("C")
//        const Index i_c = system.indexElementWithError("C");
//        const auto tot_c = b(i_c);
//
//        // The index of the mineral
//        const Index imineral = system.indexSpeciesWithError(min_reaction_dolomite.mineral());
//
//        // The current and the initial number of moles of the mineral
//        auto nm = n[imineral].val;
//        auto nm0 = n0[imineral];
//
//        VectorConstRef lna = properties.lnActivities().val;
//        const Index i_h = system.indexSpeciesWithError("H+");
//        double activity_h = std::exp(lna(i_h));
//
//        // The molar mass of the mineral (in units of kg/mol)
//        const double molar_mass = system.species(imineral).molarMass();
//
//        // Specific surface area
//        const auto ssa = min_reaction_dolomite.specificSurfaceArea();
//
//        // Specific surface area correction factor
//        auto ssa_corr_factor = 1.0;
//
//        // REM If SI too small, then No dissolution nor precipitation
//        // IF (SatIndex > -1e-15) AND ( SatIndex < 1e-15) THEN GOTO 900
//        if(lnOmega > -1e-15 && lnOmega < 1e-15)
//            return res;
//        // IF (SatRatio > 10000) THEN SatRatio = 10000
//        if(Omega > 1e4)
//            Omega = 1e4;
//        // IF (M = 0) AND (SatIndex < PARM(2)) THEN GOTO 900
//        // PARM(2): used when M0 = 0, (Nucleation) SatIndex threshold for precipitation
//        double sat_index_threshold = 1e6;
//        if(nm == 0 && lnOmega < sat_index_threshold)
//            return res;
//        // IF (M0 > 0) AND (M > 0) THEN t = (M / M0)^PARM(4)
//        // PARM(4): exponent for M/M0 for surface area correction (=1 if no data)
//        if(nm0 > 0 && nm > 0)
//            ssa_corr_factor = pow(nm / nm0, 2 / 3);
//
////        .addMechanism("logk = -7.53 mol/(m2*s); Ea = 52.2 kJ/mol")
////        .addMechanism("logk = -3.19 mol/(m2*s); Ea = 36.1 kJ/mol; a[H+] = 0.5")
//        // k25n =	2.95121E-08
//        // Ean =	52200
//        // IF (k25n>0) THEN KTn = k25n * EXP((-Ean/R) * (1/TK - 1/298.15))
//        const auto kappa_neu = 2.95121E-08 * exp(- 52.2 / R * (1.0/T - 1.0/298.15));
//        // k25a =	0.000645654
//        // Eaa =	36100
//        // na =	0.5
//        // IF (k25a > 0) THEN KTa = k25a * EXP((-Eaa/R) * (1/TK - 1/298.15)) * ACT("H+")^na
//        const auto kappa_acid = 0.000645654 * exp(- 36.1 / R * (1.0/T - 1.0/298.15)) * std::pow(activity_h, 0.5);
//        // k25b =	0
//        // Eab =	0
//        // nb =	0
//        const auto kappa_base = 0.0;
//
//        auto kappa_diss = kappa_neu + kappa_acid + kappa_base;
//        auto kappa_pre = kappa_neu + kappa_acid + kappa_base;
//        auto rate_constant_ratio = kappa_pre / kappa_diss; // PARM(5): when M0 > 0, rate constant ratio (K_precipitation / kappa_dissolution), (=1 if no data)
//
//        /*
//         * Dolomite-dis
//        -start
//        10 k25a =	0.000645654
//        20 k25n =	2.95121E-08
//        30 k25b =	0
//        40 Eaa =	36100
//        50 Ean =	52200
//        60 Eab =	0
//        70 na =	0.5
//        80 nb =	0
//        90  SatIndex = SI("Dolomite-dis")
//        100 SatRatio = SR("Dolomite-dis")
//        110 R = 8.314472 # (gas constant J/mol/Kelvin)
//        120 REM If SI too small, then No dissolution nor precipitation
//        130 IF (SatIndex > -1e-15) AND ( SatIndex < 1e-15) THEN GOTO 900
//        135 IF (SatRatio > 10000) THEN SatRatio = 10000
//        140 IF (M = 0) AND (SatIndex < PARM(2)) THEN GOTO 900
//        160 t = 1    # surface area correction factor
//        170 IF (M0 > 0) AND (M > 0) THEN t = (M / M0)^PARM(4)
//        200 REM acid mechanism
//        210   IF (k25a > 0) THEN KTa = k25a * EXP((-Eaa/R) * (1/TK - 1/298.15)) * ACT("H+")^na
//        220 REM neutral mechanism
//        230   IF (k25n>0) THEN KTn = k25n * EXP((-Ean/R) * (1/TK - 1/298.15))
//        240 REM base mechanism
//        250   KTb = 0     # Xu(2010)
//        300 REM Dissolution
//        310   moles = PARM(1) * M0 * t * (KTa+KTn+KTb) * (1 - SatRatio) * TIME
//        320 IF (SatIndex > 0) THEN GOTO 400
//        330   REM Do not dissolve more than what is available
//        340   IF (moles > M) THEN moles = M
//        350   GOTO 900
//        400 REM Precipitation when M0 > 0
//        410 IF (M0 = 0) THEN GOTO 500
//        420   moles = moles * PARM(5)
//        430   GOTO 600
//        500 REM Precipitation when M0 = 0 (secondary mineral)
//        520   surfArea = 1.0e-5  # default nucleation surface area
//        530   IF (M > 0) THEN surfArea = M * PARM(3)
//        540   moles = surfArea * (KTa+KTn+KTb)* PARM(5) * (1 - SatRatio) * TIME
//        600 REM Do not precipitate more than the elements in solution
//        610   maxMol = TOT("Ca")
//        620   IF (maxMol > TOT("C")/2) THEN maxMol =TOT("C")/2
//        630   IF (maxMol > TOT("Mg")) THEN maxMol =TOT("Mg")
//        690   IF (maxMol < -moles) THEN moles = -maxMol
//        900 SAVE moles
//        -end
//         */
//
//        // dissolution kinetics
//        // moles = PARM(1) * M0 * t * (KTa+KTn+KTb) * (1 - SatRatio) * TIME
//        res = f * ssa * nm0 * ssa_corr_factor * molar_mass * kappa_pre * (1 - Omega);
//
//        if(lnOmega > 0) // lnOmega > 0 // precipitation kinetics
//        {
//            // REM Precipitation when M0 = 0 (secondary mineral)
//            if(nm0 == 0)
//            {
//                // Default nucleation surface area
//                auto nucl_sa = 1.0e-5;
//                // IF (M > 0) THEN surfArea = M * PARM(3)
//                if(nm) nucl_sa *= nm;
//
//                // moles = surfArea * (KTa+KTn+KTb)* PARM(5) * (1 - SatRatio) * TIME
//                res = f * nucl_sa * kappa_pre * rate_constant_ratio * (1 - Omega);
//            }
//            else // IF (M0 = 0) THEN GOTO 500
//                // Set nucleation rate
//                // moles = moles * PARM(5) // this taking into account that res was initialized by dissolution kinetic rate
//                res *= f * rate_constant_ratio;
//
////            600 REM Do not precipitate more than the elements in solution
////            610   maxMol = TOT("Ca")
////            620   IF (maxMol > TOT("C")/2) THEN maxMol =TOT("C")/2
////            630   IF (maxMol > TOT("Mg")) THEN maxMol =TOT("Mg")
////            690   IF (maxMol < -moles) THEN moles = -maxMol
//            // TODO: implement if in the kinetic solver
////            // Implement upper bound in precipitation kinetics
////            auto max_mol = tot_ca;
////            if(max_mol > tot_c/2) max_mol = tot_c/2;
////            if(max_mol < -nm.val) return ChemicalScalar(num_species, 0.0);
//        }
//        else
//        {
////            330   REM Do not dissolve more than what is available
////            340   IF (moles > M) THEN moles = M
////            350   GOTO 900
//        }
//        std::cout << "lnOmega = " << lnOmega << std::endl;
//        std::cout << "nm0 = " << nm0 << std::endl;
//        std::cout << "res (Dolomite) = " << res << std::endl;
//        getchar();
//        return res;
//    };
//    reaction_dolomite.setRate(rate_func_dolomite_llnl);

    // Step **: Create the ReactionSystem instances
    ReactionSystem reactions(editor);
    //ReactionSystem reactions(system, {reaction_calcite});

    // Step **: Create the ReactionSystem instances
    Partition partition(system);
    //partition.setKineticSpecies(std::vector<std::string>{"Calcite", "Dolomite"});
    partition.setKineticSpecies(std::vector<std::string>{"Calcite"});

    // Step **: Define the initial condition (IC) of the reactive transport modeling problem
    EquilibriumProblem problem_ic(system);
    //problem_ic.setPartition(system);
    problem_ic.setTemperature(params.T, "celsius");
    problem_ic.setPressure(params.P, "bar");
    problem_ic.add("H2O", 1.0, "kg");
    problem_ic.add("O2", 1.0, "umol");
    problem_ic.add("NaCl", 0.7, "mol");
    problem_ic.add("CaCO3", 10, "mol");
    problem_ic.add("CaMg(CO3)2", 10, "mol");
    problem_ic.add("SiO2", 10, "mol");
    problem_ic.add("MgCl2", 1e-10, "mol");

    // Step **: Define the boundary condition (BC)  of the reactive transport modeling problem
    EquilibriumProblem problem_bc(system);
    problem_bc.setTemperature(params.T, "celsius");
    problem_bc.setPressure(params.P, "bar");
    problem_bc.add("H2O", 1.00, "kg");
    problem_bc.add("O2", 1.0, "umol");
    problem_bc.add("NaCl", 0.90, "mol");
    problem_bc.add("MgCl2", 0.05, "mol");
    problem_bc.add("CaCl2", 0.01, "mol");
    problem_bc.add("CO2", 0.75, "mol");

    // Step **: Calculate the equilibrium states for the IC and BC
    ChemicalState state_ic = equilibrate(problem_ic);
    ChemicalState state_bc = equilibrate(problem_bc);

    // Step **: Scale the boundary condition state
    state_bc.scaleVolume(1.0, "m3");

    // Step **: Scale the volumes of the phases in the initial condition
    state_ic.scalePhaseVolume("Aqueous", 0.1, "m3");    // 10% if the 1.0m3
    state_ic.scalePhaseVolume("Quartz", 0.81, "m3");   // 0.81 = 0.90 * 0.9 (0.9 of rock is due to 10% porosity, 0.90 is 90% quartz of the rock)
    state_ic.scalePhaseVolume("Calcite", 0.054, "m3");  // 0.054 = 0.06 * 0.9 (0.9 of rock due to 10% porosity, 0.06 is 6% calcite of the rock)
    state_ic.scalePhaseVolume("Dolomite", 0.036, "m3");  // 0.036 = 0.04 * 0.9 (0.9 of rock due to 10% porosity, 0.04 is 4% dolomite of the rock)

    //reaction_calcite.setInitialAmounts(state_ic.speciesAmounts());
    //reaction_dolomite.setInitialAmounts(state_ic.speciesAmounts());

    //if(params.use_smart_equilibrium_solver) std::cout << "state_ic = \n" << state_ic << std:: endl;
    //getchar();

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
    output.add("speciesMolality(Ca++)");
    output.add("speciesMolality(Mg++)");
    output.add("speciesMolality(HCO3-)");
    output.add("speciesMolality(CO2(aq))");
    output.add("phaseVolume(Calcite)");
    output.add("phaseVolume(Dolomite)");
    output.add("speciesMolality(Calcite)");
    output.add("speciesMolality(Dolomite)");
    output.add("speciesMolality(CO3--)");
    output.add("speciesMolality(CaCl+)");
    output.add("speciesMolality(Ca(HCO3)+)");
    output.add("speciesMolality(MgCl+)");
    output.add("speciesMolality(Mg(HCO3)+)");
    output.add("speciesMolality(OH-)");
    output.add("elementmolality(C)");
    output.add("elementmolality(Ca)");
    output.add("elementmolality(Cl)");
    output.add("elementmolality(H)");
    output.add("elementmolality(Mg)");
    output.add("elementmolality(Na)");
    output.add("elementmolality(O)");
    output.add("elementmolality(Si)");
    output.add("elementmolality(Z)");
    output.add("speciesMolality(MgCO3(aq))");
    output.add("speciesMolality(MgOH+)");
    output.add("speciesAmount(Calcite)");
    output.add("speciesAmount(Dolomite)");
    output.filename(folder + "/" + "test.txt");

    // Step **: Create RTProfiler to track the timing and results of reactive transport
    ReactiveTransportProfiler profiler;

    // Step **: Set initial time and counter of steps in time
    double t = 0.0;
    int step = 0;

    tic(TRANSPORT);

    // Reactive transport simulations in the cycle
    while(step < params.nsteps)
    {
        std::cout << "Step " << step << std::endl;

        // Perform one reactive transport time step (with profiling of some parts of the transport simulations)
        rtsolver.stepKinetics(field);

        // Update the profiler after every call to step method
        profiler.update(rtsolver.result());

        // Increment time step and number of time steps
        t += params.dt;

        step += 1;

    }
    if(params.use_smart_kinetics_solver)
        rtsolver.outputClusterInfo();

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

    //std::string tag = "../plotting-results-03.09.2020/rt-2-react-diss-no-x"; // -kin-clustering-eq-clustering";
    //std::string tag = "../plotting-results-07.09.2020/rt-2-react-diss-neg-eps"; // -kin-clustering-eq-clustering";
    std::string tag = "../plotting-results-24.10.20/rt-2-react-llnl"; // -kin-clustering-eq-clustering";

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

