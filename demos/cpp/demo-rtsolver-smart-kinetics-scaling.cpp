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
    double water_kg = 1.0; // amount of water used in the experiment

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
inline auto errroZeroSurfaceArea(const MineralReaction& reaction) -> void
{
    Exception exception;
    exception.error << "Cannot calculate the molar surface area of the mineral " << reaction.mineral() << ".";
    exception.reason << "The specific surface area of the mineral was not set in reaction " << reaction.equation() << ".";
    RaiseError(exception);
}
auto molarSurfaceArea(const MineralReaction& reaction, const ChemicalSystem& system) -> double
{
    // The temperature and pressure for the calculation of the mineral density
    // Note: These values do not matter much, since the density of the minerals is a constant function
    const double T = 298.15; // in units of kelvin
    const double P = 1.0e5;  // in units of pascal

    // The index of the mineral species
    const Index ispecies = system.indexSpecies(reaction.mineral());

    // The specific surface area of the mineral (in units of m2/kg)
    const double specific_surface_area = reaction.specificSurfaceArea();

    // The molar mass of the mineral (in units of kg/mol)
    const double molar_mass = system.species(ispecies).molarMass();

    // Check if the specific surface area of the mineral was set
    if(specific_surface_area) return specific_surface_area * molar_mass;

    // The standard partial molar volumes at 25 C and 1 bar of all species
    const ThermoVector V = system.properties(T, P).standardPartialMolarVolumes();

    // The molar volume of the mineral species (in units of m3/mol)
    const double molar_volume = V.val[ispecies];

    // The volumetric surface area of the mineral (in units of m2/m3)
    const double volumetric_surface_area = reaction.volumetricSurfaceArea();

    // Check if the volumetric surface area of the mineral was set
    if(volumetric_surface_area) return volumetric_surface_area * molar_volume;

    errroZeroSurfaceArea(reaction);

    return 0.0;
}


int main()
{
    Time start = time();

    // Step 1: Initialise auxiliary time-related constants
    int minute = 60;
    int hour = 60 * minute;

    // Step 2: Define parameters for the reactive transport simulation
    Params params = {};

    // Define discretization parameters
    params.xl = 0.0; // the x-coordinates of the left boundaries
    params.xr = 25.0; // the x-coordinates of the right boundaries
    params.ncells = 243; // the number of cells in the spacial discretization
    params.dx = (params.xr - params.xl) / params.ncells; // the time step (in units of s)
    params.dt = 1*hour; // the time step (in units of s)

    //params.nsteps_cb = 45;  // the number of steps in the reactive transport simulation of the first injection phase
    //params.nsteps_sw = 855;  // the number of steps in the reactive transport simulation of the second injection phase
    params.nsteps = 500;     // the total number of steps in the reactive transport simulation

    // Define physical and chemical parameters
    params.D = 0.0;             // the diffusion coefficient (in units of m2/s)
    params.v = 0.8e-5;          // the Darcy velocity (in units of m/s)
    params.T = 60.0;            // the temperature (in units of degC)
    params.P = 200 * 1.01325;   // the pressure (in units of bar)
    params.water_kg = 1.0;      // amount of water used in the experiment

    // Run clustering algorithm
    params.smart_method = "kin-clustering-eq-clustering";
    params.smart_equilibrium_reltol = 1e-2;
    params.smart_kinetics_tol = 1e-2;
    params.smart_kinetics_reltol = 1e-1;
    params.smart_kinetics_abstol = 1e-4;

//    // Run priority-based queue algorithm
//    params.smart_method = "kin-priority-eq-priority";
//    params.smart_equilibrium_reltol = 5e-3;
//    params.smart_kinetics_tol = 1e-3;
//    params.smart_kinetics_reltol = 1e-1;
//    params.smart_kinetics_abstol = 1e-4;
//    params.output_results = true;

    //
    // // Run nn-search algorithm
    // params.smart_method = "kin-nnsearch-eq-nnsearch";
    // params.smart_equilibrium_reltol = 1e-1;
    // params.smart_kinetics_reltol = 1e-1;
    // params.smart_kinetics_abstol = 1e-4;

    params.output_results = true;

    params.hessian = GibbsHessian::Exact;
    //params.hessian = GibbsHessian::Approximation;

    //params.activity_model = "dk-full";
    params.activity_model = "pitzer-full";
    //params.activity_model = "hkf-full";


    // Output
    outputConsole(params);

    // RTKineticsResults
    RTKineticsResults results;

    // ------------------------------------------------------------------------------------------------------------- //
    // CONVENTIONAL kinetics & SMART equilibrium
    // ------------------------------------------------------------------------------------------------------------- //
    /// Execute reactive transport with different solvers
    // params.use_smart_kinetics_solver = false; params.use_smart_equilibrium_solver = true; runReactiveTransport(params, results);

    /// ------------------------------------------------------------------------------------------------------------- //
    /// SMART kinetics & CONVENTIONAL equilibrium
    /// ------------------------------------------------------------------------------------------------------------- //

    // Execute reactive transport with different solvers
    //params.use_smart_kinetics_solver = true; params.use_smart_equilibrium_solver = false; runReactiveTransport(params, results);

    /// **************************************************************************************************************///
    /// CONVENTIONAL kinetics & CONVENTIONAL equilibrium
    /// **************************************************************************************************************///

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
    //equilibrium_options.hessian = params.hessian;

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
    StringList selected_elements = "H Cl S O Ba Ca Sr Na K Mg C Si";

    if(params.activity_model == "hkf-full"){
        // HKF full system
        editor.addAqueousPhaseWithElements(selected_elements);
    }
    else if(params.activity_model == "pitzer-full"){
        // Pitzer full system
        editor.addAqueousPhaseWithElements(selected_elements)
                .setChemicalModelPitzerHMW()
                .setActivityModelDrummondCO2();
    }
    else if(params.activity_model == "dk-full"){
        // Debye-Huckel full system
        editor.addAqueousPhaseWithElements(selected_elements)
                .setChemicalModelDebyeHuckel(dhModel);
    }
    editor.addMineralPhase("Barite");

    // Step **: Create the ChemicalSystem object using the configured editor
    ChemicalSystem system(editor);
    //std::cout << "system = \n" << system << std:: endl;


    // Barite: BaSO4
    // ########## start precipitation bloc ##########
    //	130 If (m <= 1e-5) then GoTo 170
    //	140 knu = 2.18E-09 * exp((-22000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (10 ^ (0.6 * MU^0.5))
    //	150 kpre = (-1) * knu
    //	160 rate = S * m * Mm * kpre * (SRmin - 1)
    //	170 GoTo 190
    //	#start nucleation
    //	180 rate = -1e-12
    //	190 moles = rate * Time
    //	200 REM Do not precipitate more than the elements in solution
    //	210 maxMol = TOT("Ba")
    //	220 IF (maxMol > TOT("S(6)")) THEN maxMol = TOT("S(6)")
    //	230 IF (maxMol < -moles) THEN moles = -maxMol
    //	########## end precipitation bloc ##########
    // k = std::pow(10.0, logk)
    // log = log10(2.18E-09) = -8.6615
    // Ea = 22
    /*
    MineralReaction min_reaction_barite = editor.addMineralReaction("Barite")
            .setEquation("SO4-- + Ba++ = Barite")
            .addMechanism("logk = -8.6615 mol/(m2*s); Ea = 22 kJ/mol")
            .setSpecificSurfaceArea(0.006, "m2/g");
    Reaction reaction_barite = createReaction(min_reaction_barite, system);
    reaction_barite.setName("Barite precipitation");
    ReactionRateFunction rate_func_barite = [&min_reaction_barite, &system](const ChemicalProperties& properties) -> ChemicalScalar {

        // The number of chemical species in the system
        const unsigned num_species = system.numSpecies();

        // The mineral reaction rate using specified surface area
        ChemicalScalar res(num_species);

        // The universal gas constant (in units of kJ/(mol*K))
        const double R = 8.3144621e-3;

        // The temperature and pressure of the system
        const Temperature T = properties.temperature();

        // Create a Reaction instance
        Reaction reaction(min_reaction_barite.equation(), system);

        for(const MineralMechanism& mechanism : min_reaction_barite.mechanisms()){
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

            // Constant proportional to surface area (or specific surface area)
            auto A = 0.0;
            // Calcultate surface area
            const auto sa = min_reaction_barite.surfaceArea();

            // Calculate the p and q powers of the saturation index Omega
            const auto pOmega = pow(Omega, mechanism.p);
            const auto qOmega = pow(1 - pOmega, mechanism.q);

            if(sa){
                // A is identical to the surface area
                A = sa;
            }
            else{
                const auto ssa = min_reaction_barite.specificSurfaceArea();
                const auto molar_sa = molarSurfaceArea(min_reaction_barite, system);

                // The composition of the chemical system
                const auto n = properties.composition();

                // The index of the mineral
                const Index imineral = system.indexSpeciesWithError(min_reaction_barite.mineral());

                // The number of moles of the mineral
                auto nm = n[imineral];

                // Prevent negative mole numbers here for the solution of the ODEs
                nm.val = std::max(nm.val, 0.0);

                // A is defined by the molar surface area and amount of mineral
                A = molar_sa * nm.val;

//                std::cout << "ssa = " << ssa << std:: endl;
//                std::cout << "molar_sa = " << molar_sa << std:: endl;
//                std::cout << "nm.val = " << nm.val << std:: endl;
//                getchar();
            }

            // Calculate the function f
            f = kappa * qOmega;

            // Calculate the function g
            g = ChemicalScalar(num_species, 1.0);

            // Calculate the resulting mechanism function
            res = f * g;

            // Multiply the mechanism contributions by the surface area of the mineral
            res *= A;

//            std::cout << "lnK = " << lnK << std:: endl;
//            std::cout << "lnQ = " << lnQ << std:: endl;
//            std::cout << "lnOmega = " << lnOmega << std:: endl;
//            std::cout << "Omega = " << Omega << std:: endl;
//            std::cout << "p = " << mechanism.p << std:: endl;
//            std::cout << "q = " << mechanism.q << std:: endl;
//            std::cout << "pOmega = " << pOmega << std:: endl;
//            std::cout << "qOmega = " << qOmega << std:: endl;
//            std::cout << "kappa = " << kappa << std:: endl;
//            std::cout << "sa = " << sa << std::endl;
//            std::cout << "res = " << res << std:: endl;

        }

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
    reaction_barite.setRate(rate_func_barite);
    */

    // Create reaction for the Barite mineral
    std::string eq_str_barite = "Barite = SO4-- + Ba++";
    MineralReaction min_reaction_barite = editor.addMineralReaction("Barite")
            .setEquation(eq_str_barite)
            .addMechanism("logk = -8.6615 mol/(m2*s); Ea = 22 kJ/mol")
            .setSpecificSurfaceArea(0.006, "m2/g");
    Reaction reaction_barite = createReaction(min_reaction_barite, system);
    reaction_barite.setName("Barite reaction");

    // Ionic strength function
    const auto I = ChemicalProperty::ionicStrength(system);

    ReactionRateFunction rate_func_barite_shell = [&min_reaction_barite, &reaction_barite, &system, &I](const ChemicalProperties& properties) -> ChemicalScalar {

        // The number of chemical species in the system
        const unsigned num_species = system.numSpecies();

        // The mineral reaction rate
        ChemicalScalar res(num_species);
        // Auxiliary variable for calculating mineral reaction rate
        ChemicalScalar f(num_species, 1.0);

        // The universal gas constant (in units of kJ/(mol*K))
        const double R = 8.3144621e-3;

        // The temperature and pressure of the system
        const Temperature T = properties.temperature();

        // Create a Reaction instance
        Reaction reaction(min_reaction_barite.equation(), system);


        // Calculate the saturation index of the mineral
        const auto lnK = reaction.lnEquilibriumConstant(properties);
        const auto lnQ = reaction.lnReactionQuotient(properties);
        const auto lnOmega = lnQ - lnK;

        const auto lnK_rxn = reaction_barite.lnEquilibriumConstant(properties);
        const auto lnQ_rxn = reaction_barite.lnReactionQuotient(properties);
        const auto lnOmega_rxn = lnQ_rxn - lnK_rxn;

        // Calculate the saturation index
        const auto Omega = exp(lnOmega);
        const auto Omega_rxn = exp(lnOmega_rxn);

        // The composition of the chemical system
        const auto n = properties.composition();
        // The initial composition of the chemical system
        const auto n0_rxn = reaction_barite.initialAmounts();

//        // Amount of elements
//        const auto b = system.elementAmounts(n.val);
//
//        // Calculate TOT("Ba")
//        const Index i_ba = system.indexElementWithError("Ba");
//        const auto b_ba = b(i_ba);
//        const auto tot_ba = b_ba;
//
//        // Calculate TOT("S(6)") = TOT("SO4")
//        const Index i_s = system.indexElementWithError("S");
//        const auto b_s = b(i_s);
//        const Index i_o = system.indexElementWithError("O");
//        const auto b_o = b(i_o);
//        const auto tot_so4 = b_s + 4 * b_o;

        // The index of the mineral
        const Index imineral = system.indexSpeciesWithError(min_reaction_barite.mineral());

        // The current and the initial number of moles of the mineral
        auto nm = n[imineral];
        auto nm0_rxn = n0_rxn[imineral];

        VectorConstRef lna = properties.lnActivities().val;
        const Index i_h = system.indexSpeciesWithError("H+");
        double activity_h = std::exp(lna(i_h));
        double ionic_strength = I(properties).val;

        // The molar mass of the mineral (in units of kg/mol)
        const double molar_mass = system.species(imineral).molarMass();

        // If (m <= 0) and (SRmin < 1) Then GoTo 240
        // If (SRmin = 1) Then GoTo 240
        if((nm <= 0 && Omega < 1) || (Omega == 1)) // the is no way to precipitate further
            res = 0;
        // S = 0.006 # average BET from 16zhe/did ; suggested value in m2/g
        const auto ssa = 0.006;

        // If (SRmin > 1) Then GoTo 130
        if(Omega > 1) // precipitation kinetics
        {
            /*
             * ########## start precipitation bloc ##########
            130 If (m <= 1e-5) then GoTo 170
            140 knu = 2.18E-09 * exp((-22000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (10 ^ (0.6 * MU^0.5))
            150 kpre = (-1) * knu
            160 rate = S * m * Mm * kpre * (SRmin - 1)
            170 GoTo 190
            #start nucleation
            180 rate = -1e-12
            190 moles = rate * Time
            200 REM Do not precipitate more than the elements in solution
            210 maxMol = TOT("Ba")
            220 IF (maxMol > TOT("S(6)")) THEN maxMol = TOT("S(6)")
            230 IF (maxMol < -moles) THEN moles = -maxMol
            ########## end precipitation bloc ##########
             */

//            // Implement upper bound in precipitation kinetics
//            auto max_mol = tot_ba;
//            if(max_mol > tot_so4) max_mol = tot_so4;
//            if(nm.val >= max_mol) return ChemicalScalar g(num_species, 0.0);

            // If (m <= 1e-5) then GoTo 170
            if(nm > 1e-5)
            {
                // knu = 2.18E-09 * exp((-22000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (10 ^ (0.6 * MU^0.5))
                // kpre = (-1) * knu
                const auto kappa_pre = -2.18e-9 * exp(- 22 / R * (1.0/T - 1.0/298.15))
                                                * pow(10, 0.6 * pow(ionic_strength, 0.5));
                // const auto kappa_pre = -2.18e-9 * exp(- 22 / R * (1.0/T - 1.0/298.15));

                // rate = S * m * Mm * kpre * (SRmin - 1)
                res += f * ssa * nm.val * molar_mass * kappa_pre * (Omega - 1);
            }
            else
                // Set nucleation rate
                res = -1e-12 * f;

        }
        else // dissolution kinetics
        {
            /*
             * ########## start dissolution bloc ##########
            70 k1 = 2.75E-08 * exp((-25000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("H3O+") ^ 0.03) * (10 ^ (0.6 * MU^0.5))
            80 k = k1
            90 rate = S * m * Mm * ((m/m0)^(2/3)) * k * (1 - (SRmin ^ 0.2))
            100 moles = rate * Time
            110 IF (moles > M) THEN moles = M
            120 GoTo 240
            ########## end dissolution bloc ##########
             */

            // k1 = 2.75E-08 * exp((-25000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("H3O+") ^ 0.03) * (10 ^ (0.6 * MU^0.5))
            const auto kappa = 2.75e-8 * exp(- 25 / R * (1.0/T - 1.0/298.15))
                                       * std::pow(activity_h, 0.03)
                                       * pow(10, 0.6 * pow(ionic_strength, 0.5));
            //const auto kappa = 2.75e-8 * exp(- 25 / R * (1.0/T - 1.0/298.15));
            //const auto kappa = 2.75e-8 * exp(- 25 / R * (1.0/T - 1.0/298.15)) * std::pow(activity_h, 0.03);

            // rate = S * m * Mm * ((m/m0)^(2/3)) * k * (1 - (SRmin ^ 0.2))
            res += f * ssa * nm.val * molar_mass * pow(nm.val/nm0_rxn, 2/3) * kappa * (1 - pow(Omega, 0.2));

        }
        return res;

    };
    reaction_barite.setRate(rate_func_barite_shell);

    // Step **: Create the ReactionSystem instances
    //ReactionSystem reactions(editor);
    ReactionSystem reactions(system, {reaction_barite});

    //std::cout << "system = \n" << system << std:: endl;
    //getchar();
    //std::cout << "reactions number = " << reactions.numReactions() << std:: endl;

    // Step **: Create the ReactionSystem instances
    Partition partition(system);
    partition.setKineticSpecies(std::vector<std::string>{"Barite"});

    // ************************************************************************************************************** //
    // Initial condition (IC)
    // ************************************************************************************************************** //

    // Define the initial condition *with formation water)
    EquilibriumInverseProblem problem_ic_fw(system);
    problem_ic_fw.setTemperature(params.T, "celsius");
    problem_ic_fw.setPressure(params.P, "atm");
    problem_ic_fw.add("H2O", params.water_kg, "kg");
    problem_ic_fw.add("SO4", 10 * params.water_kg, "ug");
    problem_ic_fw.add("Ca", 995 * params.water_kg, "mg");
    problem_ic_fw.add("Ba", 995 * params.water_kg, "mg");
    problem_ic_fw.add("Sr", 105 * params.water_kg, "mg");
    problem_ic_fw.add("Na", 27250 * params.water_kg, "mg");
    problem_ic_fw.add("K", 1730 * params.water_kg, "mg");
    problem_ic_fw.add("Mg", 110 * params.water_kg, "mg");
    problem_ic_fw.add("Cl", 45150 * params.water_kg, "mg");
    problem_ic_fw.add("HCO3", 1980 * params.water_kg, "mg");
    problem_ic_fw.pH(7.0, "HCl", "NaOH");

    // Equilibrate the initial condition
    ChemicalState state_ic = equilibrate(problem_ic_fw);

    //state_ic.setSpeciesAmount("Barite", 1.0, "umol");
    state_ic.setSpeciesAmount("Barite", 10, "mcmol");

    // Scale the percentage of the aqueous phase and the whole volume
    state_ic.scalePhaseVolume("Aqueous", 0.1, "m3");    // 10% if the 1.0m3
    state_ic.scaleVolume(1.0, "m3");

    // Set initial value of Barite
    reaction_barite.setInitialAmounts(state_ic.speciesAmounts());

    // Define function to evaluate ph of the chemical system
    auto evaluate_pH = ChemicalProperty::pH(system);

    // Fetch the pH of the initial sate
    ChemicalProperties props = state_ic.properties();
    ChemicalScalar pH_ic = evaluate_pH(props);
    std::cout << "ph(FW) = " << pH_ic.val << std:: endl;
    std::cout << "names b   : "; for (auto elem : system.elements()) std::cout << elem.name() << " "; std::cout << std::endl;
    std::cout << "state_ic b: " << tr(state_ic.elementAmounts()) << std::endl;
    //std::cout << "state_ic b: " << state_ic << std::endl;
    //getchar();

    // ************************************************************************************************************** //
    // Boundary condition (BC)
    // ************************************************************************************************************** //

    // Define the first boundary condition (with completion brine)
    EquilibriumProblem problem_bc_cb(system);
    problem_bc_cb.setTemperature(params.T, "celsius");
    problem_bc_cb.setPressure(params.P, "atm");
    problem_bc_cb.add("H2O", params.water_kg, "kg");
    problem_bc_cb.add("NaCl", 0.7, "mol");

//    problem_ic_fw.add("SO4", 10 * params.water_kg, "ug");
//    problem_ic_fw.add("Ca", 995 * params.water_kg, "mg");
//    problem_ic_fw.add("Ba", 995 * params.water_kg, "mg");
//    problem_ic_fw.add("Sr", 105 * params.water_kg, "mg");
//    problem_ic_fw.add("Na", 27250 * params.water_kg, "mg");
//    problem_ic_fw.add("K", 1730 * params.water_kg, "mg");
//    problem_ic_fw.add("Mg", 110 * params.water_kg, "mg");
//    problem_ic_fw.add("Cl", 45150 * params.water_kg, "mg");
//    problem_ic_fw.add("HCO3", 1980 * params.water_kg, "mg");


    // Equilibrate the initial condition
    ChemicalState state_bc_cb = equilibrate(problem_bc_cb);
    // Scale the percentage of the aqueous phase and the whole volume
    state_bc_cb.scaleVolume(1.0, "m3");

    // Fetch the pH of the initial sate
    props = state_bc_cb.properties();
    ChemicalScalar pH_bc = evaluate_pH(props);
    std::cout << "ph(CB) = " << pH_bc.val << std:: endl;

    std::cout << "state_bc_cb b: " << tr(state_bc_cb.elementAmounts()) << std::endl;
    //std::cout << "state_bc_cb: " << state_bc_cb << std::endl;
    //getchar();

    // Define the first boundary condition (with seawater)
    EquilibriumInverseProblem problem_bc_sw(system);
    problem_bc_sw.setTemperature(params.T, "celsius");
    problem_bc_sw.setPressure(params.P, "atm");
    problem_bc_sw.add("H2O", params.water_kg, "kg");
    problem_bc_sw.add("SO4--", 2710 * params.water_kg, "mg");
    problem_bc_sw.add("Ca++", 411 * params.water_kg, "mg");
    problem_bc_sw.add("Ba++", 0.01 * params.water_kg, "mg");
    problem_bc_sw.add("Sr++", 8 * params.water_kg, "mg");
    problem_bc_sw.add("Na+", 10760 * params.water_kg, "mg");
    problem_bc_sw.add("K+", 399 * params.water_kg, "mg");
    problem_bc_sw.add("Mg++", 1290 * params.water_kg, "mg");
    problem_bc_sw.add("Cl-", 19350 * params.water_kg, "mg");
    problem_bc_sw.add("HCO3-", 142 * params.water_kg, "mg");
    problem_bc_sw.pH(8.1, "HCl", "NaOH");

    // Equilibrate the initial condition
    ChemicalState state_bc_sw = equilibrate(problem_bc_sw);
    // Scale the percentage of the aqueous phase and the whole volume
    state_bc_sw.scaleVolume(1.0, "m3");

    // Fetch the pH of the initial sate
    props = state_bc_sw.properties();
    ChemicalScalar pH_sw = evaluate_pH(props);
    std::cout << "ph(SW) = " << pH_sw.val << std:: endl;

    std::cout << "state_bc_sw b: " << tr(state_bc_sw.elementAmounts()) << std::endl;
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
    reactive_transport_options.kinetics = kinetic_options;
    reactive_transport_options.smart_kinetics = smart_kinetic_options;

    // Step **: Define the reactive transport modeling
    ReactiveTransportSolver rtsolver(reactions, partition);
    rtsolver.setOptions(reactive_transport_options);
    rtsolver.setMesh(mesh);
    rtsolver.setVelocity(params.v);
    rtsolver.setDiffusionCoeff(params.D);
    rtsolver.setBoundaryState(state_bc_sw);
    rtsolver.setTimeStep(params.dt);
    rtsolver.initialize();

    // Step **: Define the quantities that should be output for every cell, every time step
    if(params.output_results)
    {
        ChemicalOutput output(rtsolver.output());
        output.add("pH");
        output.add("speciesAmount(H+)");
        output.add("speciesAmount(Cl-)");
        output.add("speciesAmount(SO4--)");
        output.add("speciesAmount(Ba++)");
        output.add("speciesAmount(Ca++)");
        output.add("speciesAmount(Sr++)");
        output.add("speciesAmount(Na+)");
        output.add("speciesAmount(Barite)");
        output.add("elementAmount(Ba)");
        output.add("elementAmount(C)");
        output.add("elementAmount(Ca)");
        output.add("elementAmount(Cl)");
        output.add("elementAmount(H)");
        output.add("elementAmount(K)");
        output.add("elementAmount(Mg)");
        output.add("elementAmount(Na)");
        output.add("elementAmount(O)");
        output.add("elementAmount(S)");
        output.add("elementAmount(Si)");
        output.add("elementAmount(Sr)");
        output.add("elementAmount(Z)");
        output.add("phaseAmount(Barite)");
        output.add("phaseMass(Barite)");
        output.add("phaseVolume(Barite)");

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
        std::cout << "Step " << step << " of " << params.nsteps << std::endl;

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

    std::string tag = "../plotting-results/rt-scaling-barite-shell-activity-h-ionic-strength"; // -kin-clustering-eq-clustering";
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