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
    double xl = 0; // the x-coordinates of the left boundaries
    double xr = 0; // the x-coordinates of the right boundaries
    double dx = 0; // the space step (in units of m)
    double dt = 0; // the time step (in units of s)

    int nsteps_cb = 0;  // the number of steps in the reactive transport simulation of the first injection phase
    int nsteps_sw = 0;  // the number of steps in the reactive transport simulation of the second injection phase
    int nsteps = 0;     // the total number of steps in the reactive transport simulation

    // Physical params
    double D = 0; // the diffusion coefficient (in units of m2/s)
    double v = 0; // the Darcy velocity (in units of m/s)
    double T = 0; // the temperature (in units of degC)
    double P = 0; // the pressure (in units of bar)
    double water_kg = 1.0;  // amount of water used in the experiment

    // Solver params
    bool use_smart_equilibrium_solver = false;
    double smart_equilibrium_reltol = 0;
    double amount_fraction_cutoff = 0;
    double mole_fraction_cutoff = 0;

    GibbsHessian hessian = GibbsHessian::Exact;

    std::string activity_model = "";

    std::string smart_method = "";
};


struct Results
{
    /// Total CPU time (in s) required by smart equilibrium scheme
    double smart_total;

    /// Total CPU time (in s) excluding the costs for the search of the closest reference states.
    double smart_total_ideal_search;

    /// Total CPU time (in s) required by smart equilibrium scheme
    /// excluding the costs for the search and storage of the closest reference states.
    double smart_total_ideal_search_store;

    /// Total CPU time (in s) required by conventional equilibrium scheme
    double conventional_total;

    /// The total time taken to perform all time steps using conventional equilibrium algorithm
    double time_reactive_transport_conventional;

    /// The total time taken to perform all time steps using smart equilibrium algorithm
    double time_reactive_transport_smart;

    /// The accumulated timing information of all equilibrium calculations.
    EquilibriumTiming equilibrium_timing;

    /// The accumulated timing information of all smart equilibrium calculations.
    SmartEquilibriumTiming smart_equilibrium_timing;

    // Rate of the smart equlibrium estimation w.r.t to the total chemical equilibrium calculation
    double smart_equilibrium_acceptance_rate = 0.0;
};

/// Forward declaration
auto mkdir(const std::string& folder) -> bool;
auto outputConsole(const Params& params) -> void;
auto makeResultsFolder(const Params& params) -> std::string;
auto runReactiveTransport(const Params& params, Results& results) -> void;

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
    params.nsteps = 500;                                   // the number of steps in the reactive transport simulation
    params.dx = (params.xr - params.xl) / params.ncells;    // the time step (in units of s)
    params.dt = 0.1*day;                                    // the time step (in units of s)

    // Define physical and chemical parameters
    params.D = 0.0;     // the diffusion coefficient (in units of m2/s)
    params.v = 3e-5;    // the Darcy velocity (in units of m/s)
    params.T = 25.0;    // the temperature (in units of degC)
    params.P = 1.01325; // the pressure (in units of bar)

    params.activity_model = "dk-full";
    //params.activity_model = "pitzer-full";
    //params.activity_model = "hkf-full";

    // Run clustering algorithm
    //params.smart_method = "eq-clustering";
    //params.smart_equilibrium_reltol = 5e-3;

    // Run clustering algorithm
    params.smart_method = "eq-priority";
    params.smart_equilibrium_reltol = 2e-3;
    params.amount_fraction_cutoff = 1e-14;
    params.mole_fraction_cutoff = 1e-14;

    // Output
    outputConsole(params);

    // Results
    Results results;

    // Execute reactive transport with different solvers
    params.use_smart_equilibrium_solver = true; runReactiveTransport(params, results);
    //params.use_smart_equilibrium_solver = false; runReactiveTransport(params, results);

    results.conventional_total = results.equilibrium_timing.solve;
    results.smart_total = results.smart_equilibrium_timing.solve;
    results.smart_total_ideal_search = results.smart_equilibrium_timing.solve
                                       - results.smart_equilibrium_timing.estimate_search
                                       - results.smart_equilibrium_timing.estimate_database_priority_update;
    results.smart_total_ideal_search_store = results.smart_equilibrium_timing.solve
                                             - results.smart_equilibrium_timing.estimate_search
                                             - results.smart_equilibrium_timing.estimate_database_priority_update
                                             - results.smart_equilibrium_timing.learn_storage;

    // Output speed-us
    std::cout << "speed up                            : "
              << results.conventional_total / results.smart_total << std::endl;
    std::cout << "speed up (with ideal search)        : "
              << results.conventional_total / results.smart_total_ideal_search << std::endl;
    std::cout << "speed up (with ideal search & store): "
              << results.conventional_total / results.smart_total_ideal_search_store << std::endl << std::endl;

    std::cout << "time_reactive_transport_conventional: " << results.time_reactive_transport_conventional << std::endl;
    std::cout << "time_reactive_transport_smart       : " << results.time_reactive_transport_smart << std::endl;
    std::cout << "reactive_transport_speedup          : " << results.time_reactive_transport_conventional / results.time_reactive_transport_smart << std::endl;

    std::cout << "total time                          : " << elapsed(start) << std::endl;

    return 0;
}
auto runReactiveTransport(const Params& params, Results& results) -> void
{
    // Step **: Create the results folder
    auto folder = makeResultsFolder(params);

    // Step **: Define chemical equilibrium solver options
    EquilibriumOptions equilibrium_options;
    equilibrium_options.hessian = params.hessian;

    // Step **: Define smart chemical equilibrium solver options
    SmartEquilibriumOptions smart_equilibrium_options;
    smart_equilibrium_options.reltol = params.smart_equilibrium_reltol;
    smart_equilibrium_options.amount_fraction_cutoff = params.amount_fraction_cutoff;
    smart_equilibrium_options.mole_fraction_cutoff = params.mole_fraction_cutoff;
    smart_equilibrium_options.smart_method = params.smart_method;


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
    StringList selected_elements = "Al C Ca Cl Fe H K Mg Na O S Si";

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
    editor.addMineralPhase("Pyrrhotite");
    editor.addMineralPhase("Pyrite");
    editor.addMineralPhase("Calcite");
    editor.addMineralPhase("Daphnite,14A"); // Chlorite
    editor.addMineralPhase("Siderite");
    editor.addMineralPhase("Kaolinite");
    editor.addMineralPhase("Quartz");


    // Step **: Create the ChemicalSystem object using the configured editor
    ChemicalSystem system(editor);
    std::cout << "system = \n" << system << std:: endl;
    getchar();

    const auto I = ChemicalProperty::ionicStrength(system);

    //Chamosite(Daphnite)
    //Fe5Al(AlSi3)O10(OH)8 + 16.000H+ = 2.000Al+3 + 5.000Fe+2 + 3.000H4SiO4 + 6.000H2O
    //     log_k    47.579
    //     delta_h -504.518    #kJ/mol        #01vid/par
    //     -analytic -2.6210061E+3  -4.2497094E-1  1.5576281E+5  9.4858884E+2  -6.610337E+6
    //     #References = LogK/DGf: Internal calculation; DHf/DHr: 01vid/par; S°: 01vid/par; Cp: 05vid/par; V°: 05vid/par;

    std::string eq_str_daphnite = "Daphnite,14A + 16*H+ = 2*Al+++ + 5*Fe++ + 3*SiO2(aq) + 9*H2O(l)";
    MineralReaction min_reaction_daphnite = editor.addMineralReaction("Daphnite,14A")
            .setEquation(eq_str_daphnite)
            .setSpecificSurfaceArea(0.0027, "m2/g");
    Reaction reaction_daphnite = createReaction(min_reaction_daphnite, system);
    reaction_daphnite.setName("Daphnite reaction");
    ReactionRateFunction rate_func_daphnite_shell = [&min_reaction_daphnite, &reaction_daphnite, &system](const ChemicalProperties& properties) -> ChemicalScalar {

        // The number of chemical species in the system
        const unsigned num_species = system.numSpecies();

        // The mineral reaction rate using specified surface area
        ChemicalScalar res(num_species, 0.0);

        // The universal gas constant (in units of kJ/(mol*K))
        const double R = 8.3144621e-3;

        // The temperature and pressure of the system
        const Temperature T = properties.temperature();

        // Auxiliary variables
        ChemicalScalar f(num_species, 1.0);

        // Create a Reaction instance
        Reaction reaction(min_reaction_daphnite.equation(), system);

        // Calculate the saturation index of the mineral
        const auto lnK = reaction.lnEquilibriumConstant(properties);
        const auto lnQ = reaction.lnReactionQuotient(properties);
        const auto lnOmega = lnQ - lnK;

        // Calculate the saturation index
        const auto Omega = exp(lnOmega);
        //std::cout << "Omega(Daphnite) = " << Omega << std::endl;

        // The composition of the chemical system
        const auto n = properties.composition();
        const auto n0_rxn = reaction_daphnite.initialAmounts();

        // The index of the mineral
        const Index imineral = system.indexSpeciesWithError(min_reaction_daphnite.mineral());

        // The number of moles of the mineral
        auto nm = n[imineral];
        auto nm0_rxn = n0_rxn[imineral];

        // Calculate activities ofr H+ and OH- species
        VectorConstRef lna = properties.lnActivities().val;
        const Index i_h = system.indexSpeciesWithError("H+");
        double activity_h = std::exp(lna(i_h));
        const Index i_oh = system.indexSpeciesWithError("OH-");
        double activity_oh = std::exp(lna(i_oh));

        // The molar mass of the mineral (in units of kg/mol)
        const double molar_mass = system.species(imineral).molarMass();

        /*
         * Chamosite(Daphnite)
        # original compilation from 15mar/cla
        # kinetic data extracted from 09gol/ben
        # surface area data extracted from (03bra/bos)
        # warning dissolution only
        # Confidence level: 4
        -start
        1 SRmin = SR("Chamosite(Daphnite)")
        10 moles = 0
        20 If (m <= 0) and (SRmin < 1) Then GoTo 250
        30 S = 0.0027 # 0.2% BET (03bra/bos); suggested value in m2/g
        40 Mm = 713.5 # molar mass in g/mol
        50 If (SRmin > 1) Then GoTo 250
        ########## start dissolution bloc ##########
        60 knu = 6.4E-17 * exp((-16000 / 8.314) * ((1 / TK) - (1 / 298.15)))
        70 k1 = 8.2E-09 * exp((-17000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("H3O+") ^ 0.28)
        80 k2 = 6.9E-09 * exp((-16000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("OH-") ^ 0.34)
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
        if(nm.val <= 0 && Omega < 1) // the is no way to precipitate further
            res = ChemicalScalar(num_species, 0.0);
        // S = 0.0027 # 0.2% BET (03bra/bos); suggested value in m2/g
        const auto ssa = 0.0027;

        if(Omega < 1) // dissolution kinetics
        {
            /*
             * ########## start dissolution bloc ##########
            60 knu = 6.4E-17 * exp((-16000 / 8.314) * ((1 / TK) - (1 / 298.15)))
            70 k1 = 8.2E-09 * exp((-17000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("H3O+") ^ 0.28)
            80 k2 = 6.9E-09 * exp((-16000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("OH-") ^ 0.34)
            90 k = knu + k1 + k2
            100 rate = S * m * Mm * ((m/m0)^(2/3)) * k * (1 - SRmin) # by default
            110 moles = rate * Time
            120 REM Do not dissolve more than what is available
            130 IF (moles > M) THEN moles = M
            ########## end dissolution bloc ##########
             */

            // knu = 6.4E-17 * exp((-16000 / 8.314) * ((1 / TK) - (1 / 298.15)))
            const auto kappa_neu = 6.4E-17 * exp(- 16.0 / R * (1.0/T - 1.0/298.15));

            // k1 = 8.2E-09 * exp((-17000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("H3O+") ^ 0.28)
            const auto kappa_acid = 8.2E-09 * exp(- 17.0 / R * (1.0/T - 1.0/298.15)) * std::pow(activity_h, 0.28);

            // OH- catalyzer (sulfide promotion)
            // k2 = 6.9E-09 * exp((-16000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("OH-") ^ 0.34)
            const auto kappa_oh = 6.9E-09 * exp(- 16.0 / R * (1.0/T - 1.0/298.15)) * std::pow(activity_oh, 0.34);

            const auto kappa = kappa_neu + kappa_acid + kappa_oh;

            // Calculate the resulting mechanism function
            // rate = S * m * Mm * ((m/m0)^(2/3)) * k * (1 - SRmin) # by default
            res += f * ssa * nm.val * molar_mass * pow(nm.val/nm0_rxn, 2/3) * kappa * (1 - Omega);

            // Do not dissolve more than what is available
            // IF (moles > M) THEN moles = M
//            double total_moles = nm.val; // current amount of mols of available minerals
//            if (res > nm.val) res += nm.val;

        }
        //std::cout << "res = " << res << std::endl;
        //getchar();
        return res;

    };
    reaction_daphnite.setRate(rate_func_daphnite_shell);

    //Siderite
    //FeCO3 + 1.000H+ = 1.000HCO3- + 1.000Fe+2
    //     log_k    -0.273
    //     delta_h  -27.862    #kJ/mol        #Internal calculation
    //     -analytic -9.0291123E+2  -1.4586221E-1  4.9931005E+4  3.2756219E+2  -2.8333834E+6
    //     #References = LogK/DGf: 04chi; DHf/DHr: Internal calculation; S°: 04chi; Cp: 04chi; V°: 78hel/del,85hel;
    std::string eq_str_siderite = "Siderite + H+ = HCO3- + Fe++";
    MineralReaction min_reaction_siderite = editor.addMineralReaction("Siderite")
            .setEquation(eq_str_siderite)
            .setSpecificSurfaceArea(0.13, "m2/g");
    Reaction reaction_siderite = createReaction(min_reaction_siderite, system);
    reaction_siderite.setName("Siderite reaction");
    ReactionRateFunction rate_func_siderite_shell = [&min_reaction_siderite, &reaction_siderite, &system](const ChemicalProperties& properties) -> ChemicalScalar {

        // The number of chemical species in the system
        const unsigned num_species = system.numSpecies();

        // The mineral reaction rate using specified surface area
        ChemicalScalar res(num_species, 0.0), res_growth(num_species, 0.0), res_nuc(num_species, 0.0);

        // The universal gas constant (in units of kJ/(mol*K))
        const double R = 8.3144621e-3;

        // The temperature and pressure of the system
        const Temperature T = properties.temperature();

        // Auxiliary variables
        ChemicalScalar f(num_species, 1.0);

        // Create a Reaction instance
        Reaction reaction(min_reaction_siderite.equation(), system);

        // Calculate the saturation index of the mineral
        const auto lnK = reaction.lnEquilibriumConstant(properties);
        const auto lnQ = reaction.lnReactionQuotient(properties);
        const auto lnOmega = lnQ - lnK;
        //std::cout << "lnK = " << lnK << std::endl;
        //std::cout << "lnQ = " << lnQ << std::endl;
        //std::cout << "lnOmega = " << lnOmega << std::endl;

        // Calculate the saturation index
        const auto Omega = exp(lnOmega);
        //std::cout << "Omega(Siderite) = " << Omega << std::endl;

        // The composition of the chemical system
        const auto n = properties.composition();
        const auto n0_rxn = reaction_siderite.initialAmounts();

        // The index of the mineral
        const Index imineral = system.indexSpeciesWithError(min_reaction_siderite.mineral());

        // The number of moles of the mineral
        auto nm = n[imineral];
        auto nm0_rxn = n0_rxn[imineral];

        // Calculate activities ofr H+ and OH- species
        VectorConstRef lna = properties.lnActivities().val;
        const Index i_h = system.indexSpeciesWithError("H+");
        double activity_h = std::exp(lna(i_h));

        // The molar mass of the mineral (in units of kg/mol)
        const double molar_mass = system.species(imineral).molarMass();

        /*
         * Siderite
            # kinetic data extracted from 09gol/ben	92gre/tom
            # surface area data extracted from 15mar/cla
            # warning dissolution only
            # Confidence level: 4
            -start
            1 SRmin = SR("Siderite")
            10 moles = 0
            20 If (m <= 0) and (SRmin < 1) Then GoTo 350
            25 If (SRmin = 1) Then GoTo 350
            30 S = 0.13 # average BET; suggested value in m2/g
            40 Mm = 115.86	# molar mass in g/mol
            45 Mv = 2.049E-05									# molar volume in m3/mol
            50 If (SRmin > 1) Then GoTo 130
            ########## start dissolution bloc ##########
            60 knu = 2.1E-09 * exp((-56000 / 8.314) * ((1 / TK) - (1 / 298.15)))
            70 k1 = 5.9E-06 * exp((-56000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("H3O+") ^ 0.60)
            80 k = knu + k1 + k2
            90 rate = S * m * Mm * ((m/m0)^(2/3)) * k * (1 - SRmin)
            100 moles = rate * Time
            110 IF (moles > M) THEN moles = M
            120 GoTo 350
            ########## end dissolution bloc ##########
            ########## start precipitation bloc ##########
            130 surf_energy = 0.06 									# interfacial energy in J m-2
            140 molvol = 3.4e-29									# molecular volume in m3
            150 u = ( 16 * 3.14159 * surf_energy^3 * molvol ^ 2 ) / ( 3 * (1.38E-23 * TK)^3 ) 	# (after Fritz et al 2009)
            160 J0 = 1E20  											# nucleation rate in nuclei kg-1 sec-1 (after Fritz et al 2009)
            170 SRcrit = exp( sqrt( u / LOG( J0 ) ) ) 				# critical saturation threshold
            180 IF ( SRmin < SRcrit ) Then GoTo 250 				# if saturation index below threshold no nucleation occurs
            190 nuclie = J0 * exp( -u / LOG(SRmin)^2 ) * time		# --- nucleation rate in nuclie sec-1 multiplied by time to account for nucleation over time
            200 IF ( nuclie < 1 ) THEN GoTo 250 					# condition that rate needs to be bigger than 1 nucl per sec (Arbitrary)
            210 nj = ( 2 * u ) / (LOG(SRmin)^3) 					# number of growth units in critical nuclie radius
            220 vol = nj * molvol 									# critical nuclie volume
            230 moles_nuc = ( -nuclie * vol ) / Mv 	# moles of nuclie formed
            240 GoTo 260
            250 moles_nuc = 0
            #### kinetic data for growth ####
            260 IF ( m = 0 ) THEN GoTo 300
            270 knu = 1.6E-11 * exp((-108000 / 8.314) * ((1 / TK) - (1 / 298.15)))
            280 kpre = (-1) * knu
            290 moles_growth = S * m * Mm * kpre * ABS(SRmin - 1) * Time
            300 moles = moles_nuc + moles_growth
            310 IF ( moles = 0 ) THEN GoTo 350
            ########## end precipitation bloc ##########
            320 maxMol = TOT("Fe(2)")
            330 IF (maxMol > TOT("C")) THEN maxMol = TOT("C")
            340 IF (maxMol < -moles) THEN moles = -maxMol
            ########## end precipitation bloc ##########
            350 Save moles
            -end
         */

        // If (m <= 0) and (SRmin < 1) Then GoTo 350
        // If (SRmin = 1) Then GoTo 350
        if((nm.val <= 0 && Omega < 1) || Omega == 1) // the is no way to precipitate further
            return res;

        // S = 0.13 # average BET; suggested value in m2/g
        const auto ssa = 0.13;
        // Mv = 2.049E-05
        const auto molar_volume = 2.049E-05;

        if(Omega < 1) // dissolution kinetics
        {
            /*
             * 60 knu = 2.1E-09 * exp((-56000 / 8.314) * ((1 / TK) - (1 / 298.15)))
                70 k1 = 5.9E-06 * exp((-56000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("H3O+") ^ 0.60)
                80 k = knu + k1 + k2
                90 rate = S * m * Mm * ((m/m0)^(2/3)) * k * (1 - SRmin)
                100 moles = rate * Time
             */

            // knu = 2.1E-09 * exp((-56000 / 8.314) * ((1 / TK) - (1 / 298.15)))
            const auto kappa_neu = 2.1E-09 * exp(- 56.0 / R * (1.0/T - 1.0/298.15));

            // k1 = 5.9E-06 * exp((-56000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("H3O+") ^ 0.60)
            const auto kappa_acid = 5.9E-06 * exp(- 56.0 / R * (1.0/T - 1.0/298.15)) * std::pow(activity_h, 0.6);

            const auto kappa = kappa_neu + kappa_acid;

            // Calculate the resulting mechanism function
            // rate = S * m * Mm * ((m/m0)^(2/3)) * k * (1 - SRmin) # by default
            res += f * ssa * nm.val * molar_mass * pow(nm.val/nm0_rxn, 2/3) * kappa * (1 - Omega);

            // Do not dissolve more than what is available
            // IF (moles > M) THEN moles = M
//            double total_moles = nm.val; // current amount of mols of available minerals
//            if (res > nm.val) res += nm.val;

        }
        else // precipitation kinetics
        {
            /*
             * ########## start precipitation bloc ##########
                130 surf_energy = 0.06 									# interfacial energy in J m-2
                140 molvol = 3.4e-29									# molecular volume in m3
                150 u = ( 16 * 3.14159 * surf_energy^3 * molvol ^ 2 ) / ( 3 * (1.38E-23 * TK)^3 ) 	# (after Fritz et al 2009)
                160 J0 = 1E20  											# nucleation rate in nuclei kg-1 sec-1 (after Fritz et al 2009)
                170 SRcrit = exp( sqrt( u / LOG( J0 ) ) ) 				# critical saturation threshold
                180 IF ( SRmin < SRcrit ) Then GoTo 250 				# if saturation index below threshold no nucleation occurs
                190 nuclie = J0 * exp( -u / LOG(SRmin)^2 ) * time		# --- nucleation rate in nuclie sec-1 multiplied by time to account for nucleation over time
                200 IF ( nuclie < 1 ) THEN GoTo 250 					# condition that rate needs to be bigger than 1 nucl per sec (Arbitrary)
                210 nj = ( 2 * u ) / (LOG(SRmin)^3) 					# number of growth units in critical nuclie radius
                220 vol = nj * molvol 									# critical nuclie volume
                230 moles_nuc = ( -nuclie * vol ) / Mv 	# moles of nuclie formed
                240 GoTo 260
                250 moles_nuc = 0
                #### kinetic data for growth ####
                260 IF ( m = 0 ) THEN GoTo 300
                270 knu = 1.6E-11 * exp((-108000 / 8.314) * ((1 / TK) - (1 / 298.15)))
                280 kpre = (-1) * knu
                290 moles_growth = S * m * Mm * kpre * ABS(SRmin - 1) * Time
                300 moles = moles_nuc + moles_growth
                310 IF ( moles = 0 ) THEN GoTo 350
                ########## end precipitation bloc ##########
             */

            // surf_energy = 0.06 # interfacial energy in J m-2
            const auto surf_energy = 0.06;

            // molvol = 3.4e-29	# molecular volume in m3
            const auto molecular_volume = 3.4e-29;

            // u = ( 16 * 3.14159 * surf_energy^3 * molvol ^ 2 ) / ( 3 * (1.38E-23 * TK)^3 ) 	# (after Fritz et al 2009)
            const auto u = 16 * 3.14159 * std::pow(surf_energy, 3.0) * std::pow(molecular_volume, 2.0) / (3 * std::pow(1.38E-23 * T.val, 3.0));

            // J0 = 1E20 # nucleation rate in nuclei kg-1 sec-1 (after Fritz et al 2009)
            const auto j0 = 1E20;

            // SRcrit = exp( sqrt( u / LOG( J0 ) ) ) # critical saturation threshold
            const auto Omega_crit = std::exp(sqrt(u / std::log10(j0)));

            if(Omega >= Omega_crit)
            {
                // knu = 1.6E-11 * exp((-108000 / 8.314) * ((1 / TK) - (1 / 298.15)))
                // kpre = (-1) * knu
                const auto kappa_growth = -1.6E-11 * exp(- 108 / R * (1.0 / T - 1.0 / 298.15));

                // rate = S * m * Mm * kpre * (SRmin - 1)
                res_growth += f * ssa * nm.val * molar_mass * kappa_growth * std::abs(Omega.val - 1);

            }
            else
            {
                // nuclie = J0 * exp( -u / LOG(SRmin)^2 ) * time # --- nucleation rate in nuclie sec-1 multiplied by time to account for nucleation over time
                const auto nuclie = j0 * exp(- u / std::pow(std::log10(Omega.val), 2.0));

                // IF ( nuclie < 1 ) THEN GoTo 250      # condition that rate needs to be bigger than 1 nucl per sec (Arbitrary)
                // nj = ( 2 * u ) / (LOG(SRmin)^3) 		# number of growth units in critical nuclie radius
                const auto nj = 2 * u / std::pow(std::log10(Omega.val), 3.0);

                // vol = nj * molvol 					# critical nuclie volume
                const auto vol = nj * molecular_volume;

                // moles_nuc = ( -nuclie * vol ) / Mv 	# moles of nuclie formed
                res_nuc += - nuclie * vol / molar_volume;
            }

            // moles = moles_nuc + moles_growth
            res += res_growth + res_nuc;

        }
        //std::cout << "res = " << res << std::endl;
        //getchar();
        return res;

    };
    reaction_siderite.setRate(rate_func_siderite_shell);


    // Kaolinite: Al2Si2O5(OH)4
    //std::string eq_str_kaolinite = "Kaolinite + 4*H+ = 3*H2O(l) + 2*HSiO3- + 2*Al+++";
    // Al2Si2O5(OH)4 + 6 H+  =  2 Al+++ + 2 SiO2 + 5 H2O
    std::string eq_str_kaolinite = "Kaolinite + 6*H+ = 2*Al+++ + 2*SiO2(aq) + 5*H2O(l)";
    MineralReaction min_reaction_kaolinite = editor.addMineralReaction("Kaolinite")
            .setEquation(eq_str_kaolinite)
            .addMechanism("logk = -13.18 mol/(m2*s); Ea = 22.2 kJ/mol")
            .addMechanism("logk = -11.31 mol/(m2*s); Ea = 65.9 kJ/mol; a[H+] = 0.777")
            .setSpecificSurfaceArea(1, "m2/kg");
    Reaction reaction_kaolinite = createReaction(min_reaction_kaolinite, system);
    reaction_kaolinite.setName("Kaolinite reaction");
    ReactionRateFunction rate_func_kaolinite = [&min_reaction_kaolinite, &reaction_kaolinite, &system, &I](const ChemicalProperties& properties) -> ChemicalScalar {

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

        // Create a Reaction instance
        Reaction reaction(min_reaction_kaolinite.equation(), system);

        // Calculate the saturation index of the mineral
        const auto lnK = reaction.lnEquilibriumConstant(properties);
        const auto lnQ = reaction.lnReactionQuotient(properties);
        const auto lnOmega = lnQ - lnK;

        const auto lnK_rxn = reaction_kaolinite.lnEquilibriumConstant(properties);
        const auto lnQ_rxn = reaction_kaolinite.lnReactionQuotient(properties);
        const auto lnOmega_rxn = lnQ_rxn - lnK_rxn;

        // Calculate the saturation index
        //auto Omega = exp(lnOmega);
        //auto Omega_rxn = exp(lnOmega_rxn);
        auto Omega = exp(lnOmega).val;
        auto Omega_rxn = exp(lnOmega_rxn).val;

        // The composition of the chemical system
        const auto n = properties.composition();
        // The initial composition of the chemical system
        const auto n0_rxn = reaction_kaolinite.initialAmounts();

        // Amount of elements
        const auto b = system.elementAmounts(n.val);

        // Calculate TOT("Si")
        const Index i_si = system.indexElementWithError("Si");
        const auto tot_si = b(i_si);

        // Calculate TOT("Al")
        const Index i_al = system.indexElementWithError("Al");
        const auto tot_al = b(i_al);

        // The index of the mineral
        const Index imineral = system.indexSpeciesWithError(min_reaction_kaolinite.mineral());

        // The current and the initial number of moles of the mineral
        auto nm = n[imineral];
        auto nm0_rxn = n0_rxn[imineral];

        VectorConstRef lna = properties.lnActivities().val;
        const Index i_h = system.indexSpeciesWithError("H+");
        double activity_h = std::exp(lna(i_h));
        double ionic_strength = I(properties).val;

        // The molar mass of the mineral (in units of kg/mol)
        const double molar_mass = system.species(imineral).molarMass();

        // Specific surface area
        const auto ssa = min_reaction_kaolinite.specificSurfaceArea();

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
        if(nm.val == 0 && lnOmega < sat_index_threshold)
            return res;
        // IF (M0 > 0) AND (M > 0) THEN t = (M / M0)^PARM(4)
        // PARM(4): exponent for M/M0 for surface area correction (=1 if no data)
        if(nm0_rxn > 0 && nm.val > 0)
            ssa_corr_factor = pow(nm.val/nm0_rxn, 2/3);

        // IF (M0 > 0) AND (M > 0) THEN t = (M / M0)^PARM(4)
        // PARM(4): exponent for M/M0 for surface area correction (=1 if no data)
        if(nm0_rxn > 0 && nm.val > 0)
            ssa_corr_factor = pow(nm.val/nm0_rxn, 2/3);

        // Acid mechanism
        // k25a =  0
        // Eaa =   0
        // na = 0
        // IF (k25a > 0) THEN KTa = k25a * EXP((-Eaa/R) * (1/TK - 1/298.15)) * ACT("H+")^na
        const auto kappa_acid = 0.0;

        // Neutral mechanism
        // k25n =   1E-13
        // Ean =    7100
        // IF (k25n>0) THEN KTn = k25n * EXP((-Ean/R) * (1/TK - 1/298.15))
        const auto kappa_neu = 1E-13 * exp(- 7.1 / R * (1.0/T - 1.0/298.15));

        // Base mechanism
        // k25b =   0
        // Eab =    0
        // nb = 0
        const auto kappa_base = 0.0;

        // Define dissolution and precipitation rates
        auto kappa_diss = kappa_neu + kappa_acid + kappa_base;
        auto kappa_pre = kappa_neu + kappa_acid + kappa_base;
        /*
         * Kaolinite
            -start
            10 k25a =	0
            20 k25n =	1E-13
            30 k25b =	0
            40 Eaa =	0
            50 Ean =	7100
            60 Eab =	0
            70 na =	0
            80 nb =	0
            90  SatIndex = SI("Kaolinite")
            100 SatRatio = SR("Kaolinite")
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
            250   IF (k25b > 0) THEN KTb = k25b * EXP((-Eab/R) * (1/TK - 1/298.15)) * ACT("H+")^nb
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
            610   maxMol = TOT("Si")/2
            620   IF (maxMol > TOT("Al")/2) THEN maxMol =TOT("Al")/2
            690   IF (maxMol < -moles) THEN moles = -maxMol
            900 SAVE moles
            -end
         */
        if(Omega <= 1) // dissolution kinetics
        {
            // moles = PARM(1) * M0 * t * (KTa+KTn+KTb) * (1 - SatRatio) * TIME
            res += f * ssa * nm0_rxn * ssa_corr_factor * molar_mass * kappa_pre * (1 - Omega);

            // Do not dissolve more than what is available
            // 320 IF (SatIndex > 0) THEN GOTO 400
            // 330   REM Do not dissolve more than what is available
            // 340   IF (moles > M) THEN moles = M
            // 350   GOTO 900
//            double total_moles = nm.val; // current amount of mols of available minerals
//            if (lnOmega <= 0 && res > nm.val)
//                res += nm.val;
        }
        else // precipitation kinetics
        {
            // REM Precipitation when M0 = 0 (secondary mineral)
            if(nm0_rxn != 0)
            {
                // Default nucleation surface area
                auto nucl_sa = 1.0e-5;
                // IF (M > 0) THEN surfArea = M * PARM(3)
                if(nm.val) nucl_sa *= nm.val;

                // moles = surfArea * (KTa+KTn+KTb)* PARM(5) * (1 - SatRatio) * TIME
                res += f * nucl_sa * kappa_pre * (1 - Omega);
            }
            else // IF (M0 = 0) THEN GOTO 500
                // Set nucleation rate
                // moles = moles * PARM(5)
                res *= f * kappa_pre / kappa_diss;

            // 600 REM Do not precipitate more than the elements in solution
            // 610   maxMol = TOT("Si")/2
            // 620   IF (maxMol > TOT("Al")/2) THEN maxMol =TOT("Al")/2
            // 690   IF (maxMol < -moles) THEN moles = -maxMol
            // 900 SAVE moles

            // TODO: implement if in the kinetic solver
//            // Implement upper bound in precipitation kinetics
//            auto max_mol = tot_si;
//            if(max_mol > tot_al) max_mol = tot_al;
//            if(max_mol < -nm.val) return ChemicalScalar(num_species, 0.0)

        }

        return res;

    };
    reaction_kaolinite.setRate(rate_func_kaolinite);

    std::string eq_str_quartz = "Quartz + H2O(l) = H+ + HSiO3-";
    //std::string eq_str_quartz = "Quartz = SiO2(aq)";
    MineralReaction min_reaction_quartz = editor.addMineralReaction("Quartz")
            .setEquation(eq_str_quartz)
            .addMechanism("logk = -13.99 mol/(m2*s); Ea = 87.7 kJ/mol")
            .setSpecificSurfaceArea(0.1, "m2/kg");
    Reaction reaction_quartz = createReaction(min_reaction_quartz, system);
    reaction_quartz.setName("Quartz reaction");
    ReactionRateFunction rate_func_quartz_llnl = [&min_reaction_quartz, &reaction_quartz, &system, &I](const ChemicalProperties& properties) -> ChemicalScalar {

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

        // Create a Reaction instance
        Reaction reaction(min_reaction_quartz.equation(), system);

        // Calculate the saturation index of the mineral
        const auto lnK = reaction.lnEquilibriumConstant(properties);
        const auto lnQ = reaction.lnReactionQuotient(properties);
        const auto lnOmega = lnQ - lnK;

        const auto lnK_rxn = reaction_quartz.lnEquilibriumConstant(properties);
        const auto lnQ_rxn = reaction_quartz.lnReactionQuotient(properties);
        const auto lnOmega_rxn = lnQ_rxn - lnK_rxn;

        // Calculate the saturation index
        //auto Omega = exp(lnOmega);
        //auto Omega_rxn = exp(lnOmega_rxn);
        auto Omega = exp(lnOmega).val;
        auto Omega_rxn = exp(lnOmega_rxn).val;

        // The composition of the chemical system
        const auto n = properties.composition();
        // The initial composition of the chemical system
        const auto n0_rxn = reaction_quartz.initialAmounts();

        // Amount of elements
        const auto b = system.elementAmounts(n.val);

        // Calculate TOT("Si")
        const Index i_si = system.indexElementWithError("Si");
        const auto tot_si = b(i_si);

        // The index of the mineral
        const Index imineral = system.indexSpeciesWithError(min_reaction_quartz.mineral());

        // The current and the initial number of moles of the mineral
        auto nm = n[imineral];
        auto nm0_rxn = n0_rxn[imineral];

        VectorConstRef lna = properties.lnActivities().val;
        const Index i_h = system.indexSpeciesWithError("H+");

        // The molar mass of the mineral (in units of kg/mol)
        const double molar_mass = system.species(imineral).molarMass();

        // Specific surface area
        const auto ssa = min_reaction_quartz.specificSurfaceArea();

        // Specific surface area correction factor
        auto ssa_corr_factor = 1.0;

        // REM If SI too small, then No dissolution nor precipitation
        // IF (SatIndex > -1e-15) AND ( SatIndex < 1e-15) THEN GOTO 900
        if(lnOmega > -1e-15 && lnOmega < 1e-15)
            return res;
        //  IF (SatIndex > 0) THEN GOTO 900	# IDISPRE=1 in TR model (i.e. only dissolution allowed)
        if(lnOmega > 0)
            return res;
        // IF (SatRatio > 10000) THEN SatRatio = 10000
        if(Omega > 1e4)
            Omega = 1e4;
        // IF (M = 0) AND (SatIndex < PARM(2)) THEN GOTO 900
        // PARM(2): used when M0 = 0, (Nucleation) SatIndex threshold for precipitation
        double sat_index_threshold = 1e6;
        if(nm.val == 0 && lnOmega < sat_index_threshold)
            return res;
        // IF (M0 > 0) AND (M > 0) THEN t = (M / M0)^PARM(4)
        // PARM(4): exponent for M/M0 for surface area correction (=1 if no data)
        if(nm0_rxn > 0 && nm.val > 0)
            ssa_corr_factor = pow(nm.val/nm0_rxn, 2/3);

        // k25a =   0
        // Eaa =    0
        // na = 0
        // IF (k25a > 0) THEN KTa = k25a * EXP((-Eaa/R) * (1/TK - 1/298.15)) * ACT("H+")^na
        const auto kappa_acid = 0.0;
        // k25n =	4.52E-14
        // Ean =	90100
        // IF (k25n>0) THEN KTn = k25n * EXP((-Ean/R) * (1/TK - 1/298.15))
        const auto kappa_neu = 4.52E-14 * exp(- 90.1 / R * (1.0/T - 1.0/298.15));
        // k25b =   0
        // Eab =    0
        // nb = 0
        const auto kappa_base = 0.0;

        // Precipitation and dissolution constant rates
        auto kappa_diss = kappa_neu + kappa_acid + kappa_base;
        auto kappa_pre = kappa_neu + kappa_acid + kappa_base;

        /*
         * Quartz
        -start
        10 k25a =	0
        20 k25n =	4.52E-14 # NB we have not implemented the rate correction below pH=3 and pH=8 (mineral record 3 in TR model), since pH stays within this interval
        30 k25b =	0
        40 Eaa =	0
        50 Ean =	90100
        60 Eab =	0
        70 na =	0
        80 nb =	0
        90  SatIndex = SI("Quartz")
        100 SatRatio = SR("Quartz")
        110 R = 8.314472 # (gas constant J/mol/Kelvin)
        120 REM If SI too small, then No dissolution nor precipitation
        130 IF (SatIndex > -1e-15) AND ( SatIndex < 1e-15) THEN GOTO 900
        132 IF (SatIndex > 0) THEN GOTO 900	# IDISPRE=1 in TR model (i.e. only dissolution allowed)
        135 IF (SatRatio > 10000) THEN SatRatio = 10000
        140 IF (M = 0) AND (SatIndex < PARM(2)) THEN GOTO 900
        160 t = 1    # surface area correction factor
        170 IF (M0 > 0) AND (M > 0) THEN t = (M / M0)^PARM(4)
        200 REM acid mechanism
        210   IF (k25a > 0) THEN KTa = k25a * EXP((-Eaa/R) * (1/TK - 1/298.15)) * ACT("H+")^na
        220 REM neutral mechanism
        230   IF (k25n>0) THEN KTn = k25n * EXP((-Ean/R) * (1/TK - 1/298.15))
        240 REM base mechanism
        250   IF (k25b > 0) THEN KTb = k25b * EXP((-Eab/R) * (1/TK - 1/298.15)) * ACT("H+")^nb
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
        610   maxMol = TOT("Si")
        690   IF (maxMol < -moles) THEN moles = -maxMol
        900 SAVE moles
        -end
         */

        if(Omega <= 1) // dissolution kinetics
        {
            // moles = PARM(1) * M0 * t * (KTa+KTn+KTb) * (1 - SatRatio) * TIME
            res += f * ssa * nm0_rxn * ssa_corr_factor * molar_mass * kappa_pre * (1 - Omega);

            // Do not dissolve more than what is available
            // 320 IF (SatIndex > 0) THEN GOTO 400
            // 330   REM Do not dissolve more than what is available
            // 340   IF (moles > M) THEN moles = M
            // 350   GOTO 900
//            double total_moles = nm.val; // current amount of mols of available minerals
//            if (lnOmega <= 0 && res > nm.val)
//                res += nm.val;
        }
        else // precipitation kinetics
        {
            // REM Precipitation when M0 = 0 (secondary mineral)
            if(nm0_rxn != 0)
            {
                // Default nucleation surface area
                auto nucl_sa = 1.0e-5;
                // IF (M > 0) THEN surfArea = M * PARM(3)
                if(nm.val) nucl_sa *= nm.val;

                // moles = surfArea * (KTa+KTn+KTb)* PARM(5) * (1 - SatRatio) * TIME
                res += f * nucl_sa * kappa_pre * (1 - Omega);
            }
            else // IF (M0 = 0) THEN GOTO 500
                // Set nucleation rate
                // moles = moles * PARM(5)
                res *= f * kappa_pre / kappa_diss;

            //     600 REM Do not precipitate more than the elements in solution
            // 610   maxMol = TOT("Si")
            // 690   IF (maxMol < -moles) THEN moles = -maxMol
            // 900 SAVE moles

            // TODO: implement if in the kinetic solver
//            // Implement upper bound in precipitation kinetics
//            auto max_mol = tot_si;
//            if(max_mol < -nm.val) return ChemicalScalar(num_species, 0.0);

        }

        return res;

    };
    reaction_quartz.setRate(rate_func_quartz_llnl);

    // Step **: Create the ReactionSystem instances
    ReactionSystem reactions(editor);
    //ReactionSystem reactions(system, {reaction_siderite});
    //ReactionSystem reactions(system, {reaction_siderite, reaction_daphnite});
    //ReactionSystem reactions(system, {reaction_siderite, reaction_daphnite});

    //std::cout << "system = \n" << system << std:: endl;
    // getchar();
    //std::cout << "reactions number = " << reactions.numReactions() << std:: endl;

    // Step **: Create the ReactionSystem instances
    Partition partition(system);
    //partition.setKineticSpecies(std::vector<std::string>{"Daphnite,14A"});
    //partition.setKineticSpecies(std::vector<std::string>{"Siderite", "Daphnite,14A"});

    /*
     * SOLUTION 0
        temp      25
        pH        7.1015
        pe        -3.6887
        redox     S(-2)/S(6)
        units     mg/l
        density   1.02368
        Alkalinity 1441.11 as HCO3
        Ca        310
        Cl        19100 charge
        Fe(2)     0.0034
        Mg        1240
        Na        11400
        S(-2)     100
        S(6)      2420
        -water    1 # kg

    EQUILIBRIUM_PHASES 1
        Pyrrhotite 0 0
        Calcite 0 0
     */
    // Step **: Define the initial condition (IC) of the reactive transport modeling problem
    EquilibriumInverseProblem problem_ic(partition);
    problem_ic.setTemperature(params.T, "celsius");
    problem_ic.setPressure(params.P, "bar");
    problem_ic.add("H2O", 1.0, "kg");
    problem_ic.add("Cl-", 19100, "mg"); // 19100 charge
    problem_ic.add("Fe++", 1240, "mg"); // 1240
    problem_ic.add("Na+", 11400, "mg"); // 11400
    problem_ic.add("2420", 157.18e-3, "mg"); // 2420 mg/l
    problem_ic.add("Mg++", 1240, "mg"); // 1240 mg/l
    problem_ic.add("Ca++", 23.838e-3, "mg"); // 310 mg/l
    problem_ic.add("HCO3-", 1441.11, "mg"); // 1441.11 mg/l
    problem_ic.add("S", 100, "mg"); // 100 mg/l
    problem_ic.add("Calcite", 0.0, "mol"); // Calcite 0 0
    problem_ic.add("Pyrrhotite", 0.0, "mol"); // Pyrrhotite 0 0
    problem_ic.pH(7.1015);

    ChemicalState state_ic = equilibrate(problem_ic);

    //std::cout << "state_ic = \n" << state_ic << std:: endl;
    //getchar();

    // Step **: Define the boundary condition (BC)  of the reactive transport modeling problem
    EquilibriumInverseProblem problem_bc(partition);
    problem_bc.setTemperature(params.T, "celsius");
    problem_bc.setPressure(params.P, "bar");
    problem_bc.add("H2O", 1.0, "kg");
    problem_bc.add("Cl-", 19100, "mg"); // 19100 charge
    problem_bc.add("Fe++", 1240, "mg"); // 1240
    problem_bc.add("Na+", 11400, "mg"); // 11400
    problem_bc.add("2420", 157.18e-3, "mg"); // 2420 mg/l
    problem_bc.add("Mg++", 1240, "mg"); // 1240 mg/l
    problem_bc.add("Ca++", 23.838e-3, "mg"); // 310 mg/l
    problem_bc.add("HCO3-", 1441.11, "mg"); // 1441.11 mg/l
    problem_bc.add("S", 100, "mg"); // 100 mg/l
    problem_bc.add("HS-", 0.0196504 / 58, "mol");
    problem_bc.add("H2S(aq)", 0.167794 / 58, "mol");
    problem_bc.pH(5.726);

    // Step **: Calculate the equilibrium states for the IC and BC
    ChemicalState state_bc = equilibrate(problem_bc);

    // Step **: Scale the boundary condition state
    state_bc.scaleVolume(1.0, "m3");

    // Set minerals amount
    // Chamosite(Daphnite)
    //    -formula  Chamosite(Daphnite)  1
    //    -m        0.1
    //    -m0       0.1
    //    -tol      1e-008
    state_ic.setSpeciesAmount("Daphnite,14A", 0.1, "mol"); // MM(Daphnite) = 713.5 g / mol
    state_ic.setSpeciesAmount("Siderite", 0.1, "mol"); // MM(Siderite) = 115.86 g / mol
    state_ic.setSpeciesAmount("Quartz", 0.1, "mol"); // MM(Quartz) = 60.083 g / mol
    state_ic.setSpeciesAmount("Kaolinite", 0.1, "mol"); // MM(Kaolinite) = 258.071 g / mol

    // Step **: Scale the volumes of the phases in the initial condition
    state_ic.scalePhaseVolume("Aqueous", 0.1, "m3");    // 10% if the 1.0m3
    state_ic.scaleVolume(1.0, "m3");
    //std::cout << "state_ic = \n" << state_ic << std:: endl;
    //getchar();

    // Set initial value of Daphnite
    reaction_daphnite.setInitialAmounts(state_ic.speciesAmounts());
    reaction_siderite.setInitialAmounts(state_ic.speciesAmounts());

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
        output.add("speciesAmount(Fe+++)");
        output.add("speciesAmount(Pyrrhotite)");
        output.add("phaseMass(Pyrrhotite)");
        output.add("speciesAmount(Calcite)");
        output.add("phaseMass(Calcite)");
        output.add("speciesAmount(Pyrite)");
        output.add("phaseMass(Pyrite)");
        output.add("speciesAmount(Daphnite,14A)");
        output.add("phaseMass(Daphnite,14A)");
        output.add("speciesAmount(Siderite)");
        output.add("phaseMass(Siderite)");
        output.add("speciesAmount(Kaolinite)");
        output.add("phaseMass(Kaolinite)");
        output.add("speciesAmount(Quartz)");
        output.add("phaseMass(Quartz)");
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

auto outputStatisticsSmartKineticsConventionalEquilibrium(const Params& params, Results& results) -> void
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

auto outputStatisticsSmartKineticsSmartEquilibrium(const Params& params, Results& results) -> void
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

auto outputStatisticsConventionalKineticsConventionalEquilibrium(const Params& params, Results& results) -> void
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

auto outputStatisticsConventionalKineticsSmartEquilibrium(const Params& params, Results& results) -> void
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

    //std::string tag = "../plotting-results-20.10.20/rt-scavenging-complex-shell-siderite"; // -kin-clustering-eq-clustering";
    //std::string tag = "../plotting-results-20.10.20/rt-scavenging-complex-shell-daphnite"; // -kin-clustering-eq-clustering";
    std::string tag = "../plotting-results-20.10.20/rt-scavenging-complex-shell-siderite-daphnite-eq-only"; // -kin-clustering-eq-clustering";
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

