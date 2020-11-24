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
    // Discretisation params
    int ncells = 0; // the number of cells in the spacial discretization
    int nsteps = 0; // the number of steps in the reactive transport simulation
    double xl = 0; // the x-coordinates of the left boundaries
    double xr = 0; // the x-coordinates of the right boundaries
    double dx = 0; // the space step (in units of m)
    double dt = 0; // the time step (in units of s)

    // Physical params
    double D = 0; // the diffusion coefficient (in units of m2/s)
    double v = 0; // the Darcy velocity (in units of m/s)
    double T = 0; // the temperature (in units of degC)
    double P = 0; // the pressure (in units of bar)

    // Solver params
    bool use_smart_equilibrium_solver = false;
    double smart_equilibrium_reltol = 0.0;
    double amount_fraction_cutoff = 0.0;
    double mole_fraction_cutoff = 0.0;

    std::string activity_model = "";

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
    int minute = 60;
    int hour = 60 * minute;
    int day = 24 * hour;
    int week = 7 * day;

    // Step 2: Define parameters for the reactive transport simulation
    Params params;

    // Define discretization parameters
    params.xl = 0.0; // the x-coordinates of the left boundaries
    /*
    params.xr = 0.1; // the x-coordinates of the right boundaries
    params.ncells = 10; // the number of cells in the spacial discretization
    */
    ///*
    params.xr = 1.0; // the x-coordinates of the right boundaries
    params.ncells = 100; // the number of cells in the spacial discretization
    //params.nsteps = 1000; // the number of steps in the reactive transport simulation
    params.nsteps = 100; // the number of steps in the reactive transport simulation
    //params.nsteps = 100; // the number of steps in the reactive transport simulation
    params.dx = (params.xr - params.xl) / params.ncells; // the time step (in units of s)
    params.dt = 30 * minute; // the time step (in units of s)

    // Define physical and chemical parameters
    params.D = 1.0e-9;     // the diffusion coefficient (in units of m2/s)
    params.v = 1.0 / week; // the Darcy velocity (in units of m/s)
    params.T = 300;                     // the temperature (in units of degC)
    params.P = 85.88;                      // the pressure (in units of bar)
// Step 1: Initialise auxiliary time-related constants
//    int minute = 60;
//    int hour = 60 * minute;
//    int day = 24 * hour;
//
//    // Step 2: Define parameters for the reactive transport simulation
//    Params params;
//
//    // Define discretization parameters
//    params.xl = 0.0; // the x-coordinates of the left boundaries
//    params.xr = 100.0; // the x-coordinates of the right boundaries
//    params.ncells = 100; // the number of cells in the spacial discretization
//    params.nsteps = 50; // the number of steps in the reactive transport simulation
//    params.dx = (params.xr - params.xl) / params.ncells; // the time step (in units of s)
//    params.dt = 0.05*day; // the time step (in units of s)
//
//    // Define physical and chemical parameters
//    params.D = 0.0;     // the diffusion coefficient (in units of m2/s)
//    params.v = 1.05e-5; // the Darcy velocity (in units of m/s)
//    params.T = 300;                     // the temperature (in units of degC)
//    params.P = 85.88;                      // the pressure (in units of bar)

    // Define parameters of the equilibrium solvers
    params.smart_equilibrium_reltol = 1e-3;  // priority-based search with potential

    //params.activity_model = "hkf-full";
    params.activity_model = "hkf-selected-species";
    //params.activity_model = "pitzer-full";
    //params.activity_model = "pitzer-selected-species";
    //params.activity_model = "dk-full";
    //params.activity_model = "dk-selected-species";

    params.amount_fraction_cutoff = 1e-14;
    params.mole_fraction_cutoff = 1e-14;

    // Output
    outputConsole(params);

    // Results
    Results results;

    // Execute reactive transport with different solvers
    //params.use_smart_equilibrium_solver = true; runReactiveTransport(params, results);
    params.use_smart_equilibrium_solver = false; runReactiveTransport(params, results);

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

    // Step **: Define smart chemical equilibrium solver options
    SmartEquilibriumOptions smart_equilibrium_options;
    smart_equilibrium_options.reltol = params.smart_equilibrium_reltol;
    smart_equilibrium_options.amount_fraction_cutoff = params.amount_fraction_cutoff;
    smart_equilibrium_options.mole_fraction_cutoff = params.mole_fraction_cutoff;

    // Step **: Construct the chemical system with its phases and species (using ChemicalEditor)
    ChemicalEditor editor;

    StringList selected_elements = "Al Cl H K Na O Si";
//    StringList selected_species = "H2O(l) H+ OH- Cl- HCl(aq) Na+ NaOH(aq) NaHSiO3(aq) NaCl(aq) NaAl(OH)4(aq) "
//                                  "K+ KOH(aq) KCl(aq) KAl(OH)4(aq) Al+++ AlOH++ Al(OH)2+ Al(OH)3(aq) Al(OH)4-";
    StringList selected_species = "H2O(l) H+ OH- Cl- HCl(aq) Na+ NaOH(aq) NaHSiO3(aq) NaCl(aq) "
                                  "K+ KOH(aq) KCl(aq) Al+++ AlOH++";

    if(params.activity_model == "hkf-full"){
        // HKF full system
        editor.addAqueousPhaseWithElements(selected_elements); // total 40 species
    }
    else if(params.activity_model == "hkf-selected-species"){
        // HKF selected species
        editor.addAqueousPhaseWithElementsOf(selected_species);
    }
    else if(params.activity_model == "pitzer-full"){
        // Pitzer full system
        editor.addAqueousPhaseWithElements(selected_elements)
                .setChemicalModelPitzerHMW()
                .setActivityModelDrummondCO2();
    }
    else if(params.activity_model == "pitzer-selected-species"){
        // Pitzer selected species
        editor.addAqueousPhaseWithElementsOf(selected_species)
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
        editor.addAqueousPhaseWithElementsOf(selected_species)
                .setChemicalModelDebyeHuckel()
                .setActivityModelDrummondCO2();
    }
    editor.addMineralPhase("Quartz"); // SiO2
    editor.addMineralPhase("Diaspore"); // AlO(OH)
    editor.addMineralPhase("Gibbsite"); // Al(OH)3
    editor.addMineralPhase("Andalusite"); // Al2SiO5
    editor.addMineralPhase("Kyanite"); // Al2SiO5
    editor.addMineralPhase("Sillimanite"); // Al2SiO5
    editor.addMineralPhase("Muscovite"); // KAl2(AlSi3)O10(OH)2
    editor.addMineralPhase("Paragonite"); // NaAl2(AlSi3)O10(OH)2
    editor.addMineralPhase("Pyrophyllite"); // Al2Si4O10(OH)2
    editor.addMineralPhase("Kaolinite"); // Al2Si2O5(OH)4
    editor.addMineralPhase("Albite"); // Na(AlSi3)O8
    editor.addMineralPhase("K-Feldspar"); // K(AlSi3)O8

    // Step **: Create the ChemicalSystem object using the configured editor
    ChemicalSystem system(editor);
    //std::cout << "system = \n" << system << std:: endl;

    // Step **: Define the initial condition (IC) of the reactive transport modeling problem
    EquilibriumProblem problem_ic(system);
    problem_ic.setTemperature(params.T, "celsius");
    problem_ic.setPressure(params.P, "bar");
    problem_ic.add("H2O", 55.51, "mol"); // H2O 55.51 M
    problem_ic.add("NaCl", 0.27, "mol"); // NaCl (aq) 0.27 M
    problem_ic.add("KCl", 0.03, "mol"); // KCl (aq)  0.03 M

//    problem_ic.add("Quartz", 10, "mol");
//    problem_ic.add("Muscovite", 10, "mol");
//    problem_ic.add("Albite", 10, "mol");
//    problem_ic.add("K-Feldspar", 10, "mol");

    // Step **: Define the boundary condition (BC)  of the reactive transport modeling problem
    EquilibriumProblem problem_bc(system);
    problem_bc.setTemperature(params.T, "celsius");
    problem_bc.setPressure(params.P, "bar");
    problem_bc.add("H2O", 55.51, "mol"); // H2O 55.51 M
    problem_bc.add("HCl", 0.1, "mol"); // HCl (aq) 0.1 M
    problem_bc.add("NaCl", 0.17, "mol"); // NaCl (aq) 0.17 M
    problem_bc.add("KCl", 0.03, "mol"); // KCl (aq)  0.03 M

    // Step **: Calculate the equilibrium states for the IC and BC
    ChemicalState state_ic = equilibrate(problem_ic);
    ChemicalState state_bc = equilibrate(problem_bc);

    // Define function to evaluate ph of the chemical system
    auto evaluate_pH = ChemicalProperty::pH(system);

    // Fetch the pH of the initial sate
    ChemicalProperties props = state_ic.properties();
    ChemicalScalar pH_ic = evaluate_pH(props);
    std::cout << "ph(IC) = " << pH_ic.val << std:: endl;
    //std::cout << "state_ic = \n" << state_ic << std:: endl;
    //
    // getchar();

    // Fetch the pH of the boundary sate
    ChemicalProperties props_bc = state_bc.properties();
    ChemicalScalar pH_bc = evaluate_pH(props_bc);
    std::cout << "ph(BC) = " << pH_bc.val << std:: endl;
    //std::cout << "state_bc = \n" << state_bc << std:: endl;
    //getchar();

    // Step **: Scale the boundary condition state
    state_bc.scaleVolume(1.0, "m3");

    // Step **: Scale the volumes of the phases in the initial condition
    state_ic.scalePhaseVolume("Aqueous", 0.1, "m3");    // 10% if the 1.0m3
    state_ic.scalePhaseVolume("Quartz", 0.3 * 0.9, "m3");    // 30% of 90% of remaining volume
    state_ic.scalePhaseVolume("Muscovite", 0.05 * 0.9, "m3");    // 5% of 90% of remaining volume
    state_ic.scalePhaseVolume("Albite", 0.33 * 0.9, "m3");    // 33% of 90% of remaining volume
    state_ic.scalePhaseVolume("K-Feldspar", 0.32 * 0.9, "m3");    // 32% of 90% of remaining volume


    std::cout << "state_ic = \n" << state_ic << std:: endl;
    getchar();
    std::cout << "state_bc = \n" << state_bc << std:: endl;
    getchar();

    // Step **: Create the mesh for the column
    Mesh mesh(params.ncells, params.xl, params.xr);

    // Step **: Create a chemical field object with every cell having state given by state_ic
    ChemicalField field(mesh.numCells(), state_ic);

    // Step **: Define the options for the reactive transport solver
    ReactiveTransportOptions reactive_transport_options;
    reactive_transport_options.use_smart_equilibrium_solver = params.use_smart_equilibrium_solver;
    reactive_transport_options.equilibrium = equilibrium_options;
    reactive_transport_options.smart_equilibrium = smart_equilibrium_options;

    // Step **: Define the reactive transport modeling
    ReactiveTransportSolver rtsolver(system);
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
    output.add("speciesMolality(Cl-)");
    output.add("speciesMolality(HCl(aq))");
    output.add("speciesMolality(Na+)");
    output.add("speciesMolality(NaCl(aq))");
    output.add("speciesMolality(OH-)");
    output.add("speciesMolality(NaOH(aq))");
    output.add("speciesMolality(NaHSiO3(aq))");
    //output.add("speciesMolality(NaAl(OH)4(aq))");
    output.add("speciesMolality(K+)");
    output.add("speciesMolality(KOH(aq))");
    output.add("speciesMolality(KCl(aq))");
    //output.add("speciesMolality(KAl(OH)4(aq))");
    output.add("speciesMolality(Al+++)");
    output.add("speciesMolality(AlOH++)");
    //output.add("speciesMolality(Al(OH)2+)");
    //output.add("speciesMolality(Al(OH)3(aq))");
    //output.add("speciesMolality(Al(OH)4-)");
    output.add("speciesMolality(Quartz)");
    output.add("speciesMolality(Diaspore)");
    output.add("speciesMolality(Gibbsite)");
    output.add("speciesMolality(Andalusite)");
    output.add("speciesMolality(Kyanite)");
    output.add("speciesMolality(Sillimanite)");
    output.add("speciesMolality(Muscovite)");
    output.add("speciesMolality(Paragonite)");
    output.add("speciesMolality(Pyrophyllite)");
    output.add("speciesMolality(Kaolinite)");
    output.add("speciesMolality(Albite)");
    output.add("speciesMolality(K-Feldspar)");

    output.filename(folder + "/" + "test.txt");

    // Step **: Create RTProfiler to track the timing and results of reactive transport
    ReactiveTransportProfiler profiler;

    // Step **: Set initial time and counter of steps in time
    double t = 0.0;
    int step = 0;

    tic(REACTIVE_TRANSPORT_STEPS);

    // Reactive transport simulations in the cycle
    while (step < params.nsteps)
    {
        // Print some progress
        std::cout << "Step " << step << " of " << params.nsteps << std::endl;

        // Perform one reactive transport time step (with profiling of some parts of the transport simulations)
        rtsolver.step(field);

        // Update the profiler after every call to step method
        profiler.update(rtsolver.result());

        // Increment time step and number of time steps
        t += params.dt;

        step += 1;
    }

    //if(params.use_smart_equilibrium_solver)
    //    rtsolver.outputClusterInfo();

    if(params.use_smart_equilibrium_solver)
        results.time_reactive_transport_smart = toc(REACTIVE_TRANSPORT_STEPS);
    else results.time_reactive_transport_conventional = toc(REACTIVE_TRANSPORT_STEPS);

    // Step **: Collect the analytics related to reactive transport performance
    auto analysis = profiler.analysis();
    auto rt_results = profiler.results();

    // Step **: Generate json output file with collected profiling data
    if(params.use_smart_equilibrium_solver)  JsonOutput(folder + "/" + "analysis-smart.json") << analysis;
    else    JsonOutput(folder + "/" + "analysis-conventional.json") << analysis;

    // Step **: Save equilibrium timing to compare the speedup of smart equilibrium solver versus conventional one
    if(params.use_smart_equilibrium_solver) {
        results.smart_equilibrium_timing = analysis.smart_equilibrium.timing;
        results.smart_equilibrium_acceptance_rate = analysis.smart_equilibrium.smart_equilibrium_estimate_acceptance_rate;

        std::cout << " smart equilibrium acceptance rate   : " << results.smart_equilibrium_acceptance_rate << " ( "
                  << results.smart_equilibrium_acceptance_rate * 100 << " % ) / "
                  << (1 - results.smart_equilibrium_acceptance_rate) * params.ncells * params.nsteps
                  << " fully evaluated GEMS out of " << params.ncells * params.nsteps
                  << " ( " << (1 - results.smart_equilibrium_acceptance_rate) * 100 << "% )" << std::endl;


    }
    else results.equilibrium_timing = analysis.equilibrium.timing;
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

    std::ostringstream reltol_stream, dt_stream;
    dt_stream << params.dt;
    reltol_stream << std::scientific << std::setprecision(1) << params.smart_equilibrium_reltol;

    std::string test_tag = "-dt-" + dt_stream.str() +
                           "-ncells-" + std::to_string(params.ncells) +
                           "-nsteps-" + std::to_string(params.nsteps) +
                           "-" + params.activity_model + "-reference";

    std::string smart_test_tag = "-dt-" + dt_stream.str() +
                                 "-ncells-" + std::to_string(params.ncells) +
                                 "-nsteps-" + std::to_string(params.nsteps) +
                                 "-reltol-" + reltol_stream.str() +
                                 "-" + params.activity_model +
                                 "-smart";

    std::string folder = "results-granite";

    folder = (params.use_smart_equilibrium_solver) ?
             folder + smart_test_tag :
             folder + test_tag;

    if (stat(folder.c_str(), &status) == -1) mkdir(folder);

    std::cout << "\nsolver                         : " << (params.use_smart_equilibrium_solver ? "smart" : "conventional") << std::endl;

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
    std::cout << "eqreltol       : " << params.smart_equilibrium_reltol << std::endl;
    std::cout << "activity model : " << params.activity_model << std::endl;

}