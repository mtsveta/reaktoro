//
// Created by skyas on 16.02.20.
//

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
    int ncells; // the number of cells in the spacial discretization
    int nsteps; // the number of steps in the reactive transport simulation
    double xl; // the x-coordinates of the left boundaries
    double xr; // the x-coordinates of the right boundaries
    double dx; // the space step (in units of m)
    double dt; // the time step (in units of s)

    // Physical params
    double D; // the diffusion coefficient (in units of m2/s)
    double v; // the Darcy velocity (in units of m/s)
    double T; // the temperature (in units of degC)
    double P; // the pressure (in units of bar)

    // Solver params
    bool use_smart_eqilibirum_solver;
    bool track_statistics;
    double smart_equlibrium_reltol;
    double smart_equlibrium_abstol;

    double tol;

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
    int second = 1;
    int minute = 60;
    int hour = 60 * minute;
    int day = 24 * hour;
    int week = 7 * day;
    int month = 30 * day;
    int year = 365 * day;

    // Step 2: Define parameters for the reactive transport simulation
    Params params;

    // Define discretization parameters
    params.xl = 0.0; // the x-coordinates of the left boundaries
    /*
    params.xr = 0.1; // the x-coordinates of the right boundaries
    params.ncells = 10; // the number of cells in the spacial discretization
    */
    ///*
    params.xr = 100.0; // the x-coordinates of the right boundaries
    params.ncells = 100; // the number of cells in the spacial discretization
    //*/
    params.nsteps = 10000; // the number of steps in the reactive transport simulation
    params.dx = (params.xr - params.xl) / params.ncells; // the time step (in units of s)
    params.dt = 0.0416*day; // the time step (in units of s)

    // Define physical and chemical parameters
    params.D = 0.0;     // the diffusion coefficient (in units of m2/s)
    params.v = 1.05e-5; // the Darcy velocity (in units of m/s)
    params.T = 25.0;                     // the temperature (in units of degC)
    params.P = 1.01325;                      // the pressure (in units of bar)

    // Define parameters of the equilibrium solvers
    params.smart_equlibrium_reltol = 1e-1;
    params.smart_equlibrium_abstol = 1e-8;
    params.tol = 8e-1;
    params.track_statistics = true;

    // Output
    outputConsole(params);

    // Results
    Results results;

    // Execute reactive transport with different solvers
    params.use_smart_eqilibirum_solver = true; runReactiveTransport(params, results);
    params.use_smart_eqilibirum_solver = false; runReactiveTransport(params, results);

    results.conventional_total = results.equilibrium_timing.solve;
    results.smart_total = results.smart_equilibrium_timing.solve;
    results.smart_total_ideal_search = results.smart_equilibrium_timing.solve
                                       - results.smart_equilibrium_timing.estimate_search;
    results.smart_total_ideal_search_store = results.smart_equilibrium_timing.solve
                                             - results.smart_equilibrium_timing.estimate_search
                                             - results.smart_equilibrium_timing.learn_storage;

    // Output speed-us
    std::cout << std::defaultfloat << "speed up                            : "
              << results.conventional_total / results.smart_total << std::endl;
    std::cout << "speed up (with ideal search)        : "
              << results.conventional_total / results.smart_total_ideal_search << std::endl;
    std::cout << "speed up (with ideal search & store): "
              << results.conventional_total / results.smart_total_ideal_search_store << std::endl << std::endl;
    std::cout << " smart equilibrium acceptance rate   : " << results.smart_equilibrium_acceptance_rate << " / "
              << (1 - results.smart_equilibrium_acceptance_rate) * params.ncells *params.nsteps
              << " fully evaluated GEMS out of " << params.ncells * params.nsteps  << std::endl;

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
    smart_equilibrium_options.reltol = params.smart_equlibrium_reltol;
    smart_equilibrium_options.abstol = params.smart_equlibrium_abstol;
    smart_equilibrium_options.tol = params.tol;

    // Step **: Construct the chemical system with its phases and species (using ChemicalEditor)
    Database database("supcrt07.xml");

    DebyeHuckelParams dhModel{};
    dhModel.setPHREEQC();

    ChemicalEditor editor(database);
    // Default chemical model (HKF extended Debye-Hückel model)
    // editor.addAqueousPhase("H2O(l) H+ OH- Na+ Cl- Ca++ Mg++ HCO3- CO2(aq) CO3--");
    // Create aqueous phase with all possible elements
    // Set a chemical model of the phase with the Pitzer equation of state
    // With an exception for the CO2, for which Drummond model is set
    editor.addAqueousPhase({"H2O(l)",  "H+", "OH-",
                            "HCO3-", "Mg(HCO3)+", "Ca(HCO3)+", "MgCO3(aq)",  "CO3--", "CaCO3(aq)" ,
                            "Ca++", "CaSO4(aq)", "CaOH+",
                            "Cl-", "FeCl++", "FeCl2(aq)", "FeCl+",
                            "Fe++", "FeOH+",  "FeOH++", "Fe+++",
                            "H2(aq)",
                            "K+", "KSO4-",
                            "Mg++", "MgSO4(aq)", "MgCO3(aq)", "MgOH+",
                            "Na+", "NaSO4-",
                            "O2(aq)",
                            "H2S(aq)", "HS-", "S5--", "S4--", "S3--", "S2--",
                            "SO4--", "NaSO4-", "MgSO4(aq)", "CaSO4(aq)", "KSO4-", "HSO4-"}).setChemicalModelDebyeHuckel(dhModel);
    editor.addMineralPhase("Pyrrhotite");
    editor.addMineralPhase("Siderite");

    // Step **: Create the ChemicalSystem object using the configured editor
    ChemicalSystem system(editor);
    //if (params.use_smart_eqilibirum_solver) std::cout << "system = \n" << system << std:: endl;

    // Step **: Define the initial condition (IC) of the reactive transport modeling problem
    EquilibriumInverseProblem problem_ic(system);
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
    problem_ic.add("Pyrrhotite", 0.0, "mol");
    problem_ic.add("Siderite", 0.5, "mol");
    problem_ic.pH(8.951);
    problem_ic.pE(8.676);

    // Step **: Define the boundary condition (BC)  of the reactive transport modeling problem
    EquilibriumInverseProblem problem_bc(system);
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
    problem_bc.add("Pyrrhotite", 0.0, "mol");
    problem_bc.add("Siderite", 0.0, "mol");
    problem_bc.add("HS-", 0.0196504, "mol");
    problem_bc.add("H2S(aq)", 0.167794, "mol");
    problem_bc.pH(5.726);
    problem_bc.pE(8.220);

    // Step **: Calculate the equilibrium states for the IC and BC
    ChemicalState state_ic = equilibrate(problem_ic);
    ChemicalState state_bc = equilibrate(problem_bc);

    // Step **: Scale the boundary condition state
    state_bc.scaleVolume(1.0, "m3");

    // Step **: Scale the volumes of the phases in the initial condition
    state_ic.scalePhaseVolume("Aqueous", 0.1, "m3");    // 10% if the 1.0m3
    state_ic.scaleVolume(1.0, "m3");
    //state_ic.scalePhaseVolume("Quartz", 0.882, "m3");   // 0.882 = 0.98 * 0.9 (0.9 is 90% of 1.0m3, 0.98 is 98% quartz of the rock)
    //state_ic.scalePhaseVolume("Calcite", 0.018, "m3");  // 0.018 = 0.02 * 0.9 (0.9 is 90% of 1.0m3, 0.02 is 2% calcite of the rock)

    // Step **: Create the mesh for the column
    Mesh mesh(params.ncells, params.xl, params.xr);

    // Step **: Create a chemical field object with every cell having state given by state_ic
    ChemicalField field(mesh.numCells(), state_ic);

    // Step **: Define the options for the reactive transport solver
    ReactiveTransportOptions reactive_transport_options;
    reactive_transport_options.use_smart_equilibrium_solver = params.use_smart_eqilibirum_solver;
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
    output.add("speciesMolality(HS-)");
    output.add("speciesMolality(S2--)");
    output.add("speciesMolality(SO4--)");
    output.add("speciesMolality(HSO4-)");
    output.add("speciesMolality(H2S(aq))");
    output.add("phaseAmount(Pyrrhotite)");
    output.add("phaseAmount(Siderite)");
    output.filename(folder + "/" + "state.txt");

    // Step **: Create RTProfiler to track the timing and results of reactive transport
    ReactiveTransportProfiler profiler;

    // Step **: Set initial time and counter of steps in time
    double t = 0.0;
    int step = 0;

    tic(0);

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

    toc(0, (params.use_smart_eqilibirum_solver ? results.time_reactive_transport_smart : results.time_reactive_transport_conventional) );

    // Step **: Collect the analytics related to reactive transport performance
    auto analysis = profiler.analysis();
    auto rt_results = profiler.results();

    // Step **: Generate json output file with collected profiling data
    if(params.use_smart_eqilibirum_solver)  JsonOutput(folder + "/" + "analysis-smart.json") << analysis;
    else    JsonOutput(folder + "/" + "analysis-conventional.json") << analysis;

    // Step **: Save equilibrium timing to compare the speedup of smart equilibrium solver versus conventional one
    if(params.use_smart_eqilibirum_solver) {
        results.smart_equilibrium_timing = analysis.smart_equilibrium.timing;
        results.smart_equilibrium_acceptance_rate = analysis.smart_equilibrium.smart_equilibrium_estimate_acceptance_rate;
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

    std::ostringstream tol_stream, dt_stream;
    dt_stream << params.dt;
    tol_stream << std::scientific << std::setprecision(1) << params.tol;

    std::string test_tag = "-dt-" + dt_stream.str() +
                           "-ncells-" + std::to_string(params.ncells) +
                           "-nsteps-" + std::to_string(params.nsteps) + "-reference";

    std::string smart_test_tag = "-dt-" + dt_stream.str() +
                                 "-ncells-" + std::to_string(params.ncells) +
                                 "-nsteps-" + std::to_string(params.nsteps) +
                                 "-tol-" + tol_stream.str() + "-smart";      // name of the folder with results
    std::string folder = "results-pitzer-geometric-nnsearch";
    folder = (params.use_smart_eqilibirum_solver) ?
             folder + smart_test_tag :
             folder + test_tag;
    //std::string folder = "results-pitzer-" + test_tag;
    if (stat(folder.c_str(), &status) == -1) mkdir(folder.c_str());

    std::cout << "\nsolver                         : " << (params.use_smart_eqilibirum_solver == true ? "smart" : "conventional") << std::endl;

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
    std::cout << "eqabstol  : " << params.smart_equlibrium_abstol << std::endl;
    std::cout << "eqreltol  : " << params.smart_equlibrium_reltol << std::endl;
    std::cout << "tol       : " << params.tol << std::endl;

}