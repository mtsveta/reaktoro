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
    double smart_equilibrium_reltol = 0;
    double amount_fraction_cutoff = 0;
    double mole_fraction_cutoff = 0;

    std::string activity_model = "";
};

struct Results
{
    /// Total CPU time (in s) required by smart equilibrium scheme
    double smart_total = 0.0;

    /// Total CPU time (in s) excluding the costs for the search of the closest reference states.
    double smart_total_ideal_search = 0.0;

    /// Total CPU time (in s) required by smart equilibrium scheme
    /// excluding the costs for the search and storage of the closest reference states.
    double smart_total_ideal_search_store = 0.0;

    /// Total CPU time (in s) required by conventional equilibrium scheme
    double conventional_total = 0.0;

    /// The total time taken to perform all time steps using conventional equilibrium algorithm
    double time_reactive_transport_conventional = 0.0;

    /// The total time taken to perform all time steps using smart equilibrium algorithm
    double time_reactive_transport_smart = 0.0;

    /// The accumulated timing information of all equilibrium calculations.
    EquilibriumTiming equilibrium_timing = {};

    /// The accumulated timing information of all smart equilibrium calculations.
    SmartEquilibriumTiming smart_equilibrium_timing = {};

    // Rate of the smart equilibrium estimation w.r.t to the total chemical equilibrium calculation
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
    params.xr = 1.0; // the x-coordinates of the right boundaries
    params.ncells = 100; // the number of cells in the spacial discretization
    params.nsteps = 1000; // the number of steps in the reactive transport simulation
    params.dx = (params.xr - params.xl) / params.ncells; // the time step (in units of s)
    params.dt = 30 * minute; // the time step (in units of s)

    // Define physical and chemical parameters
    params.D = 1.0e-9;     // the diffusion coefficient (in units of m2/s)
    params.v = 1.0 / week; // the Darcy velocity (in units of m/s)
    params.T = 60.0;                     // the temperature (in units of degC)
    params.P = 100;                      // the pressure (in units of bar)

    // Define parameters of the equilibrium solvers
    params.smart_equilibrium_reltol = 0.001;

    // Define the activity model for the aqueous species
    params.activity_model = "hkf-full";
    //params.activity_model = "hkf-selected-species";
    //params.activity_model = "pitzer-full";
    //params.activity_model = "pitzer-selected-species";
    //params.activity_model = "dk-full";
    //params.activity_model = "dk-selected-species";

    // Define activity model depending on the parameter
    params.amount_fraction_cutoff = 1e-14;
    params.mole_fraction_cutoff = 1e-14;

    // Output
    outputConsole(params);

    // Results
    Results results;

    // Execute reactive transport with different solvers
    params.use_smart_equilibrium_solver = true; runReactiveTransport(params, results);
    params.use_smart_equilibrium_solver = false; runReactiveTransport(params, results);

    // Collect the time spent for total simulation (excluding search and store procedures costs)
    results.conventional_total = results.equilibrium_timing.solve;
    results.smart_total = results.smart_equilibrium_timing.solve;
    results.smart_total_ideal_search = results.smart_equilibrium_timing.solve
                                        - results.smart_equilibrium_timing.estimate_search
                                        - results.smart_equilibrium_timing.estimate_database_priority_update;
    results.smart_total_ideal_search_store = results.smart_equilibrium_timing.solve
                                                - results.smart_equilibrium_timing.estimate_search
                                                - results.smart_equilibrium_timing.estimate_database_priority_update
                                                - results.smart_equilibrium_timing.learning_storage;

    // Output speed-us
    std::cout << "speed up                            : "
              << results.conventional_total / results.smart_total << std::endl;
    std::cout << "speed up (with ideal search)        : "
              << results.conventional_total / results.smart_total_ideal_search << std::endl;
    std::cout << "speed up (with ideal search & store): "
              << results.conventional_total / results.smart_total_ideal_search_store << std::endl << std::endl;
    std::cout << " smart equilibrium acceptance rate   : " << results.smart_equilibrium_acceptance_rate << " / "
              << (1 - results.smart_equilibrium_acceptance_rate) * params.ncells *params.nsteps
              << " fully evaluated GEMS out of " << params.ncells * params.nsteps  << std::endl;
    // Output reactive transport times and speedup
    std::cout << "time_reactive_transport_conventional: " << results.time_reactive_transport_conventional << std::endl;
    std::cout << "time_reactive_transport_smart       : " << results.time_reactive_transport_smart << std::endl;
    std::cout << "reactive_transport_speedup          : " << results.time_reactive_transport_conventional / results.time_reactive_transport_smart << std::endl;
    // Output total time
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
    smart_equilibrium_options.amount_fraction_cutoff = params.amount_fraction_cutoff;
    smart_equilibrium_options.mole_fraction_cutoff = params.mole_fraction_cutoff;

    // Step **: Construct the chemical system with its phases and species (using ChemicalEditor)
    ChemicalEditor editor;

    // Depending on the activity model, define it using ChemicalEditor
    if(params.activity_model == "hkf-full"){
        // HKF full system
        editor.addAqueousPhaseWithElements("H O Na Cl Ca Mg C");
    }
    else if(params.activity_model == "hkf-selected-species"){
        // HKF selected species
        editor.addAqueousPhase("H2O(l) H+ OH- Na+ Cl- Ca++ Mg++ HCO3- CO2(aq) CO3-- CaCl+ Ca(HCO3)+ MgCl+ Mg(HCO3)+");
    }
    else if(params.activity_model == "pitzer-full"){
        // Pitzer full system
        editor.addAqueousPhaseWithElements("H O Na Cl Ca Mg C")
                .setChemicalModelPitzerHMW()
                .setActivityModelDrummondCO2();
    }
    else if(params.activity_model == "pitzer-selected-species"){
        // Pitzer selected species
        editor.addAqueousPhase("H2O(l) H+ OH- Na+ Cl- Ca++ Mg++ HCO3- CO2(aq) CO3-- CaCl+ Ca(HCO3)+ MgCl+ Mg(HCO3)+")
                .setChemicalModelPitzerHMW()
                .setActivityModelDrummondCO2();
    }
    else if(params.activity_model == "dk-full"){
        // Debye-Huckel full system
        editor.addAqueousPhaseWithElements("H O Na Cl Ca Mg C")
                .setChemicalModelDebyeHuckel()
                .setActivityModelDrummondCO2();
    }
    else if(params.activity_model == "dk-selected-species"){
        // Debye-Huckel selected species
        editor.addAqueousPhase("H2O(l) H+ OH- Na+ Cl- Ca++ Mg++ HCO3- CO2(aq) CO3-- CaCl+ Ca(HCO3)+ MgCl+ Mg(HCO3)+")
                .setChemicalModelDebyeHuckel()
                .setActivityModelDrummondCO2();
    }
    editor.addMineralPhase("Quartz");
    editor.addMineralPhase("Calcite");
    editor.addMineralPhase("Dolomite");

    // Step **: Create the ChemicalSystem object using the configured editor
    ChemicalSystem system(editor);

    // Step **: Create the Partition of inert and equilibrium species
    Partition partition(system);
    partition.setInertSpecies({"Quartz"});

    // Step **: Define the initial condition (IC) of the reactive transport modeling problem
    EquilibriumProblem problem_ic(system);
    problem_ic.setTemperature(params.T, "celsius");
    problem_ic.setPressure(params.P, "bar");
    problem_ic.add("H2O",   1.0, "kg");
    problem_ic.add("O2",    1.0, "umol");
    problem_ic.add("NaCl",  0.7, "mol");
    problem_ic.add("CaCO3", 10,  "mol");
    problem_ic.add("SiO2",  10,  "mol");
    problem_ic.add("MgCl2", 1e-10, "mol");

    // Step **: Define the boundary condition (BC)  of the reactive transport modeling problem
    EquilibriumProblem problem_bc(system);
    problem_bc.setTemperature(params.T, "celsius");
    problem_bc.setPressure(params.P, "bar");
    problem_bc.add("H2O",   1.00, "kg");
    problem_bc.add("O2",    1.0, "umol");
    problem_bc.add("NaCl",  0.90, "mol");
    problem_bc.add("MgCl2", 0.05, "mol");
    problem_bc.add("CaCl2", 0.01, "mol");
    problem_bc.add("CO2",   0.75, "mol");

    // Step **: Calculate the equilibrium states for the IC and BC
    ChemicalState state_ic = equilibrate(problem_ic);
    ChemicalState state_bc = equilibrate(problem_bc);

    // Step **: Scale the boundary condition state
    state_bc.scaleVolume(1.0, "m3");

    // Step **: Scale the volumes of the phases in the initial condition
    state_ic.scalePhaseVolume("Aqueous", 0.1, "m3");    // 10% if the 1.0m3
    state_ic.scalePhaseVolume("Quartz", 0.882, "m3");   // 0.882 = 0.98 * 0.9 (0.9 is 90% of 1.0m3, 0.98 is 98% quartz of the rock)
    state_ic.scalePhaseVolume("Calcite", 0.018, "m3");  // 0.018 = 0.02 * 0.9 (0.9 is 90% of 1.0m3, 0.02 is 2% calcite of the rock)

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
    ReactiveTransportSolver rtsolver(partition);
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
    output.add("speciesAmount(Dolomite)");
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
        // Print the progress of calculations
        std::cout << "Step " << step << " of " << params.nsteps << std::endl;

        // Perform one reactive transport time step (with profiling of some parts of the transport simulations)
        rtsolver.step(field);

        // Update the profiler after every call to step method
        profiler.update(rtsolver.result());

        // Increment time step and number of time steps
        t += params.dt;

        step += 1;
    }

    if(params.use_smart_equilibrium_solver)
    {
        rtsolver.outputClusterInfo();
        results.time_reactive_transport_smart = toc(REACTIVE_TRANSPORT_STEPS);
    }
        
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

        std::cout << "smart equilibrium acceptance rate   : " << results.smart_equilibrium_acceptance_rate << " / "
                  << (1 - results.smart_equilibrium_acceptance_rate) * params.ncells *params.nsteps
                  << " fully evaluated GEMS out of " << params.ncells * params.nsteps  << std::endl;

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

    std::string folder = "results-dolomitization";
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

