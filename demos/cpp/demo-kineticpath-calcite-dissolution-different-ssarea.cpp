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

#if defined _WIN32      // for creating a new folder
#include <windows.h>
#ifdef __MINGW32__
#include <sys/stat.h>
#endif
#else
#include <sys/stat.h>
#endif

#include <algorithm>

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

struct Params{

    // Discretisation params
    double t0;      // starting time of simulations
    double tfinal;  // final time of simulation

};
/// Forward declaration
auto mkdir(const std::string &) -> bool;
auto makeResultsFolder(const Params &) -> std::string;
auto runKinetics(const Params &) -> void;

int main()
{
    int minute(60);

    Params params = {};

    params.t0 = 0;
    params.tfinal = 30 * minute;

    runKinetics(params);
}

auto runKinetics(const Params & params) -> void{

    auto folder = makeResultsFolder(params);
    std::vector<double> surface_areas = {1, 10, 100, 1000, 10000};
    std::vector<ChemicalState> chemical_states;
    const unsigned int states_num = surface_areas.size();

    ChemicalEditor editor;
    editor.addAqueousPhaseWithElements("H O Na Cl Ca Mg C");
    editor.addMineralPhase("Calcite");

    MineralReaction reaction = editor.addMineralReaction("Calcite");
    reaction.setEquation("Calcite = Ca++ + CO3--");
    reaction.addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol");
    reaction.addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0");

    for(unsigned i = 0; i < states_num; ++i)
    {
        reaction.setSpecificSurfaceArea(surface_areas[i], "cm2/g");

        ChemicalSystem system(editor);

        ReactionSystem reactions(editor);

        Partition partition(system);
        partition.setKineticSpecies({"Calcite"});

        EquilibriumProblem problem(partition);
        problem.setTemperature(60, "celsius");
        problem.setPressure(100, "bar");
        problem.add("H2O",   1.0, "kg");
        problem.add("O2",    1.0, "umol");
        problem.add("NaCl",  0.9, "mol");
        problem.add("MgCl2", 0.05, "mol");
        problem.add("CaCl2", 0.01, "mol");
        problem.add("CO2",   0.75, "mol");

        // Create initial state by equilibration
        chemical_states.emplace_back(ChemicalState(equilibrate(problem)));
        // Set amount of Calcite
        chemical_states[i].setSpeciesMass("Calcite", 100, "g");

        chemical_states[i].scalePhaseVolume("Aqueous", 1, "m3");
        chemical_states[i].scalePhaseVolume("Calcite", 1, "m3");

        KineticPath path(reactions, partition);

        ChemicalOutput output = path.output();
        output.add("t");
        output.add("pH");
        output.add("speciesMass(Calcite units=g)", "Calcite");
        output.add("speciesMolality(HCO3-)", "HCO3- [molal]");
        output.add("speciesMolality(Ca++)", "Ca++  [molal]");

        std::cout << "**********************************************************************************" << std::endl;
        std::cout << "Surface area: " << surface_areas[i] << " cm2/g" << std::endl;

        output.filename(folder + "/path-sa-" +  std::to_string(surface_areas[i]) + ".txt");

        double initial_calcite = chemical_states[i].speciesAmount("Calcite");
        //std::cout << "initial Calcite: " << initial_calcite << std::endl;
        //std::cout << "initial CO2    : " << chemical_states[i].speciesAmount("CO2(aq)") << std::endl;
        // std::cout << "initial pH     : " << ChemicalProperty::pH(chemical_states[i].system())(chemical_states[i].properties()) << std::endl;
        path.solve(chemical_states[i], params.t0, params.tfinal);
        double result_calcite = chemical_states[i].speciesAmount("Calcite");
        double dissolved_calcite = initial_calcite - result_calcite;
        std::cout << "Dissol. Calcite: " << dissolved_calcite << std::endl;

    }

}

/// Make directory for Windows and Linux
auto mkdir(const std::string & folder) -> bool
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
auto makeResultsFolder(const Params & params) -> std::string
{
    struct stat status = {0};               // structure to get the file status

    std::ostringstream tol_stream, t0_stream, tfinal_stream;
    t0_stream << params.t0;
    tfinal_stream << params.tfinal;

    std::string test_tag = "-t0-" + t0_stream.str() +
                           "-tfinal-" + tfinal_stream.str();      // name of the folder with results
    std::string folder = "../kinetics-perturbed-co2" + test_tag;
    if (stat(folder.c_str(), &status) == -1) mkdir(folder);

    return folder;
}

/*

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{

    std::string activity_model = "hkf-full";

    // Step **: Construct the chemical system with its phases and species (using ChemicalEditor)
    ChemicalEditor editor;

    // Step **: Add aqueous phase, default chemical model (HKF extended Debye-Hückel model)
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
    // Step **: Add mineral phase
    editor.addMineralPhase("Dolomite");
    editor.addMineralPhase("Calcite");
    editor.addMineralPhase("Quartz");

    MineralReaction reaction = editor.addMineralReaction("Calcite");
    reaction.setEquation("Calcite = Ca++ + CO3--");
    reaction.addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol"); // neutral
    reaction.addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0"); // acidic
    reaction.setSpecificSurfaceArea(1, "cm2/g");

    // Step **: Create the ChemicalSystem object using the configured editor
    ChemicalSystem system(editor);

    // Step **: Create the ReactionSystem instances
    ReactionSystem reactions(editor);

    // Step **: Create the ReactionSystem instances
    Partition partition(system);
    partition.setInertSpecies({"Quartz"});
    partition.setKineticSpecies({"Calcite"});

    // Step **: Define the initial condition (IC) of the reactive transport modeling problem
    EquilibriumProblem problem_ic(partition);
    problem_ic.setTemperature(params.T, "celsius");
    problem_ic.setPressure(params.P, "bar");
    problem_ic.add("H2O",   1.0, "kg");
    problem_ic.add("O2",    1.0, "umol");
    problem_ic.add("NaCl",  0.7, "mol");
    problem_ic.add("CaCO3", 10,  "mol");
    problem_ic.add("SiO2",  10,  "mol");
    problem_ic.add("MgCl2", 1e-10, "mol");

    std::vector<double> surface_areas = {1, 10, 100, 1000, 10000};

    for(auto area : surface_areas)
    {
        problem_ic.add("NaCl",  0.90, "mol");
        problem_ic.add("MgCl2", 0.05, "mol");
        problem_ic.add("CaCl2", 0.01, "mol");
        problem_ic.add("CO2",   0.75, "mol");

        // Step **: Calculate the equilibrium states for the IC and BC
        ChemicalState state_ic = equilibrate(problem_ic);//

        // Step **: Scale the volumes of the phases in the initial condition
        state_ic.scalePhaseVolume("Aqueous", 0.1, "m3");    // 10% if the 1.0m3
        state_ic.scalePhaseVolume("Quartz", 0.81, "m3");   // 0.81 = 0.90 * 0.9 (0.9 of rock is due to 10% porosity, 0.90 is 90% quartz of the rock)
        state_ic.scalePhaseVolume("Calcite", 0.054, "m3");  // 0.054 = 0.06 * 0.9 (0.9 of rock due to 10% porosity, 0.06 is 6% calcite of the rock)
        state_ic.scalePhaseVolume("Dolomite", 0.036, "m3");  // 0.036 = 0.04 * 0.9 (0.9 of rock due to 10% porosity, 0.04 is 4% calcite of the rock)


        // Change the surface area of the mineral reaction
        reaction.setSurfaceArea(area, "cm2/g");

        KineticPath path(reactions, partition);

        ChemicalPlot plot = path.plot();
        plot.x("time(units=minute)");
        plot.y("phaseMass(Calcite units=g)", "Calcite");
        plot.xlabel("Time [minute]");
        plot.ylabel("Mass [g]");
        plot.legend("right center");
        plot.title("Surface area is" + std::to_string(area) + "cm2/g");

        ChemicalOutput output = path.output();
        output.add("time(units=minute)");
        output.add("phaseMass(Calcite units=g)", "Calcite");
        output.filename("path-sa-" +  std::to_string(area) + ".txtcmake .")

        // The initial chemical state conditions for both kinetic path calculations below
        ChemicalState state = state_ic;

        // Perform the kinetic path calculation with 1 m2 surface area for 100 minutes
        path.solve(state, 0, 30, "minutes");
    }
}


*/
