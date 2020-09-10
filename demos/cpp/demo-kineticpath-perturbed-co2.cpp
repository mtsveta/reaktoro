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

    ChemicalEditor editor;
    //editor.addAqueousPhaseWithElementsOf("H2O CO2 CaCO3");
    editor.addAqueousPhaseWithElementsOf("H2O HCl CaCO3");
    editor.addMineralPhase("Calcite");

    MineralReaction reaction = editor.addMineralReaction("Calcite");
    reaction.setEquation("Calcite = Ca++ + CO3--");
    reaction.addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol");
    reaction.addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0");
    reaction.setSpecificSurfaceArea(10, "cm2/g");

    ChemicalSystem system(editor);
    std::cout << system << std::endl;

    ReactionSystem reactions(editor);

    Partition partition(system);
    partition.setKineticSpecies({"Calcite"});

    EquilibriumProblem problem(partition);
    problem.setTemperature(60, "celsius");
    problem.setPressure(100, "bar");
    problem.add("H2O", 1, "kg");

    // Create a vector of chemical states
    std::vector<ChemicalState> chemical_states;
    const unsigned int states_num = 5;
    for(unsigned i = 0; i < states_num; ++i){
        // Add CO2
        problem.add("CO2", 0.1, "mol");
        // Create initial state by equilibration
        chemical_states.emplace_back(ChemicalState(equilibrate(problem)));
        // Set amount of Calcite
        chemical_states[i].setSpeciesMass("Calcite", 100, "g");
    }

    KineticPath path(reactions, partition);

    ChemicalOutput output = path.output();
    output.add("t");
    output.add("pH");
    output.add("speciesMass(Calcite units=g)", "Calcite");
    output.add("speciesMolality(HCO3-)", "HCO3- [molal]");
    output.add("speciesMolality(Ca++)", "Ca++  [molal]");

    for(unsigned i = 0; i < states_num; ++i){

        std::cout << "**********************************************************************************" << std::endl;
        std::cout << i + 1 << std::endl;
        std::cout << "**********************************************************************************" << std::endl;

        output.filename(folder + "/path" + std::to_string(i + 1));
        double initial_calcite = chemical_states[i].speciesAmount("Calcite");
        std::cout << "initial Calcite: " << initial_calcite << std::endl;
        std::cout << "initial CO2    : " << chemical_states[i].speciesAmount("CO2(aq)") << std::endl;
        std::cout << "initial pH     : " << ChemicalProperty::pH(chemical_states[i].system())(chemical_states[i].properties()) << std::endl;
        path.solve(chemical_states[i], params.t0, params.tfinal);
        double result_calcite = chemical_states[i].speciesAmount("Calcite");
        double dissolved_calcite = initial_calcite - result_calcite;
        std::cout << "dissol. Calcite: " << dissolved_calcite << std::endl;

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