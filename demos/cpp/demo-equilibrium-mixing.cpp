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

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    ChemicalEditor editor;
    editor.addAqueousPhaseWithElements("H O C Ca Cl Na K Mg S Si");
    editor.addMineralPhase("Calcite");
    editor.addMineralPhase("Dolomite");

    ChemicalSystem system(editor);
    std::cout << "system = " << system << std::endl;
    getchar();

    // Amount of water
    double water_kg = 1.00;
    double T = 25.0;
    double P = 1.0;

    // Define the state with pure water
    EquilibriumInverseProblem problem_pw(system);
    problem_pw.setTemperature(T, "celsius");
    problem_pw.setPressure(P, "atm");
    problem_pw.add("H2O", water_kg, "kg");
    problem_pw.pH(7.0, "HCl", "NaOH");
           
    ChemicalState state_pw = equilibrate(problem_pw);
    std::cout << "state_pw = " << state_pw << std::endl;
    getchar();
    
    // Define the state with seawater
    EquilibriumInverseProblem problem_sw(system);
    problem_sw.setTemperature(T, "celsius");
    problem_sw.setPressure(P, "atm");
    problem_sw.add("H2O", water_kg, "kg");
    problem_sw.add("SO4--", 2712.0 * water_kg, "mg");
    problem_sw.add("Ca++", 412.3 * water_kg, "mg");
    problem_sw.add("Na+", 10768.0 * water_kg, "mg");
    problem_sw.add("K+", 399.1 * water_kg, "mg");
    problem_sw.add("Mg++", 1290 * water_kg, "mg");
    problem_sw.add("Cl-", 19353.0 * water_kg, "mg");
    problem_sw.add("HCO3-", 141.682 * water_kg, "mg");
    problem_sw.add("Si", 4.28 * water_kg, "mg");
    problem_sw.pH(8.22, "HCl", "NaOH");
    
    ChemicalState state_sw = equilibrate(problem_sw);
    std::cout << "state_sw = " << state_sw << std::endl;
    getchar();
  
    // Conduct the mixing pure water and seawater
    EquilibriumSolver solver(system);
    EquilibriumResult res = solver.solve(state_pw, T, P, state_sw.elementAmounts());

    std::cout << "state_pw (after mixing with pure water) = " << state_pw << std::endl;
    getchar();

}
