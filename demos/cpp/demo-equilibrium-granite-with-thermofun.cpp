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

#include <ThermoFun/ThermoFun.h>

int main()
{


    ThermoFun::Database database("databases/thermofun/aq17-thermofun.json");

    ChemicalEditor editor(database);
    editor.setTemperatures({300}, "celsius");
    editor.setPressures({85.88}, "bar");

//    StringList selected_species = "H2O@ H+ OH- Cl- HCl@ Na+ NaOH@ NaHSiO3@ NaCl@ NaAl(OH)4@ "
//                                  "K+ KOH@ KCl@ KAl(OH)4@ Al+++ AlOH++ Al(OH)2+ Al(OH)3@ Al(OH)4-";
    // Missing:KAl(OH)4@
    StringList selected_species = "H2O@ H+ OH- Cl- HCl@ Na+ NaOH@ NaHSiO3@ NaCl@ NaAl(OH)4@ "
                                  "K+ KOH@ KCl@ Al+3 AlOH+2 Al(OH)2+ Al(OH)3@ Al(OH)4-";

    // Missing: NaAl(OH)4(aq), KAl(OH)4(aq), Al(OH)2+, Al(OH)3(aq),  Al(OH)4-
    //editor.addAqueousPhaseWithElements(selected_elements); // total 40 species
    editor.addAqueousPhase(selected_species); // total 25 species
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
    //editor.addMineralPhase("K-Feldspar"); // K(AlSi3)O8 // Can be replaced by Sanidine or Microcline?
    editor.addMineralPhase("Microcline");

    ChemicalSystem system(editor);
    std::cout << system << std::endl;

    EquilibriumProblem problem(system);
    problem.setTemperature(300, "celsius");
    problem.setPressure(84.75, "atm"); // is Psat a water saturation at the T = 300? 8588 *1e3 Pa / 85.88 bar / 84.75 atm
    // Residual fluid
    problem.add("H2O", 55.51, "mol"); // H2O 55.51 M
    problem.add("NaCl", 0.27, "mol"); // NaCl (aq) 0.27 M
    problem.add("KCl", 0.03, "mol"); // KCl (aq)  0.03 M
    // Initial amounts of minerals?

    ChemicalState state = equilibrate(problem);

    state.scalePhaseVolume("Aqueous", 0.1, "m3");    // 10% if the 1.0m3
    state.scalePhaseVolume("Quartz", 0.3 * 0.9, "m3");    // 30% of 90% of remaining volume
    state.scalePhaseVolume("Muscovite", 0.05 * 0.9, "m3");    // 5% of 90% of remaining volume
    state.scalePhaseVolume("Albite", 0.33 * 0.9, "m3");    // 33% of 90% of remaining volume
    //state.scalePhaseVolume("K-Feldspar", 0.32 * 0.9, "m3");    // 32% of 90% of remaining volume
    state.scalePhaseVolume("Sanidine", 0.32 * 0.9, "m3");    // 32% of 90% of remaining volume

    std::cout << state << std::endl;
}
