// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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
    editor.addAqueousPhaseWithElementsOf("H2O HCl CaCO3");
    editor.addMineralPhase("Calcite");

    // Define a mineral reaction
    MineralReaction reaction = editor.addMineralReaction("Calcite");
    reaction.setEquation("Calcite = Ca++ + CO3--");
    reaction.addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol");
    reaction.addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0");
    reaction.setSurfaceArea(1, "m2");

    ChemicalSystem system(editor);
    ReactionSystem reactions(editor);

    Partition partition(system);
    partition.setKineticPhases({"Calcite"});

    EquilibriumProblem problem(partition);
    problem.add("H2O", 1, "kg");
    problem.add("HCl", 1, "mmol");

    ChemicalState state = equilibrate(problem);

    state.setSpeciesMass("Calcite", 100, "g");

    state.scalePhaseVolume("Aqueous", 1, "m3");
    state.scalePhaseVolume("Calcite", 1, "m3");

    KineticPath path(reactions, partition);

    ChemicalPlot plot1 = path.plot();
    plot1.x("time(units=minute)");
    plot1.y("elementMolality(Ca units=mmolal)", "Ca");
    plot1.xlabel("Time [minute]");
    plot1.ylabel("Concentration [mmolal]");
    plot1.legend("right center");
    plot1.title("Surface area is 1 m2");

    ChemicalOutput output1 = path.output();
    output1.add("t");
    output1.add("pH");
    output1.add("speciesMass(Calcite units=g)", "Calcite");
    output1.add("speciesMolality(HCO3-)", "HCO3- [molal]");
    output1.add("speciesMolality(Ca++)", "Ca++  [molal]");
    output1.filename("path-sa-" +  std::to_string(1) + ".txt");

    double initial_calcite = state.speciesAmount("Calcite");

    // The initial chemical state conditions for both kinetic path calculations below
    ChemicalState state01 = state;
    ChemicalState state02 = state;

    // Perform the kinetic path calculation with 1 m2 surface area for 100 minutes
    path.solve(state01, 0, 30, "minutes");

    double result1_calcite = state01.speciesAmount("Calcite");

    // Change the surface area of the mineral reaction to 10 m2
    reaction.setSurfaceArea(10, "m2");

    ChemicalPlot plot2 = path.plot();
    plot2.x("time(units=minute)");
    plot2.y("elementMolality(Ca units=mmolal)", "Ca");
    plot2.xlabel("Time [minute]");
    plot2.ylabel("Concentration [mmolal]");
    plot2.legend("right center");
    plot2.title("Surface area is 10 m2");


    ChemicalOutput output2 = path.output();
    output2.add("t");
    output2.add("pH");
    output2.add("speciesMass(Calcite units=g)", "Calcite");
    output2.add("speciesMolality(HCO3-)", "HCO3- [molal]");
    output2.add("speciesMolality(Ca++)", "Ca++  [molal]");
    output2.filename("path-sa-" +  std::to_string(10) + ".txt");


    // Perform the kinetic path calculation with 10 m2 surface area for 100 minutes
    path.solve(state02, 0, 30, "minutes");

    double result2_calcite = state02.speciesAmount("Calcite");

    std::cout << "dissol. Calcite, case 1: " << initial_calcite - result1_calcite << std::endl;
    std::cout << "dissol. Calcite, case 2: " << initial_calcite - result2_calcite << std::endl;


}
