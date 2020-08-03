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

    int t0 = 0;
    int t1 = 30; // in minutes

    std::string result_filename = "kinetics-scavenging-tfinal-" + std::to_string(t1) + ".txt";

    Database database("supcrt07.xml");

    DebyeHuckelParams dhModel{};
    dhModel.setPHREEQC();

    ChemicalEditor editor(database);
    //editor.addAqueousPhaseWithElements("H Cl S O Ca Na K Mg C Si Al");
    editor.addAqueousPhaseWithElements("C Ca Cl Fe H K Mg Na O S").setChemicalModelDebyeHuckel(dhModel);
    editor.addMineralPhase("Siderite"); // FeCO3
    editor.addMineralPhase("Pyrite");   // FeS2
    editor.addMineralPhase("Hematite"); // Fe2O3

    editor.addMineralReaction("Hematite")
            .setEquation("Hematite + 4*H+ = 2*H2O(l) + 2*Fe++ + 0.5*O2(aq)")
            .addMechanism("logk = -14.60 mol/(m2*s); Ea = 66.2 kJ/mol")
            .addMechanism("logk = -9.39 mol/(m2*s); Ea = 66.2 kJ/mol; a[H+] = 1.0")
            .setSpecificSurfaceArea(10, "cm2/g");

    editor.addMineralReaction("Pyrite")
        .setEquation("Pyrite + H2O(l)  =  0.25*H+ + 0.25*SO4-- + Fe++ + 1.75*HS-")
        .addMechanism("logk = -4.55 mol/(m2*s); Ea = 56.9 kJ/mol; a[H+] = 0.500")
        .setSpecificSurfaceArea(10, "cm2/g");

    ChemicalSystem system(editor);
    ReactionSystem reactions(editor);

    Partition partition(system);
    partition.setKineticPhases(std::vector<std::string>{"Hematite", "Pyrite"});

    double T = 25.0 + 273.15;       // temperature (in units of celsius)
    double P = 1 * 1.01325 * 1e5;   // pressure (in units of atm)
    EquilibriumProblem problem(partition);
    problem.setTemperature(T);
    problem.setPressure(P);
    problem.add("H2O", 58.0, "kg");
    problem.add("Cl-", 1122.3e-3, "kg");
    problem.add("Na+", 624.08e-3, "kg");
    problem.add("SO4--", 157.18e-3, "kg");
    problem.add("Mg++", 74.820e-3, "kg");
    problem.add("Ca++", 23.838e-3, "kg");
    problem.add("K+", 23.142e-3, "kg");
    problem.add("HCO3-", 8.236e-3, "kg");
    problem.add("O2(aq)", 58e-12, "kg");
    problem.add("HS-", 0.0196504, "mol");
    problem.add("H2S(aq)", 0.167794, "mol");
    problem.add("Siderite", 0.0, "mol");
    problem.add("Pyrite", 0.0, "mol");
    problem.add("Hematite", 0.0, "mol");

    ChemicalState state0 = equilibrate(problem);
    //state0.output("shell-kinetics-benckmark-initial.txt");
    //std::cout << "state0 = " << state0 << std::endl;

    // Adding 0.5 mol of each mineral
    //state0.setSpeciesMass("Siderite", 115.86 * 0.5, "g");
    state0.setSpeciesMass("Hematite", 55.845 * 0.5, "g");

    KineticPath path(reactions, partition);

    ChemicalOutput output = path.output();
    output.filename(result_filename);
    output.add("time(units=minute)");
    output.add("pH");
    output.add("speciesMolality(H+)");
    output.add("speciesMolality(HS-)");
    output.add("speciesMolality(S2--)");
    output.add("speciesMolality(CO3--)");
    output.add("speciesMolality(HSO4-)");
    output.add("speciesMolality(H2S(aq))");
    output.add("speciesMolality(Fe++)");
    output.add("speciesMolality(Fe+++)");
    output.add("speciesMolality(Siderite)");
    output.add("speciesMolality(Pyrite)");
    output.add("speciesMolality(Hematite)");

    path.solve(state0, t0, t1, "minute");
}
