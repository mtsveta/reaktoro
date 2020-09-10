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

#include <algorithm>

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    double t0 = 0;
    double dt = 0.1; // in minutes
    int n = 25;

    //std::string activity_model = "pitzer";
    std::string activity_model = "hkf";

    Database database("supcrt07.xml");

    DebyeHuckelParams dhModel{};
    dhModel.setPHREEQC();

    ChemicalEditor editor(database);
    if(activity_model=="pitzer")
        editor.addAqueousPhaseWithElements("H Cl S O Ca Na K Mg C Si Al")
                .setChemicalModelDebyeHuckel(dhModel)
                .setChemicalModelPitzerHMW()
                .setActivityModelDrummondCO2();
    else if(activity_model=="hkf")
        editor.addAqueousPhaseWithElements("H Cl S O Ca Na K Mg C Si Al")
                .setChemicalModelDebyeHuckel(dhModel);

    std::vector<std::string> minerals = {"Halite", "Calcite", "Dolomite", "K-Feldspar", "Quartz", "Kaolinite"};
    auto num_minerals = minerals.size();

    for(const auto& mineral : minerals)
        editor.addMineralPhase(mineral);

    editor.addMineralReaction("Halite")
            .setEquation("Halite = Na+ + Cl-")
            .addMechanism("logk = -0.21 mol/(m2*s); Ea = 7.4 kJ/mol")
            .setSpecificSurfaceArea(10, "cm2/g");

    editor.addMineralReaction("Calcite")
        .setEquation("Calcite = Ca++ + CO3--")
        .addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol")
        .addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0")
        .setSpecificSurfaceArea(10, "cm2/g");

    editor.addMineralReaction("Dolomite")
        .setEquation("Dolomite = Ca++ + Mg++ + 2*CO3--")
        .addMechanism("logk = -7.53 mol/(m2*s); Ea = 52.2 kJ/mol")
        .addMechanism("logk = -3.19 mol/(m2*s); Ea = 36.1 kJ/mol; a[H+] = 0.5")
        .setSpecificSurfaceArea(10, "cm2/g");

    // Quartz + 2*H2O = H4SiO4(aq)
    editor.addMineralReaction("Quartz")
        .setEquation("Quartz + H2O(l) = H+ + HSiO3-")
        .addMechanism("logk = -13.99 mol/(m2*s); Ea = 87.7 kJ/mol")
        .setSpecificSurfaceArea(10, "cm2/g");

    // K-Feldspar: K(AlSi3)O8 # potassium feldspar, orthoclase
    editor.addMineralReaction("K-Feldspar")
        .setEquation("K-Feldspar + H2O(l) + H+ = Al+++ + 3*HSiO3- + K+")
        .addMechanism("logk = -12.41 mol/(m2*s); Ea = 38.0 kJ/mol")
        .addMechanism("logk = -10.06 mol/(m2*s); Ea = 51.7 kJ/mol; a[H+] = 0.5")
        .setSpecificSurfaceArea(10, "cm2/g");

    // Kaolinite: Al2Si2O5(OH)4
    editor.addMineralReaction("Kaolinite")
        .setEquation("Kaolinite + 4*H+ = 3*H2O(l) + 2*HSiO3- + 2*Al+++")
        .addMechanism("logk = -13.18 mol/(m2*s); Ea = 22.2 kJ/mol")
        .addMechanism("logk = -11.31 mol/(m2*s); Ea = 65.9 kJ/mol; a[H+] = 0.777")
        .setSpecificSurfaceArea(10, "cm2/g");

    ChemicalSystem system(editor);
    ReactionSystem reactions(editor);

    std::vector<Partition> partitions(num_minerals + 1);
    std::fill(partitions.begin(), partitions.end(), Partition(system));

    for(int i = 0; i < num_minerals + 1; i++){
        auto kinetic_minerals = std::vector<std::string>(minerals.begin(), minerals.end() - i);
        std::cout << "Kinetic minerals: "; for(const auto& mineral : kinetic_minerals) std::cout << mineral << ", "; std::cout << std::endl;

        // Initialize current partition with kinetic minerals
        partitions[i].setKineticPhases(kinetic_minerals);
        auto num_kinetic_minerals = kinetic_minerals.size();

        double T = 61.0 + 273.15;       // temperature (in units of celsius)
        double P = 1 * 1.01325 * 1e5;   // pressure (in units of atm)
        EquilibriumInverseProblem problem(partitions[i]);
        problem.setTemperature(T);
        problem.setPressure(P);
        problem.add("H2O", 1, "kg");
        problem.add("Na", 6.27 , "mol"); // 6.27 mol / kgw
        problem.add("Cl", 6.27, "mol");  // 6.27 mol / kgw
        problem.add("Mg", 1e-4, "mol");  // 0.0001 mol / kgw
        problem.add("Ca", 1e-4, "mol");  // 0.0001 mol / kgw
        problem.add("C", 1e-4, "mol");  // 0.0001 mol / kgw
        problem.add("Si", 1e-4, "mol");  // 0.0001 mol / kgw
        problem.add("Al", 1e-4, "mol");  // 6.27 mol / kgw
        problem.add("K", 1e-9, "mol");
        problem.pH(7.0);

        ChemicalState state0 = equilibrate(problem);

        state0.setSpeciesMass("Calcite", 100.0869 * 10, "g"); //  molar mass of CaCO3 = 100.0869 g/mol
        state0.setSpeciesMass("Quartz", 60.08 * 10, "g"); // molar mass of SiO2 = 60.08 g/mol
        state0.setSpeciesMass("Halite", 58.44 * 10, "g"); // molar mass of NaCl = 58.44 g/mol
        state0.setSpeciesMass("Dolomite", 184.4008 * 10, "g"); // molar mass of CaMg(CO3)2 = 184.4008 g/mol
        state0.setSpeciesMass("Kaolinite", 258.1604 * 10, "g"); // molar mass of Al2Si2O5(OH)4 =  258.1604 g/mol
        state0.setSpeciesMass("K-Feldspar", 278.3315 * 10, "g"); // molar mass of K(AlSi3)O8 = 278.3315 g/mol

        double xl = 0.0;
        double xr = 25.0;
        Index ncells = 25;

        // Step **: Create the mesh for the column
        Mesh mesh(ncells, xl, xr);
        // Step **: Create a chemical field object with every cell having state given by state_ic
        ChemicalField field(mesh.numCells(), state0);

        KineticPath path(reactions, partitions[i]);

        std::string result_filename = "kinetics-benchmark-kin-" + std::to_string(num_kinetic_minerals)
                                                           + "-eq-" + std::to_string(num_minerals - num_kinetic_minerals)
                                                           + "-dt-" + std::to_string(dt)
                                                           + "-activity-model-" + activity_model +
                                                           + "-uniform.txt";
        ChemicalOutput output = path.output();
        output.filename(result_filename);
        output.add("time(units=second)");
        output.add("speciesMolality(Na+ units=mmolal)", "Na+ [mmol]");
        output.add("speciesMolality(Cl- units=mmolal)", "Cl- [mmol]");
        output.add("speciesMolality(Mg++ units=mmolal)", "Mg++ [mmol]");
        output.add("speciesMolality(Ca++ units=mmolal)", "Ca++ [mmol]");
        output.add("speciesMolality(CO3-- units=mmolal)", "CO3-- [mmol]");
        output.add("speciesMolality(K+ units=mmolal)", "K+ [mmol]");
        output.add("speciesMolality(Al+++ units=mmolal)", "Al+++ [mmol]");
        output.add("speciesMolality(Calcite units=molal)", "Calcite [mol]");
        output.add("speciesMolality(Dolomite units=molal)", "Dolomite [mol]");
        output.add("speciesMolality(Halite units=molal)", "Halite [mol]");
        output.add("speciesMolality(K-Feldspar units=molal)", "K-Feldspar [mol]");
        output.add("speciesMolality(Quartz units=molal)", "Quartz [mol]");
        output.add("speciesMolality(Kaolinite units=molal)", "Kaolinite [mol]");
        output.add("pH");

        tic(KINETIC_PATH);

        //path.solve(state0, t0, t1, "second");
        path.solve(state0, t0, dt, n, "second");

        auto kinetic_path_time = toc(KINETIC_PATH);
        std::cout << "----------------------------------------" << std::endl;
        std::cout << "# kin : " << num_kinetic_minerals
                  << ", # eq  : " << num_minerals - num_kinetic_minerals
                  << ", time  : " << kinetic_path_time << " s" << std::endl;
        std::cout << "*******************************************************************************" << std::endl;
    }

}
