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
    int n = 1;
    Index ncells = 1;

    double initial_sa_m2mol = 1.0; // defined in PHREEQC, m2/mol
    double initial_sa_m2kg;

    double molar_mass_halite = 58.44; // g / mol
    double molar_mass_calcite = 100.0869; // g / mol
    double molar_mass_dolomite = 184.40; // g / mol
    double molar_mass_quartz = 60.083; // g / mol
    double molar_mass_kfeldspar = 278.3315; // g / mol
    double molar_mass_kaolinite = 258.071; // g / mol // or 258.1604

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

    initial_sa_m2kg = initial_sa_m2mol / molar_mass_halite / 1e-3;
    editor.addMineralReaction("Halite")
            .setEquation("Halite = Na+ + Cl-")
            .addMechanism("logk = -0.21 mol/(m2*s); Ea = 7.4 kJ/mol")
            .setSpecificSurfaceArea(initial_sa_m2kg, "m2/kg"); // conversion from

    initial_sa_m2kg = initial_sa_m2mol / molar_mass_calcite / 1e-3;
    editor.addMineralReaction("Calcite")
            .setEquation("Calcite = Ca++ + CO3--")
            .addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol")
            .addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0")
            .setSpecificSurfaceArea(initial_sa_m2kg, "m2/kg");

    initial_sa_m2kg = initial_sa_m2mol / molar_mass_dolomite / 1e-3;
    editor.addMineralReaction("Dolomite")
            .setEquation("Dolomite = Ca++ + Mg++ + 2*CO3--")
            .addMechanism("logk = -7.53 mol/(m2*s); Ea = 52.2 kJ/mol")
            .addMechanism("logk = -3.19 mol/(m2*s); Ea = 36.1 kJ/mol; a[H+] = 0.5")
            .setSpecificSurfaceArea(initial_sa_m2kg, "m2/kg");

    // Quartz + 2*H2O = H4SiO4(aq)
    initial_sa_m2kg = initial_sa_m2mol / molar_mass_quartz / 1e-3;
    editor.addMineralReaction("Quartz")
            .setEquation("Quartz + H2O(l) = H+ + HSiO3-")
            .addMechanism("logk = -13.99 mol/(m2*s); Ea = 87.7 kJ/mol")
            .setSpecificSurfaceArea(initial_sa_m2kg, "m2/kg");

    // K-Feldspar: K(AlSi3)O8 # potassium feldspar, orthoclase
    initial_sa_m2kg = initial_sa_m2mol / molar_mass_kfeldspar / 1e-3;
    editor.addMineralReaction("K-Feldspar")
            .setEquation("K-Feldspar + H2O(l) + H+ = Al+++ + 3*HSiO3- + K+")
            .addMechanism("logk = -12.41 mol/(m2*s); Ea = 38.0 kJ/mol")
            .addMechanism("logk = -10.06 mol/(m2*s); Ea = 51.7 kJ/mol; a[H+] = 0.5")
            .setSpecificSurfaceArea(initial_sa_m2kg, "m2/kg");

    // Kaolinite: Al2Si2O5(OH)4
    initial_sa_m2kg = initial_sa_m2mol / molar_mass_kaolinite / 1e-3;
    editor.addMineralReaction("Kaolinite")
            .setEquation("Kaolinite + 4*H+ = 3*H2O(l) + 2*HSiO3- + 2*Al+++")
            .addMechanism("logk = -13.18 mol/(m2*s); Ea = 22.2 kJ/mol")
            .addMechanism("logk = -11.31 mol/(m2*s); Ea = 65.9 kJ/mol; a[H+] = 0.777")
            .setSpecificSurfaceArea(initial_sa_m2kg, "m2/kg");

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


        EquilibriumOptions options;
        options.hessian = GibbsHessian::ApproximationDiagonal;

        EquilibriumInverseSolver solver(partitions[i]);
        solver.setOptions(options);
        ChemicalState state0(system);
        EquilibriumResult result = solver.solve(state0, problem);

        std::cout << "# iter (in initial state) = " << result.optimum.iterations << std::endl;

        double initial_mineral_amount_mol = 10.0;
        state0.setSpeciesAmount("Calcite", initial_mineral_amount_mol, "mol");
        state0.setSpeciesAmount("Quartz", initial_mineral_amount_mol, "mol");
        state0.setSpeciesAmount("Halite", initial_mineral_amount_mol, "mol");
        state0.setSpeciesAmount("Dolomite", initial_mineral_amount_mol, "mol");
        state0.setSpeciesAmount("Kaolinite", initial_mineral_amount_mol, "mol");
        state0.setSpeciesAmount("K-Feldspar", initial_mineral_amount_mol, "mol");

        KineticOptions kinetic_options;
        kinetic_options.equilibrium = options;

        KineticPath path(reactions, partitions[i]);
        path.setOptions(kinetic_options);

        std::string result_filename = "kinetics-benchmark-kin-" + std::to_string(num_kinetic_minerals)
                                                           + "-eq-" + std::to_string(num_minerals - num_kinetic_minerals)
                                                           + "-dt-" + std::to_string(dt)
                                                           + "-n-" +  std::to_string(n) +
                                                           + "-ncells-" +  std::to_string(ncells) +
                                                           + "-activity-model-" + activity_model +
                                                           + "-uniform-sa-1.0-m2mol.txt";
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

        path.solve(state0, t0, dt, n, "second");

        auto kinetic_path_time = toc(KINETIC_PATH);
        std::cout << "----------------------------------------" << std::endl;
        std::cout << "# kin : " << num_kinetic_minerals
                  << ", # eq  : " << num_minerals - num_kinetic_minerals
                  << ", time  : " << kinetic_path_time << " s" << std::endl;
        std::cout << "*******************************************************************************" << std::endl;
    }

}
