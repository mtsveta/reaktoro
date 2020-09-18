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
    double t0 = 0;
    double dt = 0.1; // in seconds
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

    double initial_mineral_amount_mol = 10.0;

    // Step **: Define chemical equilibrium solver options
    EquilibriumOptions equilibrium_options;
    equilibrium_options.hessian = GibbsHessian::Exact;

    // Step **: Define smart chemical equilibrium solver options
    SmartEquilibriumOptions smart_equilibrium_options;
    smart_equilibrium_options.reltol = 1e-3;
    smart_equilibrium_options.learning.hessian = GibbsHessian::Exact;
    smart_equilibrium_options.smart_method = "kin-clustering-eq-clustering";

    // Step **: Define chemical kinetic solver options
    KineticOptions kinetic_options;
    kinetic_options.equilibrium = equilibrium_options;
    kinetic_options.smart_equilibrium = smart_equilibrium_options;
    kinetic_options.use_smart_equilibrium_solver = false;

    // Step **: Define smart chemical kinetic solver options
    SmartKineticOptions smart_kinetic_options;
    smart_kinetic_options.reltol = 1e-1;
    smart_kinetic_options.abstol = 1e-4;
    smart_kinetic_options.tol = 1e-3;
    smart_kinetic_options.learning = kinetic_options;
    smart_kinetic_options.learning.equilibrium = equilibrium_options;
    smart_kinetic_options.smart_equilibrium = smart_equilibrium_options;
    smart_kinetic_options.use_smart_equilibrium_solver = false;
    smart_kinetic_options.smart_method = "kin-clustering-eq-clustering";


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

    Partition partition(system);
    partition.setKineticPhases({"Halite", "Calcite", "Dolomite", "K-Feldspar", "Quartz", "Kaolinite"});

    double T = 61.0 + 273.15;       // temperature (in units of celsius)
    double P = 1 * 1.01325 * 1e5;   // pressure (in units of atm)
    EquilibriumInverseProblem problem(partition);
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

    //problem.pE(4.0); // the result of PHREEQC has pE = 4.0

    EquilibriumInverseSolver solver(partition);
    ChemicalState state0(system);
    EquilibriumResult result = solver.solve(state0, problem);

    //std::cout << "state0 = \n" << state0 << std::endl;

    ChemicalProperties props = state0.properties();
    auto evaluate_pH = ChemicalProperty::pH(system);
    auto evaluate_alkalinity = ChemicalProperty::alkalinity(system);
    auto pH = evaluate_pH(props);
    auto alkalinity = evaluate_alkalinity(props);

    std::cout << "# iter = " << result.optimum.iterations << std::endl;
    std::cout << "pH                 = " << pH.val << std::endl;
    std::cout << "Alkalinity [eq/L]  = " << alkalinity << std::endl;

    /*
    state0 =
    =================================================================================================================================================================================================================================
    Temperature [K]          Temperature [C]          Pressure [Pa]            Pressure [bar]
    ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    334.15                   61                       101325                   1.01325
    =================================================================================================================================================================================================================================
    Element                  Amount [mol]             Aqueous [mol]            Halite [mol]             Calcite [mol]            Dolomite [mol]           K-Feldspar [mol]         Quartz [mol]             Kaolinite [mol]          Dual Potential [kJ/mol]
    ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Al                       0.0001                   0.0001                   0                        0                        0                        1e-10                    0                        2e-10                    -430.233
    C                        0.0001                   0.0001                   0                        1e-10                    2e-10                    0                        0                        0                        49.173
    Ca                       0.0001                   0.0001                   0                        1e-10                    1e-10                    0                        0                        0                        -498.783
    Cl                       6.27                     6.27                     1e-10                    0                        0                        0                        0                        0                        -172.304
    H                        111.017                  111.017                  0                        0                        0                        0                        0                        4e-10                    -2.35331
    K                        1.1e-09                  1e-09                    0                        0                        0                        1e-10                    0                        0                        -303.023
    Mg                       0.0001                   0.0001                   0                        0                        1e-10                    0                        0                        0                        -396.675
    Na                       6.27                     6.27                     1e-10                    0                        0                        0                        0                        0                        -218.465
    O                        55.5084                  55.5084                  0                        3e-10                    6e-10                    8e-10                    2e-10                    9e-10                    -235.731
    S                        5.7e-19                  5.7e-19                  0                        0                        0                        0                        0                        0                        0
    Si                       0.000100001              0.0001                   0                        0                        0                        3e-10                    1e-10                    2e-10                    -390.838
    Z                        0.00017249               0.00017249               0                        0                        0                        0                        0                        0                        -42.4272
    =================================================================================================================================================================================================================================
    Species                  Amount [mol]             Mole Fraction [mol/mol]  Activity Coefficient [-] Activity [-]             Potential [kJ/mol]       Dual Potential [kJ/mol]
    ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Al+++                    9.89996e-08              1.50498e-09              3.93858e-07              3.89923e-14              -557.515                 2.80635e-13
    AlO+                     3.85672e-08              5.86294e-10              0.292194                 1.12693e-08              -708.392                 7.20374e-13
    AlO2-                    9.60796e-05              1.46059e-06              0.292196                 2.80744e-05              -859.268                 2.89164e-16
    AlOH++                   7.36002e-09              1.11886e-10              0.00401324               2.95379e-11              -753.172                 3.77483e-12
    CO(aq)                   1.83254e-10              2.78581e-12              1                        1.83257e-10              -186.558                 1.51608e-10
    CO2(aq)                  3.76868e-06              5.72911e-08              3                        1.13062e-05              -422.289                 7.37202e-15
    CO3--                    1.01827e-06              1.54796e-08              0.0358628                3.65184e-08              -573.166                 2.72844e-14
    Ca(HCO3)+                1.79306e-08              2.72579e-10              0.292194                 5.23929e-09              -1201.58                 1.54946e-12
    Ca(HSiO3)+               1.08558e-10              1.65029e-12              0.292194                 3.17205e-11              -1641.59                 2.55925e-10
    Ca++                     1.9264e-05               2.92849e-07              0.36756                  7.08075e-06              -583.637                 1.44221e-15
    CaCO3(aq)                1.27512e-09              1.93842e-11              1                        1.27514e-09              -1156.8                  2.17884e-11
    CaCl+                    6.36868e-05              9.6816e-07               0.292194                 1.86092e-05              -713.513                 4.36241e-16
    CaCl2(aq)                1.70291e-05              2.58875e-07              1                        1.70293e-05              -843.39                  1.63149e-15
    CaOH+                    8.21547e-10              1.24891e-11              0.292194                 2.40055e-10              -779.294                 3.38176e-11
    CaSO4(aq)                1e-20                    1.52019e-22              1                        1.00001e-20              -1437.8                  0
    Cl-                      4.00233                  0.0608431                0.788966                 3.15775                  -129.876                 6.94164e-21
    ClO-                     1.50591e-22              2.28927e-24              0.292196                 4.40026e-23              -181.116                 184.492
    ClO2-                    5.91892e-23              8.99788e-25              0.292196                 1.72951e-23              -131.95                  469.389
    ClO3-                    4.10539e-23              6.24097e-25              0.292196                 1.19959e-23              -160.33                  676.74
    ClO4-                    3.0521e-23               4.63977e-25              0.292196                 8.91822e-24              -162.516                 910.285
    H+                       5.23784e-08              7.96251e-10              1.90916                  1e-07                    -44.7805                 5.30424e-13
    H2(aq)                   0.00075                  1.14014e-05              1                        0.000750011              -4.70661                 3.70437e-17
    H2O(l)                   55.5077                  0.843822                 0.967394                 0.816309                 -240.438                 5.00521e-22
    H2O2(aq)                 1.40968e-22              2.14298e-24              1                        1.4097e-22               -279.083                 197.086
    H2S(aq)                  1e-20                    1.52019e-22              1                        1.00001e-20              -160.773                 0
    H2S2O3(aq)               1e-20                    1.52019e-22              1                        1.00001e-20              -670.529                 0
    H2S2O4(aq)               1e-20                    1.52019e-22              1                        1.00001e-20              -752.688                 0
    HAlO2(aq)                3.77545e-06              5.7394e-08               1                        3.7755e-06               -904.049                 7.3588e-15
    HCO3-                    9.51723e-05              1.4468e-06               0.523504                 4.98238e-05              -617.946                 2.91921e-16
    HCl(aq)                  5.06501e-08              7.69978e-10              1                        5.06508e-08              -174.657                 5.48524e-13
    HClO(aq)                 1.49649e-22              2.27495e-24              1                        1.49651e-22              -224.735                 185.653
    HClO2(aq)                5.52827e-23              8.40401e-25              1                        5.52834e-23              -143.561                 502.559
    HO2-                     1.26278e-22              1.91967e-24              0.292196                 3.68985e-23              -211.376                 220.013
    HS-                      1e-20                    1.52019e-22              0.292196                 2.922e-21                -121.66                  0
    HS2O3-                   1e-20                    1.52019e-22              0.292196                 2.922e-21                -668.203                 0
    HS2O4-                   1e-20                    1.52019e-22              0.292196                 2.922e-21                -751.617                 0
    HSO3-                    1e-20                    1.52019e-22              0.292196                 2.922e-21                -664.121                 0
    HSO4-                    1e-20                    1.52019e-22              0.522344                 5.22351e-21              -890.08                  0
    HSO5-                    1e-20                    1.52019e-22              0.292196                 2.922e-21                -776.845                 0
    HSiO3-                   1.11635e-06              1.69706e-08              0.292196                 3.26197e-07              -1057.96                 2.48872e-14
    K+                       9.87448e-10              1.50111e-11              0.539955                 5.33185e-10              -345.45                  2.8136e-11
    KCl(aq)                  1.25523e-11              1.90818e-13              1                        1.25524e-11              -475.326                 2.21337e-09
    KHSO4(aq)                1e-20                    1.52019e-22              1                        1.00001e-20              -1155.13                 0
    KOH(aq)                  2.222e-16                3.37787e-18              1                        2.22204e-16              -541.107                 0.000125035
    KSO4-                    1e-20                    1.52019e-22              0.292196                 2.922e-21                -1168.49                 0
    Mg(HCO3)+                2.08489e-08              3.16943e-10              0.292194                 6.09202e-09              -1099.48                 1.33258e-12
    Mg(HSiO3)+               1.85483e-10              2.81969e-12              0.292194                 5.41978e-11              -1539.49                 1.49786e-10
    Mg++                     1.72811e-05              2.62705e-07              0.473801                 8.1879e-06               -481.529                 1.6077e-15
    MgCO3(aq)                4.89943e-10              7.44806e-12              1                        4.89949e-10              -1054.7                  5.67062e-11
    MgCl+                    8.26936e-05              1.2571e-06               0.292194                 2.41629e-05              -611.406                 3.35973e-16
    MgOH+                    3.82816e-09              5.81953e-11              0.591789                 2.26549e-09              -677.186                 7.25748e-12
    MgSO4(aq)                1e-20                    1.52019e-22              1                        1.00001e-20              -1336.91                 0
    Na+                      4.00248                  0.0608453                0.788914                 3.15766                  -260.892                 6.94139e-21
    NaCl(aq)                 2.26749                  0.0344701                1                        2.26752                  -390.768                 1.22527e-20
    NaHSiO3(aq)              3.03392e-05              4.61213e-07              1                        3.03396e-05              -1318.85                 9.15739e-16
    NaOH(aq)                 1.66518e-06              2.53139e-08              1                        1.6652e-06               -456.549                 1.66846e-14
    NaSO4-                   1e-20                    1.52019e-22              0.292196                 2.922e-21                -1145.5                  0
    O2(aq)                   8.11732e-23              1.23399e-24              1                        8.11743e-23              -129.197                 342.265
    OH-                      1.5076e-06               2.29184e-08              0.528428                 7.96672e-07              -195.657                 1.84284e-14
    S2--                     1e-20                    1.52019e-22              0.0040133                4.01335e-23              -64.2461                 0
    S2O3--                   1e-20                    1.52019e-22              0.0040133                4.01335e-23              -667.777                 0
    S2O4--                   1e-20                    1.52019e-22              0.0040133                4.01335e-23              -746.566                 0
    S2O5--                   1e-20                    1.52019e-22              0.0040133                4.01335e-23              -937.405                 0
    S2O6--                   1e-20                    1.52019e-22              0.0040133                4.01335e-23              -1113.91                 0
    S2O8--                   1e-20                    1.52019e-22              0.0040133                4.01335e-23              -1266.9                  0
    S3--                     1e-20                    1.52019e-22              0.0040133                4.01335e-23              -71.5219                 0
    S3O6--                   1e-20                    1.52019e-22              0.0040133                4.01335e-23              -1106.01                 0
    S4--                     1e-20                    1.52019e-22              0.0040133                4.01335e-23              -77.5273                 0
    S4O6--                   1e-20                    1.52019e-22              0.0040133                4.01335e-23              -1192.88                 0
    S5--                     1e-20                    1.52019e-22              0.0040133                4.01335e-23              -82.2767                 0
    S5O6--                   1e-20                    1.52019e-22              0.0040133                4.01335e-23              -1107.1                  0
    SO2(aq)                  1e-20                    1.52019e-22              1                        1.00001e-20              -435.356                 0
    SO3--                    1e-20                    1.52019e-22              0.0040133                4.01335e-23              -628.241                 0
    SO4--                    1e-20                    1.52019e-22              0.0319362                3.19367e-22              -882.113                 0
    SiO2(aq)                 6.85442e-05              1.042e-06                1                        6.85451e-05              -862.3                   4.05327e-16
    Halite                   1e-10                    1                        1                        1                        -386.836                 3.93251
    Calcite                  1e-10                    1                        1                        1                        -1132.71                 24.0922
    Dolomite                 1e-10                    1                        1                        1                        -2172.27                 39.2262
    K-Feldspar               1e-10                    1                        1                        1                        -3754.4                  37.2138
    Quartz                   1e-10                    1                        1                        1                        -857.834                 4.46604
    Kaolinite                1e-10                    1                        1                        1                        -3796.98                 -23.8434
    =================================================================================================================================================================================================================================
    Phase                    Amount [mol]             Stability                Stability Index [-]      Mass [kg]                Volume [m3]              Density [kg/m3]          Molar Volume [m3/mol]    Volume Fraction [m3/m3]
    ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Aqueous                  65.7813                  stable                   -9.64327e-17             1.36645                  0.00117003               1167.87                  1.77867e-05              1
    Halite                   1e-10                    under stable             -0.614721                5.84428e-12              2.7015e-15               2163.35                  2.7015e-05               2.30891e-12
    Calcite                  1e-10                    under stable             -3.76604                 1.00087e-11              3.6934e-15               2709.89                  3.6934e-05               3.15666e-12
    Dolomite                 1e-10                    under stable             -6.13176                 1.84401e-11              6.4365e-15               2864.92                  6.4365e-05               5.50112e-12
    K-Feldspar               1e-10                    under stable             -5.81718                 2.78332e-11              1.0887e-14               2556.55                  0.00010887               9.30486e-12
    Quartz                   1e-10                    under stable             -0.698121                6.00843e-12              2.2688e-15               2648.29                  2.2688e-05               1.93909e-12
    Kaolinite                1e-10                    super stable              3.72715                 2.5816e-11               9.952e-15                2594.06                  9.952e-05                8.50574e-12
    =================================================================================================================================================================================================================================
    Ionic Strength [molal]   pH                       pE                       Reduction Potential [V]  Alkalinity [eq/L]
    ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    4.00271                  7                        -0.284906                -0.01889                 0.000159835
    =================================================================================================================================================================================================================================
    */

    // Differences         of    PHREEQC   vs Reaktoro
    // Ionic Strength [molal]    6.27         4.00271
    // Total alkalinity (eq/kg)  4.949e-04    0.000159835
    // # iterations              11           74
    // pe                        4.0          -0.284906
    // Because of the added water in Reaktoro:
    // H(0)                      2.711e-25    111.017
    // O(0)                      5.309e-34    55.5084
    // even though
    // H2O                       5.553e+01    55.5077
    //
    // Which species can be crucial for the comparison?

    state0.setSpeciesAmount("Calcite", initial_mineral_amount_mol, "mol");
    state0.setSpeciesAmount("Quartz", initial_mineral_amount_mol, "mol");
    state0.setSpeciesAmount("Halite", initial_mineral_amount_mol, "mol");
    state0.setSpeciesAmount("Dolomite", initial_mineral_amount_mol, "mol");
    state0.setSpeciesAmount("Kaolinite", initial_mineral_amount_mol, "mol");
    state0.setSpeciesAmount("K-Feldspar", initial_mineral_amount_mol, "mol");

    KineticPathOptions kinetic_path_options;
    kinetic_path_options.use_smart_equilibrium_solver = kinetic_options.use_smart_equilibrium_solver;
    kinetic_path_options.use_smart_kinetic_solver = false;
    kinetic_path_options.equilibrium = equilibrium_options;
    kinetic_path_options.smart_equilibrium = smart_equilibrium_options;
    kinetic_path_options.kinetics = kinetic_options; // TODO: think about better structure of the kinetic and equilibrium options
    kinetic_path_options.smart_kinetics = smart_kinetic_options;

    KineticPath path(reactions, partition);
    path.setOptions(kinetic_path_options);

    std::string result_filename = "kinetics-simple-benchmark-6-kin-0-eq-dt-" + std::to_string(dt) +
                                  (kinetic_path_options.use_smart_kinetic_solver ? "-smart-kin" : "-conv-kin") +
                                  (kinetic_path_options.use_smart_equilibrium_solver ? "-smart-eq"  : "-conv-eq") +
                                  "-n-" +  std::to_string(n) +
                                  "-ncells-" +  std::to_string(ncells) +
                                  "-activity-model-" + activity_model +
                                  "-uniform-sa-1.0-m2mol.txt";

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

    ChemicalState state0_tmp = state0;

    path.solve(state0, t0, dt, n, "second");

    auto kinetic_path_time = toc(KINETIC_PATH);

    double delta_halite = state0.speciesAmount("Halite") - initial_mineral_amount_mol;
    double delta_calcite = state0.speciesAmount("Calcite") - initial_mineral_amount_mol;
    double delta_dolomite = state0.speciesAmount("Dolomite") - initial_mineral_amount_mol;
    double delta_kfeldspar = state0.speciesAmount("K-Feldspar") - initial_mineral_amount_mol;
    double delta_quartz = state0.speciesAmount("Quartz") - initial_mineral_amount_mol;
    double delta_kaolinite = state0.speciesAmount("Kaolinite") - initial_mineral_amount_mol;

    std::cout << "delta_halite = " << delta_halite << std::endl;
    std::cout << "delta_calcite = " << delta_calcite << std::endl;
    std::cout << "delta_dolomite = " << delta_dolomite << std::endl;
    std::cout << "delta_kfeldspar = " << delta_kfeldspar << std::endl;
    std::cout << "delta_quartz = " << delta_quartz << std::endl;
    std::cout << "delta_kaolinite = " << delta_kaolinite << std::endl << std::endl;

    std::cout << "total kinetic path time = " << kinetic_path_time << " s" << std::endl;

    //    Differences       PHREEQC   vs Reaktoro
    // --------------------------------------------------
    //    Halite           -2.442e-01    -0.606381
    //    Calcite          -4.359e-06   -4.3867e-06
    //    Dolomite-dis     -1.068e-06   -1.22795e-06
    //    K-Feldspar       -2.139e-12   -2.28084e-12
    //    Quartz           -1.776e-12   -3.73035e-13
    //    Kaolinite         1.361e-09    7.28797e-10
}
