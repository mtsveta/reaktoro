#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

#include <ThermoFun/ThermoFun.h>

int main()
{
    ThermoFun::Database database("databases/thermofun/cemdata18-thermofun.json");

    double T = 293.15;
    double P = 100000.0;
    double ncells = 100;
    double xl = 0.0;
    double xr = 1.0;

    ChemicalEditor editor(database);
    editor.setTemperatures({T}, "kelvin");
    editor.setPressures({P}, "pascal");
    editor.addAqueousPhase({"Al(SO4)+", "Al(SO4)2-", "Al+3", "AlO+", "AlO2-", "AlO2H@", "AlOH+2", "AlSiO5-3", "Ca(CO3)@", "Ca(HCO3)+", "Ca(SO4)@", "Ca+2", "CaOH+", "Ca(HSiO3)+", "CaSiO3@", "Fe(CO3)@", "Fe(HCO3)+", "Fe(HSO4)+", "Fe(SO4)@", "Fe+2", "FeCl+", "FeOH+", "Fe(HSO4)+2", "Fe(SO4)+", "Fe(SO4)2-", "Fe+3", "FeCl+2", "FeCl2+", "FeCl3@", "FeO+", "FeO2-", "FeO2H@", "FeOH+2", "K(SO4)-", "K+", "KOH@", "Mg(CO3)@", "Mg(HCO3)+", "Mg+2", "MgOH+", "MgSO4@", "Mg(HSiO3)+", "Na(CO3)-", "Na(HCO3)@", "Na(SO4)-", "Na+", "NaOH@", "HSiO3-", "Si4O10-4", "SiO2@", "SiO3-2", "CO2@", "CO3-2", "HCO3-", "CH4@", "ClO4-", "Cl-", "H2@", "N2@", "O2@", "S2O3-2", "HSO3-", "SO3-2", "HSO4-", "SO4-2", "H2S@", "HS-", "OH-", "H+", "H2O@"});
    editor.addGaseousPhase({"CO2", "CH4", "H2", "N2", "O2", "H2S", "H2O"});

    // "FeCO3(pr)",
    vector<std::string> minerals = {"monocarbonate", "C2AClH5", "C3AFS0.84H4.32", "C3FS0.84H4.32", "CSHQ-JenD", "CSHQ-JenH", "CSHQ-TobD", "CSHQ-TobH", "KSiOH", "NaSiOH", "straetlingite", "straetlingite7", "C4AH13", "monosulphate12", "tricarboalu03", "ettringite03_ss", "M075SH", "M15SH", "AlOHmic", "Kln", "C12A7", "C3A", "C3S", "C4AF", "CA", "CA2", "C2AH7.5", "C3AH6", "C4AH11", "C4AH13", "C4AH19", "CAH10", "monosulphate10.5", "monosulphate12", "monosulphate14", "monosulphate16", "monosulphate9", "chabazite", "zeoliteP_Ca", "straetlingite5.5", "monocarbonate9", "hemicarbonat10.5", "hemicarbonate", "hemicarbonate9", "monocarbonate", "C4AsClH12", "ettringite13", "ettringite9", "Arg", "Cal", "C3FH6", "C4FH13", "Fe-hemicarbonate", "Femonocarbonate", "Lim", "Portlandite", "Anh", "Gp", "hemihydrate", "Sd", "Mag", "FeOOHmic", "Py", "Tro", "Melanterite", "K2SO4", "syngenite", "hydrotalcite", "Brc", "Na2SO4", "natrolite", "zeoliteX", "zeoliteY", "Sulfur", "Qtz", "Amor-Sl"};
    for (auto minerial : minerals) editor.addMineralPhase(minerial);

    // Step **: Create the ChemicalSystem object using the configured editor
    ChemicalSystem system(editor);
    std::cout << "system = \n" << system << std:: endl;

    /*
    // Step **: Define the initial condition (IC) of the reactive transport modeling problem
    EquilibriumProblem problem_ic(system);
    problem_ic.setTemperature(T, "celsius");
    problem_ic.setPressure(P, "bar");
    problem_ic.add("H2O",   1.0, "kg");
    problem_ic.add("NaCl",  0.7, "mol");
    problem_ic.add("CaCO3", 10,  "mol");
    problem_ic.add("SiO2",  10,  "mol");
    */
    /*
    EquilibriumProblem problem(system);
    problem.add("CO2",             	0.001,  "g");
    problem.add("CaCO3",           	1,	    "g");
    problem.add("MgSiO3",          	1,	    "g");
    problem.add("NaCl",            	5,      "g");
    problem.add("NaAlSi3O8",        37,	    "g");
    problem.add("KAl3Si3O10(OH)2",  13,	    "g");
    problem.add("SiO2",          	30,	    "g");
    problem.add("KAlSi3O8",        	20,	    "g");
    */

    /*
    // Step **: Define the boundary condition (BC) of the reactive transport modeling problem with seawater
    EquilibriumProblem problem_bc(system);
    problem_bc.setTemperature(T, "celsius");
    problem_bc.setPressure(P, "bar");
    problem_bc.add("H2O",   1.00, "kg");
    problem_bc.add("NaCl",  0.90, "mol");
    problem_bc.add("MgCl2", 0.05, "mol");
    problem_bc.add("CaCl2", 0.01, "mol");
    problem_bc.add("CO2",   0.75, "mol");

    // Step **: Calculate the equilibrium states for the IC and BC
    ChemicalState state_ic = equilibrate(problem_ic);
    ChemicalState state_bc = equilibrate(problem_bc);

    state_ic.output("state_ic.txt");
    state_bc.output("state_bc.txt");

    // Step **: Scale the boundary condition state
    state_bc.scaleVolume(1.0, "m3");

    // Step **: Scale the volumes of the phases in the initial condition
    state_ic.scalePhaseVolume("Aqueous", 0.1, "m3");    // porosity 10% if the 1.0m3
    state_ic.scalePhaseVolume("Quartz", 0.882, "m3");   // 0.882 = 0.98 * 0.9 (0.9 is 90% of 1.0m3, 0.98 is 98% quartz of the rock)
    state_ic.scalePhaseVolume("Calcite", 0.018, "m3");  // 0.018 = 0.02 * 0.9 (0.9 is 90% of 1.0m3, 0.02 is 2% calcite of the rock)

    // Step **: Create the mesh for the column
    Mesh mesh(ncells, xl, xr);

    // Step **: Create a chemical field object with every cell having state given by state_ic
    ChemicalField field(mesh.numCells(), state_ic);
     */
}
