#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

#include <ThermoFun/ThermoFun.h>

int main()
{
    //ThermoFun::Database database("~/polybox/reaktoro/databases/thermofun/cemdata18-thermofun.json");
    ThermoFun::Database database("../databases/thermofun/cemdata18-thermofun.json");

    double T = 293.15;
    double P = 100000.0;
    double ncells = 100;
    double xl = 0.0;
    double xr = 1.0;

    ChemicalEditor editor(database);
    editor.setTemperatures({T}, "kelvin");
    editor.setPressures({P}, "pascal");
    editor.addAqueousPhase({"Al(SO4)+", "Al(SO4)2-", "Al+3", "AlO+", "AlO2-", "AlO2H@", "AlOH+2", "AlSiO5-3", "Ca(CO3)@", "Ca(HCO3)+", "Ca(SO4)@", "Ca+2", "CaOH+", "Ca(HSiO3)+", "CaSiO3@", "Fe(CO3)@", "Fe(HCO3)+", "Fe(HSO4)+", "Fe(SO4)@", "Fe+2", "FeCl+", "FeOH+", "Fe(HSO4)+2", "Fe(SO4)+", "Fe(SO4)2-", "Fe+3", "FeCl+2", "FeCl2+", "FeCl3@", "FeO+", "FeO2-", "FeO2H@", "FeOH+2", "K(SO4)-", "K+", "KOH@", "Mg(CO3)@", "Mg(HCO3)+", "Mg+2", "MgOH+", "MgSO4@", "Mg(HSiO3)+", "Na(CO3)-", "Na(HCO3)@", "Na(SO4)-", "Na+", "NaOH@", "HSiO3-", "Si4O10-4", "SiO2@", "SiO3-2", "CO2@", "CO3-2", "HCO3-", "CH4@", "ClO4-", "Cl-", "H2@", "N2@", "O2@", "S2O3-2", "HSO3-", "SO3-2", "HSO4-", "SO4-2", "H2S@", "HS-", "OH-", "H+", "H2O@"});
    editor.addGaseousPhase({"CO2(g)", "CH4(g)", "H2(g)", "N2(g)", "O2(g)", "H2S(g)", "H2O(g)"});
    std::vector<std::string> minerals = {"monocarbonate", "C2AClH5", "C3AFS0.84H4.32", "C3FS0.84H4.32", "CSHQ-JenD", "CSHQ-JenH", "CSHQ-TobD", "CSHQ-TobH", "KSiOH", "NaSiOH", "straetlingite", "straetlingite7", "C4AH13", "monosulphate12", "tricarboalu03", "ettringite03_ss", "M075SH", "M15SH", "AlOHmic", "Kln", "C12A7", "C3A", "C3S", "C4AF", "CA", "CA2", "C2AH7.5", "C3AH6", "C4AH11", "C4AH13", "C4AH19", "CAH10", "monosulphate10.5", "monosulphate12", "monosulphate14", "monosulphate16", "monosulphate9", "chabazite", "zeoliteP_Ca", "straetlingite5.5", "monocarbonate9", "hemicarbonat10.5", "hemicarbonate", "hemicarbonate9", "monocarbonate", "C4AsClH12", "ettringite13", "ettringite9", "Arg", "Cal", "C3FH6", "C4FH13", "Fe-hemicarbonate", "Femonocarbonate", "Lim", "Portlandite", "Anh", "Gp", "hemihydrate", "FeCO3(pr)", "Sd", "Mag", "FeOOHmic", "Py", "Tro", "Melanterite", "K2SO4", "syngenite", "hydrotalcite", "Brc", "Na2SO4", "natrolite", "zeoliteX", "zeoliteY", "Sulfur", "Qtz", "Amor-Sl"};
    for (auto minerial : minerals) editor.addMineralPhase(minerial);

    // Step **: Create the ChemicalSystem object using the configured editor
    ChemicalSystem system(editor);
    Index num_elements = system.numElements();
    Index num_species = system.numSpecies();
    std::cout << "system = \n" << system << std:: endl;
    std::cout << "Chemical system contains: \n" << " - " << num_elements << " elements, \n"
                                                << " - " << num_species << " species, \n"
                                                << " - " << system.numPhases() << " phases." << std::endl;
    EquilibriumSolver solver(system);
    // From nSWp-dbr-0-0000.dat
    // bIC: Bulk composition of reactive subsystem (main GEM input), moles of ICs [nICb]
    Vector b_sw(num_elements);
    b_sw << 1e-12, 8.39047799118904e-07, 4.60274559011631e-06, 0.00024470481263518, 1e-12, 0.0480862507085869, 4.5682055764741e-06, 2.37779781110928e-05, 0.000209584727160819, 5.4730517594269e-08, 0.239750374257753, 1.26338032622899e-05, 0.107827168643251, 0;
    /*
    Vector n_sw(num_species);
    n_sw << 2.46100221484366e-21, 7.49133660286559e-22, 5.75242310773387e-21, 6.90681743697294e-17, 9.90800739368853e-13, 9.12927398011268e-15, 8.56314982070215e-19, 4.72498975561693e-20, 5.94922644003287e-21, 8.54248693442392e-09, 1.28390023290633e-08, 3.45691884021378e-07, 4.23557433459937e-06, 2.05232883296596e-11, 7.63063575495622e-11, 1.0525861918408e-12, 4.39471302392405e-25, 3.63567913737693e-25, 4.80923510908748e-32, 9.91586748163519e-25, 1.39895482718739e-23, 1.90008627164411e-24, 1.53430533643884e-25, 2.37937181783215e-30, 3.39737043871988e-23, 2.2047927969019e-24, 5.7209131449752e-23, 1.19305043639265e-28, 3.75830040054004e-34, 6.33502555926846e-23, 2.9355392326128e-23, 8.30349031726048e-25, 6.72188840987796e-14, 1.85696519699273e-13, 7.47079270448143e-13, 5.20123301714824e-18, 1.2431602548248e-19, 9.06591037658476e-08, 4.47754516670466e-06, 1.3060035882638e-12,
            2.61236812901651e-08, 5.98734301787103e-08, 2.16483206816277e-05, 2.55807225290953e-09, 2.04024260738186e-06, 8.03111666803934e-10, 5.65266946690013e-11, 4.20224888575075e-08, 7.27152956534961e-08, 2.95494185377378e-06, 0.000206514926295945, 1.21226589402462e-10, 2.82491060715618e-09, 2.28299242396247e-19, 6.4881254166339e-08, 8.93964212306933e-14, 4.14429201643602e-09, 2.30155215456548e-08, 5.89771600313447e-07, 0, 2.89784766978446e-27, 0.00024470481263518, 0, 2.73652587971345e-08, 1.13662409601231e-07, 0, 0, 0, 6.60908355733222e-13, 7.20226715243868e-06, 0, 0, 0, 8.34730158698885e-10, 2.67834682848817e-12, 0.0240427541324881, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1078271, 0;
    */
    ChemicalState state_sw(system);
    solver.solve(state_sw, T, P, b_sw);
    std::cout << "state_sw = \n" << state_sw << std:: endl;

    // bSP: Output bulk composition of the equilibrium solid part of the system, moles
    Vector output_b_sw = state_sw.elementAmounts();
    Vector output_n_sw = state_sw.speciesAmounts();
    std::cout << "output_b_sw = \n" << output_b_sw << std:: endl;
    std::cout << "output_n_sw = \n" << output_n_sw << std:: endl;
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
