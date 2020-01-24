from reaktoro import *
import os
import thermofun.PyThermoFun as thermofun

T = 293.15
P = 100000.0

def load_data():

    data_folder = 'demos/python/data-files/'
    file_with_species = 'species.txt'
    file_with_species_codes = 'species-codes.txt'
    file_with_elements = 'elements.txt'
    file_with_phases = 'phases.txt'

    def get_data_from_file(filepath):
        items = []
        with open(filepath,'r') as file:
            for line in file:
                for word in line.split():
                    items.append(word.replace("'", "")) #items.append(re.sub(r'[^\w]', '', word))
        return items

    species = get_data_from_file(data_folder + file_with_species)
    species_codes = get_data_from_file(data_folder + file_with_species_codes)
    elements = get_data_from_file(data_folder + file_with_elements)
    phases = get_data_from_file(data_folder + file_with_phases)

    #print(species)
    #print(species_codes)
    #print(elements)
    #print(phases)

    # Fetch all the aqueous species
    aqueous_species = [] # codes 'S', 'T', 'W'
    gaseous_species = [] # codes 'G', 'M', 'I'
    minerals = [] # code 'O'
    for item in zip(species, species_codes):
        if item[1] == 'S' or item[1] == 'T' or item[1] == 'W': aqueous_species.append(item[0])
        if item[1] == 'G' or item[1] == 'M' or item[1] == 'I': gaseous_species.append(item[0])
        if item[1] == 'O': minerals.append(item[0])

    return aqueous_species, gaseous_species, minerals

# Load the database
# Plot all result files
database_path = 'databases/thermofun/'
database_files = sorted(os.listdir(database_path))
print('Loaded databases: ', database_files)

databases = [thermofun.Database(database_path + file) for file in database_files[0:4]]
print(databases)

# Load species from the date files
aqueous_species, gaseous_species, minerals = load_data()
print("Aqueous species: \n ", aqueous_species)
print("Gaseous species: \n ", gaseous_species)
print("Minerals species: \n ",minerals)

editor = ChemicalEditor(databases[3])
editor.setTemperatures(T, "kelvin")
editor.setPressures(P, "pascal")
editor.addAqueousPhase(aqueous_species)

# System with
# 14 elements,
# 160 speceis,
# 68 phases.

editor.addAqueousPhase([
    "Al(OH)2+", "Al(OH)3@", "Al(OH)4-", "Al+3", "AlH3SiO4+2", "AlOH+2",
    "Ca+2", "CaCO3@", "CaCl+", "CaCl2@", "CaHCO3+", "CaHSiO3+", "CaOH+", "CaSiO3@", "K+", "KAlO2@",
    "KCl@", "KOH@", "KCO3-", "KHCO3@", "Mg+2", "MgCO3@", "MgCl+", "MgCl2@", "MgHCO3+", "MgHSiO3+",
    "MgOH+", "MgSiO3@", "Na+", "NaAl(OH)4@", "NaCO3-", "NaCl@", "NaHCO3@", "NaHSiO3@",
    "NaOH@", "HSiO3-", "SiO2@", "CO@", "CO2@", "CO3-2", "HCO3-", "CH4@", "Cl-",
    "HCl@", "H2@", "O2@", "OH-", "H+", "H2O@"])

editor.addMineralPhase("Albite")
editor.addMineralPhase("Andalusite")
editor.addMineralPhase("Calcite")
editor.addMineralPhase("Corundum")
editor.addMineralPhase("Diopside")
editor.addMineralPhase("Dolomite")
editor.addMineralPhase("Enstatite")
editor.addMineralPhase("Grossular")
editor.addMineralPhase("Margarite")
editor.addMineralPhase("Microcline")
editor.addMineralPhase("Muscovite")
editor.addMineralPhase("Pargasite-Mg")
editor.addMineralPhase("Phlogopite")
editor.addMineralPhase("Quartz")
editor.addMineralPhase("Sanidine")
editor.addMineralPhase("Sillimanite")
editor.addMineralPhase("Zoisite")

system = ChemicalSystem(editor)

problem = EquilibriumProblem(system)
problem.add("H2O",             	1000,	"g")
problem.add("CO2",             	0.001,  "g")
problem.add("CaCO3",           	1,	    "g")
problem.add("MgSiO3",          	1,	    "g")
problem.add("NaCl",            	5,      "g")
problem.add("NaAlSi3O8",        37,	    "g")
problem.add("KAl3Si3O10(OH)2",  13,	    "g")
problem.add("SiO2",          	30,	    "g")
problem.add("KAlSi3O8",        	20,	    "g")
problem.setTemperature(T, "Kelvin")
problem.setPressure(P, "Pascal")

state = equilibrate(problem)

state.output("result.txt")

'''
# Step 7.3: Define the initial condition of the reactive transport modeling problem
problem_ic = EquilibriumProblem(system)
problem_ic.setTemperature(T)
problem_ic.setPressure(P)
problem_ic.add('H2O', 1.0, 'kg')
problem_ic.add('NaCl', 0.7, 'mol')
problem_ic.add('CaCO3', 10, 'mol')
problem_ic.add('SiO2', 10, 'mol')

# Step 7.4: Define the boundary condition of the reactive transport modeling problem
problem_bc = EquilibriumProblem(system)
problem_bc.setTemperature(T)
problem_bc.setPressure(P)
problem_bc.add('H2O', 1.0, 'kg')
problem_bc.add('NaCl', 0.90, 'mol')
problem_bc.add('MgCl2', 0.05, 'mol')
problem_bc.add('CaCl2', 0.01, 'mol')
problem_bc.add('CO2', 0.75, 'mol')

options = EquilibriumOptions()
options.smart.reltol = 0.1
options.smart.abstol = 0.1

# Step 7.5: Calculate the equilibrium states for the initial and boundary conditions
state_ic = equilibrate(problem_ic)
state_bc = equilibrate(problem_bc)

# Step 7.6: Scale the volumes of the phases in the initial condition
state_ic.scalePhaseVolume('Aqueous', 0.1, 'm3')
state_ic.scalePhaseVolume('Quartz', 0.882, 'm3')
state_ic.scalePhaseVolume('Calcite', 0.018, 'm3')

# Scale the boundary condition state to 1 m3
state_bc.scaleVolume(1.0, 'm3')

'''