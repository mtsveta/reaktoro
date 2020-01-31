from reaktoro import *
import os
import thermofun.PyThermoFun as thermofun
import numpy as np
import sys

T = 293.15 # in kevin
P = 100000.0 # pascal

T_celsius = 19.85 # in celsius
P_bar = 1 # in bar

def load_data():

    data_folder = 'demos/python/data-files/'
    file_with_species = 'species-sw.txt'
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
    gaseous_species = [] # codes 'G'
    minerals = [] # code 'O', 'M', 'I', J
    for item in zip(species, species_codes):
        # 'S' = solute, 'T' = ?, 'W' = water
        if item[1] == 'S' or item[1] == 'T' or item[1] == 'W': aqueous_species.append(item[0])
        # 'G' = gaseous
        if item[1] == 'G': gaseous_species.append(item[0])
        # 'G' = gaseous, 'I' = ?, 'J' = ?, 'O' = ?,
        if item[1] == 'M' or item[1] == 'I' or item[1] == 'J' or item[1] == 'O': minerals.append(item[0])

    return aqueous_species, gaseous_species, minerals

# Load the database
database_path = 'databases/thermofun/cemdata18-thermofun.json'
database = thermofun.Database(database_path)

# Load species from the date files
aqueous_species = ['Al(SO4)+', 'Al(SO4)2-', 'Al+3', 'AlO+', 'AlO2-', 'AlO2H@', 'AlOH+2', 'AlHSiO3+2', 'AlSiO5-3', 'Ca(CO3)@', 'Ca(HCO3)+', 'Ca(SO4)@', 'Ca+2', 'CaOH+', 'Ca(HSiO3)+', 'CaSiO3@', 'Fe(CO3)@', 'Fe(HCO3)+', 'Fe(HSO4)+', 'Fe(SO4)@', 'Fe+2', 'FeCl+', 'FeOH+', 'Fe(HSO4)+2', 'Fe(SO4)+', 'Fe(SO4)2-', 'Fe+3', 'Fe2(OH)2+4', 'Fe3(OH)4+5', 'FeCl+2', 'FeCl2+', 'FeCl3@', 'FeO+', 'FeO2-', 'FeO2H@', 'FeOH+2', 'FeHSiO3+2', 'K(SO4)-', 'K+', 'KOH@', 'Mg(CO3)@', 'Mg(HCO3)+', 'Mg+2', 'MgOH+', 'MgSO4@', 'Mg(HSiO3)+', 'MgSiO3@', 'Na(CO3)-', 'Na(HCO3)@', 'Na(SO4)-', 'Na+', 'NaOH@', 'HSiO3-', 'Si4O10-4', 'SiO2@', 'SiO3-2', 'CO2@', 'CO3-2', 'HCO3-', 'CH4@', 'ClO4-', 'Cl-', 'H2@', 'N2@', 'O2@', 'S2O3-2', 'HSO3-', 'SO3-2', 'HSO4-', 'SO4-2', 'H2S@', 'HS-', 'S-2', 'OH-', 'H+', 'H2O@']
gaseous_species = ['CO2(g)', 'CH4(g)', 'H2(g)', 'N2(g)', 'O2(g)', 'H2S(g)', 'H2O(g)']
minerals = ['monocarbonate', 'C2AClH5', 'C3AFS0.84H4.32', 'C3FS0.84H4.32', 'CSHQ-JenD', 'CSHQ-JenH', 'CSHQ-TobD', 'CSHQ-TobH', 'KSiOH', 'NaSiOH', 'straetlingite', 'straetlingite7', 'C4AH13', 'monosulphate12', 'tricarboalu03', 'ettringite03_ss', 'M075SH', 'M15SH', 'AlOHmic', 'Kln', 'C12A7', 'C3A', 'C3S', 'C4AF', 'CA', 'CA2', 'C2AH7.5', 'C3AH6', 'C4AH11', 'C4AH13', 'C4AH19', 'CAH10', 'monosulphate10.5', 'monosulphate12', 'monosulphate14', 'monosulphate16', 'monosulphate9', 'chabazite', 'zeoliteP_Ca', 'straetlingite5.5', 'monocarbonate9', 'hemicarbonat10.5', 'hemicarbonate', 'hemicarbonate9', 'monocarbonate', 'C4AsClH12', 'ettringite13', 'ettringite9', 'Arg', 'Cal', 'C3FH6', 'C4FH13', 'Fe-hemicarbonate', 'Femonocarbonate', 'Lim', 'Portlandite', 'Anh', 'Gp', 'hemihydrate', 'FeCO3(pr)', 'Sd', 'FeOOHmic', 'Py', 'Tro', 'Melanterite', 'K2SO4', 'syngenite', 'hydrotalcite', 'Brc', 'Na2SO4', 'zeoliteX', 'zeoliteY', 'Sulfur', 'Qtz', 'Amor-Sl']
print("Aqueous species: ", aqueous_species)
print("\nGaseous species: ", gaseous_species)
print("\nMinerals species: ", minerals)

editor = ChemicalEditor(database)
editor.setTemperatures([T], "kelvin")
editor.setPressures([P], "pascal")
editor.addAqueousPhase(aqueous_species).setChemicalModelHKF()
'''
params = DebyeHuckelParams()
params.setPHREEQC()
editor.addAqueousPhase(aqueous_species).setChemicalModelDebyeHuckel(params)
'''
'''
params = DebyeHuckelParams()
params.bneutral(0.122657)
params.aion(3.67)
params.bion(0.122657)
editor.addAqueousPhase(aqueous_species).setChemicalModelDebyeHuckel(params)
'''
'''
params = DebyeHuckelParams()
params.aion(3.67)
params.bion(0.122657)
params.bneutral(0.122657)
aqueous_phase.setChemicalModelDebyeHuckel(params)
'''
editor.addGaseousPhase(gaseous_species)

# Add mineral species
phase_indices = np.array([2, 2, 6, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
#'''
# Add solid solutions and pure minerals
index_minerals = 0
for i in range(0, len(phase_indices)):
    #print(minerals[index_minerals:index_minerals+phase_indices[i]])
    editor.addMineralPhase(minerals[index_minerals:index_minerals+phase_indices[i]])
    index_minerals = index_minerals + phase_indices[i]
#'''
'''
# Adding only pure minerals
index_minerals = 0
for i in range(0, len(phase_indices)):
    if phase_indices[i] == 1:
        print(minerals[index_minerals:index_minerals+phase_indices[i]])
        editor.addMineralPhase(minerals[index_minerals:index_minerals+phase_indices[i]])
    index_minerals = index_minerals + phase_indices[i]
'''
# Create chemical system
system = ChemicalSystem(editor)
# Print system charachteristics
num_elements = system.numElements()
num_species = system.numSpecies()
print(f"System with with {num_elements} elements, {num_species} species, {system.numPhases()} phases")
# Must coinside with: 'Al' 'C' 'Ca' 'Cl' 'Fe' 'H' 'K' 'Mg' 'Na' 'Nit' 'O' 'S' 'Si' 'Zz'
# Our list of elements: List with 14 elements: Al C Ca Cl Fe H K Mg N Na O S Si Z
# TODO: swap b value for 8th and 9th elements!

print(f"List with {(system.numElements())} elements:", end=" ")
for element in system.elements(): print(element.name(), end=" ")
print("")

# Create EquilibriumSolver for the chemical system
solver = EquilibriumSolver(system)

# -------------------------------------------------------------------------------------------------------------------- #
# Sea water (SW)
# -------------------------------------------------------------------------------------------------------------------- #

# bIC: Bulk composition of reactive subsystem (main GEM input), moles of ICs [nICb]
"""
b_sw = np.array([1e-12, 8.39047799118904e-07, 4.60274559011631e-06, 0.00024470481263518,
                 1e-12, 0.0480862507085869, 4.5682055764741e-06, 2.37779781110928e-05,
                 0.000209584727160819, 5.4730517594269e-08, 0.239750374257753, 1.26338032622899e-05,
                 0.107827168643251, 0])
"""
# Our list of elements: List with 14 elements:
# Al C Ca Cl
# Fe H K Mg
# N Na O S
# Si Z
# Replace by the array with swapped 8th and 9th elements
b_sw = np.array([1e-12, 8.39047799118904e-07, 4.60274559011631e-06, 0.00024470481263518,
                 1e-12, 0.0480862507085869 + 20, 4.5682055764741e-06, 2.37779781110928e-05,
                 5.4730517594269e-08, 0.000209584727160819, 0.239750374257753 + 10, 1.26338032622899e-05,
                 0.107827168643251, 0])

state_sw = ChemicalState(system)
solver.solve(state_sw, T, P, b_sw)
state_sw.output("state_sw_equilibrium_solver.txt")

print("Species in seawater : n (in mol)")
for species in system.species():
    name = species.name()
    amount = state_sw.speciesAmount(name)
    # Output according to the threshold
    print(f"{name:>30} = {amount}")


problem_sw = EquilibriumInverseProblem(system)
problem_sw.setTemperature(T, "kelvin")
problem_sw.setPressure(P, "pascal")
problem_sw.setElementInitialAmounts(b_sw)
problem_sw.fixSpeciesAmount("Qtz", 0.107827, "mol")

state_sw = equilibrate(problem_sw)
state_sw.output("state_sw_equilibrate.txt")