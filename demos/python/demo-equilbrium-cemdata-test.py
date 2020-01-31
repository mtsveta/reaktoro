from reaktoro import *
import os
import thermofun.PyThermoFun as thermofun
import numpy as np
import sys

T = 293.15 # in kevin
P = 100000.0 # pascal

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
aqueous_species, gaseous_species, minerals = load_data()
print("Aqueous species: ", aqueous_species)
print("\nGaseous species: ", gaseous_species)
print("\nMinerals species: ", minerals)

editor = ChemicalEditor(database)
editor.setTemperatures([T], "kelvin")
editor.setPressures([P], "pascal")


'''

params = DebyeHuckelParams()
params.setPHREEQC()
editor.addAqueousPhase(aqueous_species).setChemicalModelDebyeHuckel(params)

params = DebyeHuckelParams()
params.bneutral(0.122657)
params.aion(3.67)
params.bion(0.122657)
editor.addAqueousPhase(aqueous_species).setChemicalModelDebyeHuckel(params)
'''

editor.addAqueousPhase(aqueous_species).setChemicalModelHKF()
editor.addGaseousPhase(gaseous_species)

# Add mineral species
phase_indices = np.array([2, 2, 6, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
#'''
# Add solid solutions and pure minerals
index_minerals = 0
for i in range(0, len(phase_indices)):
    #print(minerals[index_minerals:index_minerals+phase_indices[i]])
    editor.addMineralPhase(minerals[index_minerals:index_minerals+phase_indices[i]])
    index_minerals = index_minerals + phase_indices[i]
#'''

# Create chemical system
system = ChemicalSystem(editor)
# Print system charachteristics
num_elements = system.numElements()
num_species = system.numSpecies()
print(f"System with with {num_elements} elements, {num_species} species, {system.numPhases()} phases")
# Must coinside with: 'Al' 'C' 'Ca' 'Cl' 'Fe' 'H' 'K' 'Mg' 'Na' 'Nit' 'O' 'S' 'Si' 'Zz'
# Our list of elements: List with 14 elements: Al C Ca Cl Fe H K Mg N Na O S Si Z
# TODO: swap b value for 8th and 9th elements!

'''
print(f"List with {(system.numElements())} elements:", end=" ")
for element in system.elements(): print(element.name(), end=" ")
print("")

engine = thermofun.ThermoEngine(database)
for species in minerals:
    print(f"{species:>15} : {engine.thermoPropertiesSubstance((T), P, species).gibbs_energy.val:>10}")
'''

problem = EquilibriumProblem(system)
problem.setTemperature(T, "kelvin")
problem.setPressure(P, "pascal")
problem.add("H2O", 1, "kg")
problem.add("CaCO3", 10, "g")
problem.add("KOH", 5, "g")
problem.add("SiO2", 10, "g")
problem.add("CaSO4", 5, "g")
problem.add("Al2O3", 5, "g")
problem.add("Fe2O3", 5, "g")
problem.add("H2", 1, "g")

state = equilibrate(problem)
state.output("state.txt")

print("Species in cement : n (in mol)")
for species in system.species():
    name = species.name()
    amount = state.speciesAmount(name)
    # Output according to the threshold
    if amount > 1e-5: print(f"{name:>30} = {amount}")
