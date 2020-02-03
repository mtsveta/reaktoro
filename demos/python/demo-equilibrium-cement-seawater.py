from reaktoro import *
import thermofun.PyThermoFun as thermofun
import numpy as np

T = 293.15 # in kevin
P = 100000.0 # pascal

def load_data():

    data_folder = 'demos/python/data-files/'
    file_with_species = 'species.txt'
    file_with_species_codes = 'species-codes.txt'

    def get_data_from_file(filepath):
        items = []
        with open(filepath,'r') as file:
            for line in file:
                for word in line.split():
                    items.append(word.replace("'", "")) #items.append(re.sub(r'[^\w]', '', word))
        return items

    species = get_data_from_file(data_folder + file_with_species)
    species_codes = get_data_from_file(data_folder + file_with_species_codes)

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
print("Gaseous species: ", gaseous_species)
print("Minerals species: ", minerals)

# Create chemical species
editor = ChemicalEditor(database)
editor.setTemperatures([T], "kelvin")
editor.setPressures([P], "pascal")
editor.addAqueousPhase(aqueous_species).setChemicalModelHKF()
editor.addGaseousPhase(gaseous_species)

# Add mineral species
phase_indices = np.array([2, 2, 6, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                          1, 1, 1, 1])
# Add solid solutions (using provided above indices) and pure minerals
index_minerals = 0
for i in range(0, len(phase_indices)):
    #print(minerals[index_minerals:index_minerals+phase_indices[i]])
    editor.addMineralPhase(minerals[index_minerals:index_minerals+phase_indices[i]])
    index_minerals = index_minerals + phase_indices[i]

# Create chemical system and print its characteristics
system = ChemicalSystem(editor)
num_elements = system.numElements()
num_species = system.numSpecies()
print(f"System with with {num_elements} elements, {num_species} species, {system.numPhases()} phases")

# Must coincide with GEMs: 'Al' 'C' 'Ca' 'Cl' 'Fe' 'H' 'K' 'Mg' 'Na' 'Nit' 'O' 'S' 'Si' 'Zz'
# Our list with 14 elements: Al C Ca Cl Fe H K Mg N Na O S Si Z
# Have to swap b value for 8th and 9th elements!
print(f"List with {(system.numElements())} elements:", end=" ")
for element in system.elements(): print(element.name(), end=" ")
print("")

# -------------------------------------------------------------------------------------------------------------------- #
# Cement
# -------------------------------------------------------------------------------------------------------------------- #
# Our list of elements: List with 14 elements:
# Al C Ca Cl
# Fe H K Mg
# N Na O S
# Si Z
# With swapped values 0.00044118313690989 for Na and 0 for N
b_cement = np.array([0.00216091995294773, 0.00131648456336122, 0.0265187935777118, 3.45606512042692e-08,
                     0.00106005014319273, 0.158129994752999, 0.000679052483286759, 0.000733976294881496,
                     0, 0.00044118313690989, 0.353701411160683, 0.00102638869260499,
                     0.118139969440399, 0])
# Create element dictionary and output the names of species and corresponding amounts
b_cement_dict = {}
i = 0
for elements, amount in zip(system.elements(), b_cement):
    b_cement_dict[elements.name()] = amount
    print(f"{elements.name():>2} : {amount}")

problem_cement = EquilibriumInverseProblem(system)
problem_cement.setTemperature(T, "kelvin")
problem_cement.setPressure(P, "pascal")
problem_cement.setElementInitialAmounts(b_cement)
problem_cement.fixSpeciesAmount("Qtz", 0.107827, "mol")
#problem_cement.fixSpeciesAmount("C3AFS0.84H4.32", 0.0, "mol")

state_cement = equilibrate(problem_cement)
state_cement.output("state_cement.txt")

print("Species in cement : n (in mol)")
for species in system.species():
    name = species.name()
    amount = state_cement.speciesAmount(name)
    # Output according to the threshold
    print(f"{name:>30} = {amount}")

# -------------------------------------------------------------------------------------------------------------------- #
# Sea water (SW)
# -------------------------------------------------------------------------------------------------------------------- #

# bIC: Bulk composition of reactive subsystem (main GEM input), moles of ICs [nICb]
# Our list of elements: List with 14 elements:
# Al C Ca Cl
# Fe H K Mg
# N Na O S
# Si Z
# Replace by the array with swapped 8th and 9th elements
b_sw = np.array([1e-12, 8.39047799118904e-07, 4.60274559011631e-06, 0.00024470481263518,
                 1e-12, 0.0480862507085869 + 2, 4.5682055764741e-06, 2.37779781110928e-05,
                 5.4730517594269e-08, 0.000209584727160819, 0.239750374257753 + 1, 1.26338032622899e-05,
                 0.107827168643251, 0])

problem_sw = EquilibriumInverseProblem(system)
problem_sw.setTemperature(T, "kelvin")
problem_sw.setPressure(P, "pascal")
problem_sw.setElementInitialAmounts(b_sw)
problem_sw.fixSpeciesAmount("Qtz", 0.107827, "mol")

state_sw = equilibrate(problem_sw)
state_sw.output("state_sw_equilibrate.txt")

print("Species in cement : n (in mol)")
for species in system.species():
    name = species.name()
    amount = state_sw.speciesAmount(name)
    # Output according to the threshold
    print(f"{name:>30} = {amount}")
