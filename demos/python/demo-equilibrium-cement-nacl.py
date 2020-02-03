from reaktoro import *
import thermofun.PyThermoFun as thermofun
import numpy as np

T = 293.15 # in kevin
P = 100000.0 # pascal

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
aqueous_species = ['Al(SO4)+', 'Al(SO4)2-', 'Al+3', 'AlO+', 'AlO2-', 'AlO2H@', 'AlOH+2', 'AlHSiO3+2', 'AlSiO5-3',
                   'Ca(CO3)@', 'Ca(HCO3)+', 'Ca(SO4)@', 'Ca+2', 'CaOH+', 'Ca(HSiO3)+', 'CaSiO3@', 'Fe(CO3)@',
                   'Fe(HCO3)+', 'Fe(HSO4)+', 'Fe(SO4)@', 'Fe+2', 'FeCl+', 'FeOH+', 'Fe(HSO4)+2', 'Fe(SO4)+',
                   'Fe(SO4)2-', 'Fe+3', 'Fe2(OH)2+4', 'Fe3(OH)4+5', 'FeCl+2', 'FeCl2+', 'FeCl3@', 'FeO+', 'FeO2-',
                   'FeO2H@', 'FeOH+2', 'FeHSiO3+2', 'K(SO4)-', 'K+', 'KOH@', 'Mg(CO3)@', 'Mg(HCO3)+', 'Mg+2', 'MgOH+',
                   'MgSO4@', 'Mg(HSiO3)+', 'MgSiO3@', 'Na(CO3)-', 'Na(HCO3)@', 'Na(SO4)-', 'Na+', 'NaOH@', 'HSiO3-',
                   'Si4O10-4', 'SiO2@', 'SiO3-2', 'CO2@', 'CO3-2', 'HCO3-', 'CH4@', 'ClO4-', 'Cl-', 'H2@', 'N2@', 'O2@',
                   'S2O3-2', 'HSO3-', 'SO3-2', 'HSO4-', 'SO4-2', 'H2S@', 'HS-', 'S-2', 'OH-', 'H+', 'H2O@']
gaseous_species = ['CO2(g)', 'CH4(g)', 'H2(g)', 'N2(g)', 'O2(g)', 'H2S(g)', 'H2O(g)']

# Here, natrolite and Mag are removed from the minerals
minerals = ['monocarbonate', 'C2AClH5', 'C3AFS0.84H4.32', 'C3FS0.84H4.32', 'CSHQ-JenD', 'CSHQ-JenH', 'CSHQ-TobD',
            'CSHQ-TobH', 'KSiOH', 'NaSiOH', 'straetlingite', 'straetlingite7', 'C4AH13', 'monosulphate12',
            'tricarboalu03', 'ettringite03_ss', 'M075SH', 'M15SH', 'AlOHmic', 'Kln', 'C12A7', 'C3A', 'C3S', 'C4AF',
            'CA', 'CA2', 'C2AH7.5', 'C3AH6', 'C4AH11', 'C4AH13', 'C4AH19', 'CAH10', 'monosulphate10.5', 'monosulphate12',
            'monosulphate14', 'monosulphate16', 'monosulphate9', 'chabazite', 'zeoliteP_Ca', 'straetlingite5.5',
            'monocarbonate9', 'hemicarbonat10.5', 'hemicarbonate', 'hemicarbonate9', 'monocarbonate', 'C4AsClH12',
            'ettringite13', 'ettringite9', 'Arg', 'Cal', 'C3FH6', 'C4FH13', 'Fe-hemicarbonate', 'Femonocarbonate', 'Lim',
            'Portlandite', 'Anh', 'Gp', 'hemihydrate', 'FeCO3(pr)', 'Sd', 'FeOOHmic', 'Py', 'Tro', 'Melanterite', 'K2SO4',
            'syngenite', 'hydrotalcite', 'Brc', 'Na2SO4', 'zeoliteX', 'zeoliteY', 'Sulfur', 'Qtz', 'Amor-Sl']
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
                          1, 1])
# Add solid solutions and pure minerals
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

# Must coincide with: 'Al' 'C' 'Ca' 'Cl' 'Fe' 'H' 'K' 'Mg' 'Na' 'Nit' 'O' 'S' 'Si' 'Zz'
# Our list of elements: List with 14 elements: Al C Ca Cl Fe H K Mg N Na O S Si Z
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
state_cement.output("state_cement_without_magnetite_natrolite.txt")

print("Species in cement : n (in mol)")
for species in system.species():
    name = species.name()
    amount = state_cement.speciesAmount(name)
    # Output according to the threshold
    print(f"{name:>30} = {amount}")

# -------------------------------------------------------------------------------------------------------------------- #
# NaCl-brine
# -------------------------------------------------------------------------------------------------------------------- #

# Our list of elements: List with 14 elements:
# Al C Ca Cl
# Fe H K Mg
# N Na O S
# Si Z
# Replace by the array with swapped 8th and 9th elements
# Initial value
b_nacl = np.array([1.43733259290492e-12, 1.21929100962448e-11, 1.43733259290492e-12, 0.000238186684843354,
                   1.43733259290492e-12, 0.0480657296087116, 1.43733259290492e-12, 1.43733259290492e-12,
                   5.34889415162815e-08, 0.000238186684843354, 0.239687078800178, 1.43733259290492e-12,
                   0.107827100001437, 0])

# With 1 kg of water
'''
b_nacl = np.array([1.0e-12, 6.736879e-12, 1.0e-12, 0.5,
                   1.0e-12, 111.01675, 1.e-12, 1.e-12,
                   1.e-12, 0.5, 55.724028, 1e-12,
                   0.107827100001437, 0])
'''
problem_nacl = EquilibriumInverseProblem(system)
problem_nacl.setTemperature(T, "kelvin")
problem_nacl.setPressure(P, "pascal")
problem_nacl.setElementInitialAmounts(b_nacl)
problem_nacl.fixSpeciesAmount("Qtz", 0.107827, "mol")

state_nacl = equilibrate(problem_nacl)
state_nacl.output("state_nacl_without_magnetite_natrolite.txt")

print("Species in NaCl-brine : n (in mol)")
for species in system.species():
    name = species.name()
    amount = state_nacl.speciesAmount(name)
    # Output according to the threshold
    print(f"{name:>30} = {amount}")

