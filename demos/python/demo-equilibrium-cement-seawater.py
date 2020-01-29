from reaktoro import *
import os
import thermofun.PyThermoFun as thermofun
import numpy as np
import sys

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
print("Aqueous species: \n ", aqueous_species)
print("Gaseous species: \n ", gaseous_species)
print("Minerals species: \n ", minerals)

editor = ChemicalEditor(database)
editor.setTemperatures([T], "kelvin")
editor.setPressures([P], "pascal")
editor.addAqueousPhase(aqueous_species)
editor.addGaseousPhase(gaseous_species)


# Add mineral species
phase_indices = np.array([2, 2, 6, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

# Add solid solutions and pure minerals
index_minerals = 0
for i in range(0, len(phase_indices)):
    #print(minerals[index_minerals:index_minerals+phase_indices[i]])
    editor.addMineralPhase(minerals[index_minerals:index_minerals+phase_indices[i]])
    index_minerals = index_minerals + phase_indices[i]
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

'''
# Print formula matrix
matrix = system.formulaMatrix()
np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(precision=1)
np.set_printoptions(suppress=True)
print(f"Formula matrix of the size {matrix.shape[0]} x {matrix.shape[1]}:\n", matrix.T)
'''

# Print system charachteristics
num_elements = system.numElements()
num_species = system.numSpecies()
print(f"System with with {num_elements} elements, {num_species} species, {system.numPhases()} phases")

# For system together with solid solutions:
# 14 elements,
# 160 speceis,
# 68 phases.

# For system only with pure minerals:
# 14 elements,
# 142 speceis,
# 61 phases.

# Create EquilibriumSolver for the chemical system
solver = EquilibriumSolver(system)

# Data arrays which are provided neither in DCH nor in DBR files
# B: Full total bulk composition (vector b), moles [nIC] (will be partially re-written from DBR files)
b_cement =np.array([0.00216091995294773, 0.00131648456336122, 0.0265187935777118, 3.45606512042692e-08, 0.00106005014319273, 0.158129994752999, 0.000679052483286759, 0.000733976294881496, 0.00044118313690989, 2.63859873267357e-08, 0.353701411160683, 0.00102638869260499, 0.118139969440399, 0])

# Create chemical state for the cement
state_cement = ChemicalState(system)
solver.solve(state_cement, T, P, b_cement)
state_cement.output("state_cement.txt")

print("Species in cement : n (in mol)")
for species in system.species():
    name = species.name()
    amount = state_cement.speciesAmount(name)
    # Output according to the threshold
    if amount >= 1e-10: print(f"{name:>13} = {amount}")

'''
molar_masses = np.array([0.123046142578125, 0.219110744476318, 0.0269815406799316, 0.0429809408187866, 0.0589803409576416, 0.0599882909059525, 0.0439888907670975, 0.104073191761971, 0.13506404209137, 0.100086999893188, 0.101094949841499, 0.136142601013184, 0.0400779991149902, 0.0570853492021561, 0.117169650197029, 0.116161700248718, 0.115854001998901, 0.116861951947212, 0.152917553067207, 0.151909603118896, 0.0558450012207031, 0.0912980003356934, 0.072852351307869, 0.152917553067207, 0.151909603118896, 0.24797420501709, 0.0558450012207031, 0.145704702615738, 0.235564404010773, 0.0912980003356934, 0.126750999450684, 0.162203998565674, 0.0718444013595581, 0.0878438014984131, 0.0888517514467239, 0.072852351307869, 0.132936652302742, 0.135162902832031, 0.0390983009338379, 0.0561056510210037,
                         0.084314001083374, 0.0853219510316849, 0.0243050003051758, 0.0413123503923416, 0.120369602203369, 0.101396651387215, 0.100388701438904, 0.08299880027771, 0.0840067502260208, 0.119054401397705, 0.0229897994995117, 0.0399971495866775, 0.0770916510820389, 0.272336004257202, 0.060084300994873, 0.076083701133728, 0.0440096006393433, 0.0600090007781982, 0.0610169507265091, 0.0160426001548767, 0.0994505996704102, 0.0354529991149902, 0.0020158998966217, 0.0280133991241455, 0.03199880027771, 0.112132203102112, 0.0810731517076492, 0.0800652017593384, 0.0970725518465042, 0.0960646018981934, 0.0340829012393951, 0.0330749512910843, 0.0320670013427734, 0.0170073500871658, 0.00100794994831085, 0.0180153000354767, 0.0440096006393433, 0.0160426001548767, 0.0020158998966217, 0.0280133991241455,
                         0.03199880027771, 0.0340829012393951, 0.0180153000354767, 0.568448779821396, 0.280665238618851, 0.427353849067688, 0.456217309608459, 0.169212552442741, 0.173886048006892, 0.119821407543999, 0.124494053130627, 0.0436727457165718, 0.0356184949994087, 0.418322781562805, 0.400307481527328, 0.560469779253006, 0.622519680976868, 0.382314694950822, 0.418370296070817, 0.225663452744484, 0.165579151749611, 0.0780035909414291, 0.258160483837128, 1.38665776348114, 0.270193479537964, 0.228316498756409, 0.48595908164978, 0.158038681030273, 0.259999962806702, 0.349230830550194, 0.378285279750824, 0.524439179182053, 0.560469779253006, 0.668561579465866, 0.33819168138504, 0.595496730923653, 0.622519680976868, 0.658550281047821, 0.694580881118774, 0.568473780870438,
                         0.506467685222626, 0.359276133179665, 0.373284531474113, 0.532418179750442, 0.537436329483986, 0.564459279537201, 0.510413379430771, 0.568448779821396, 0.609940379142761, 0.912820183038712, 0.840758982896805, 0.100086999893188, 0.100086999893188, 0.436012200832367, 0.618196700334549, 0.58615560054779, 0.644191000938416, 0.0560773992538452, 0.0740926992893219, 0.136142601013184, 0.172173201084137, 0.145150251030922, 0.115854001998901, 0.115854001998901, 0.231532604217529, 0.0888517514467239, 0.11997900390625, 0.0879120025634766, 0.278016703367233, 0.174261203765869, 0.328419104814529, 0.443331883907318, 0.0583197004795074, 0.142044200897217, 0.380223783969879, 0.425845893621445, 0.548399885177612, 0.0320670013427734, 0.060084300994873, 0.060084300994873])
phases_names = ['aq_gen', 'gas_gen', 'Friedels-ss', 'C3(AF)S0.84H', 'CSHQ', 'straetlingite', 'SO4_OH_AFm', 'SO4_CO3_AFt', 'MSH', 'Al(OH)3mic', 'Kaolinite', 'Mayenite', 'Aluminate', 'Alite', 'Ferrite', 'CA', 'CA2', 'C2AH75', 'C3AH6', 'C4AH11', 'C4AH13', 'C4AH19', 'CAH10', 'C4AsH105', 'C4AsH12', 'C4AsH14', 'C4AsH16', 'C4AsH9', 'Chabazite', 'ZeoliteP', 'C2ASH55', 'C4AcH9', 'C4Ac0.5H105', 'C4Ac0.5H12', 'C4Ac0.5H9', 'C4AcH11', 'Kuzels', 'C6AsH13', 'C6AsH9', 'Aragonite', 'Calcite', 'C3FH6', 'C4FH13', 'C4Fc05H10', 'C4FcH12', 'lime', 'Portlandite', 'Anhydrite', 'Gypsum', 'hemihydrate', 'Fe-carbonate', 'Siderite', 'Magnetite', 'Ferrihydrite-mc', 'Pyrite', 'Troilite', 'Melanterite', 'arcanite', 'syngenite', 'OH-hydrotalcite', 'Brucite', 'thenardite', 'Natrolite', 'ZeoliteX', 'ZeoliteY', 'Sulphur', 'Quartz', 'Silica-amorph']
'''

# Define equilibrium state for the seawater
problem_sw = EquilibriumProblem(system)
b_sw = np.array([1e-12, 8.39047799118904e-07, 4.60274559011631e-06, 0.00024470481263518, 1e-12, 0.0480862507085869, 4.5682055764741e-06, 2.37779781110928e-05, 0.000209584727160819, 5.4730517594269e-08, 0.239750374257753, 1.26338032622899e-05, 0.107827168643251, 0])
state_sw = ChemicalState(system)
solver.solve(state_sw, T, P, b_sw)
'''
print("Species in seawater : n (in mol)")
for species in system.species():
    name = species.name()
    amount = state_sw.speciesAmount(name)
    # Output according to the threshold
    print(f"{name:>13} = {amount}")
state_sw.output("state_sw.txt")
'''
'''
n_sw = np.array([2.46100221484366e-21, 7.49133660286559e-22, 5.75242310773387e-21, 6.90681743697294e-17, 9.90800739368853e-13, 9.12927398011268e-15, 8.56314982070215e-19, 4.72498975561693e-20, 5.94922644003287e-21, 8.54248693442392e-09, 1.28390023290633e-08, 3.45691884021378e-07, 4.23557433459937e-06, 2.05232883296596e-11, 7.63063575495622e-11, 1.0525861918408e-12, 4.39471302392405e-25, 3.63567913737693e-25, 4.80923510908748e-32, 9.91586748163519e-25, 1.39895482718739e-23, 1.90008627164411e-24, 1.53430533643884e-25, 2.37937181783215e-30, 3.39737043871988e-23, 2.2047927969019e-24, 5.7209131449752e-23, 1.19305043639265e-28, 3.75830040054004e-34, 6.33502555926846e-23, 2.9355392326128e-23, 8.30349031726048e-25, 6.72188840987796e-14, 1.85696519699273e-13, 7.47079270448143e-13, 5.20123301714824e-18, 1.2431602548248e-19, 9.06591037658476e-08, 4.47754516670466e-06, 1.3060035882638e-12,
                 2.61236812901651e-08, 5.98734301787103e-08, 2.16483206816277e-05, 2.55807225290953e-09, 2.04024260738186e-06, 8.03111666803934e-10, 5.65266946690013e-11, 4.20224888575075e-08, 7.27152956534961e-08, 2.95494185377378e-06, 0.000206514926295945, 1.21226589402462e-10, 2.82491060715618e-09, 2.28299242396247e-19, 6.4881254166339e-08, 8.93964212306933e-14, 4.14429201643602e-09, 2.30155215456548e-08, 5.89771600313447e-07, 0, 2.89784766978446e-27, 0.00024470481263518, 0, 2.73652587971345e-08, 1.13662409601231e-07, 0, 0, 0, 6.60908355733222e-13, 7.20226715243868e-06, 0, 0, 0, 8.34730158698885e-10, 2.67834682848817e-12, 0.0240427541324881, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1078271, 0])

problem_sw.setTemperature(T, "kelvin")
problem_sw.setPressure(P, "pascal")

for species, n in zip(system.species(), n_sw):
    problem_sw.add(species.name(), n, "mol")
'''
#


'''
problem_naclh = EquilibriumProblem(system)
n_nacl = np.array([7.08390205407382e-23, 3.05792140115698e-30, 9.46429572449386e-16, 3.38222265993579e-14, 1.16939763075074e-12, 2.25272922927293e-13, 7.89337370586244e-15, 9.27837865538291e-21, 2.92542199589053e-28, 2.4941732333023e-21, 7.44386960014572e-20, 1.84774555939163e-20, 1.4373321272283e-12, 3.70236949679406e-19, 2.93246636232349e-23, 2.03906048121015e-26, 1.97089343351131e-28, 3.23458112235902e-27, 7.82524797623901e-35, 8.13304052568418e-29, 7.28469372559797e-21, 1.03067270957778e-21, 4.24726410837069e-24, 3.59392121966752e-32, 2.72617255684425e-26, 2.50890769361049e-34, 2.62393470512793e-19, 6.15265151833959e-24, 2.06779162516089e-28, 3.18992129323225e-19, 1.53978234556082e-19, 4.30494795195437e-21, 9.15847575258294e-13, 6.11763364772575e-15, 5.1403007957578e-13, 1.33655574185217e-15, 6.8053150064198e-22, 4.12700321827989e-21, 1.43733256763097e-12, 2.11469419530701e-20,
                   1.49375478413927e-21, 6.79050879060711e-20, 1.43732347098944e-12, 9.03112256091319e-18, 2.13334759017737e-20, 6.03825917427575e-23, 2.14274044876232e-25, 3.75831596989159e-14, 1.35902158952379e-12, 4.83202429846381e-13, 0.000238186675912991, 7.05055407480896e-12, 2.98958603678772e-15, 0, 1.43433453497021e-12, 4.29224068516673e-21, 1.38112208616446e-12, 1.66762916903346e-14, 9.39850682283554e-12, 0, 9.89464075658034e-30, 0.000238186684843353, 0, 2.67444707581408e-08, 6.97817229238136e-09, 0, 0, 0, 1.92970529701553e-18, 9.54128189344437e-13, 0, 0, 0, 4.1403441039378e-11, 5.52155873495934e-11, 0.0240328647467665, 0, 0, 0, 0,
                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1078271, 0])

problem_nacl.setTemperature(T, "kelvin")
problem_nacl.setPressure(P, "pascal")

for species, n in zip(system.species(), n_nacl):
    problem_nacl.add(species.name(), n, "mol")
b_nacl = np.array([1.43733259290492e-12, 1.21929100962448e-11, 1.43733259290492e-12, 0.000238186684843354, 1.43733259290492e-12, 0.0480657296087116, 1.43733259290492e-12, 1.43733259290492e-12, 0.000238186684843354, 5.34889415162815e-08, 0.239687078800178, 1.43733259290492e-12, 0.107827100001437, 0])

solver = EquilibriumSolver(system)
state_nacl = ChemicalState(system)
solver.solve(state_nacl, T, P, b_nacl)

'''

'''


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