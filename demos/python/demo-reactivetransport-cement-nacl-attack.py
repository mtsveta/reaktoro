# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright (C) 2014-2018 Allan Leal
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library. If not, see <http://www.gnu.org/licenses/>.

# Step 1: Import the reaktoro Python package (and other packages)
from reaktoro import *
import os
import thermofun.PyThermoFun as thermofun
from demos.python.plotting.plotting import *
from natsort import natsorted

#------------------------------------------------------------------------------#
# Problems parameters
#------------------------------------------------------------------------------#

# Step 2: Initialise auxiliary time-related constants
second = 1
minute = 60
hour = 60 * minute
day = 24 * hour
week = 7 * day
year = 365 * day

# Step 3: Define parameters for the reactive transport simulation
xl = 0.0          # the x-coordinate of the left boundary
xr = 1.0          # the x-coordinate of the right boundary
ncells = 100      # the number of cells in the discretization
nsteps = 600      # the number of steps in the reactive transport simulation
D  = 2.3e-9       # the diffusion coefficient (in units of m2/s)
v  = 1.0/week           # the fluid pore velocity (in units of m/s)
dt = 30*minute    # the time step (30 minutes in units of s)
T = 293.15        # in kevin
P = 100000.0      # in pascal

# Premiability 1 Darcy
# Dispersivity 1e-3
# Tortuosity

dx = (xr - xl)/ncells
alpha = v*dt/dx
ndigits = len(str(nsteps))
xcells = np.linspace(xl, xr, ncells)  # the x-coordinates of the plots
results_folder = 'results-nacl'
videos_folder = "videos-nacl"
figures_folder = "figures-nacl"
CFL = dt * v / dx
print("CFL = ", CFL)

# Options for the figure plotting
plot_at_selected_steps = [1, 10, 30, 60, 120, 240, 360, 480, 540, 600]  # the time steps at which the results are plotted

# Step 4: The list of quantities to be output for each mesh cell, each time step
output_quantities = """
    pH
    Eh
    ionicStrength
    elementMolality(C)
    elementMolality(Cl)
    elementMolality(Ca)
    elementMolality(Fe)
    elementMolality(K)
    elementMolality(Mg)    
    elementMolality(Na) 
    elementMolality(S)
    elementMolality(Si) 
    elementMolality(Al)     
    speciesMolality(H+)
    speciesMolality(Cl-)
    speciesMolality(Ca+2)
    speciesMolality(Mg+2)
    speciesMolality(Na+)
    speciesMolality(HCO3-)
    speciesMolality(CO2@)
    phaseVolume(Cal)
    phaseVolume(hydrotalcite)
    phaseVolume(Portlandite)
    phaseVolume(C4AH11)
    phaseVolume(CSHQ-JenD-CSHQ-JenH-CSHQ-TobD-CSHQ-TobH-KSiOH-NaSiOH)
    phaseVolume(C3AFS0.84H4.32-C3FS0.84H4.32)
    phaseVolume(Brc)
    phaseVolume(tricarboalu03-ettringite03_ss)
    elementMolality(H) 
    elementMolality(O)
""".split()

indx_pH = 0
indx_Eh = 1
indx_SI = 2
indx_C = 3
indx_Cl = 4
indx_Ca = 5
indx_Fe = 6
indx_K = 7
indx_Mg = 8
indx_Na = 9
indx_S = 10
indx_Si = 11
indx_Al = 12
indx_Hcation = 13
indx_Clanoion = 14
indx_Cacation = 15
indx_Mgcation = 16
indx_Nacation = 17
indx_HCO3anoion = 18
indx_CO2aq = 19
indx_Cal = 20
indx_hydrotalcite = 21
indx_Portlandite = 22
indx_C4AH11 = 23
indx_CSHQ = 24
indx_C3AFS = 25
indx_Brc = 26
indx_ettringite = 27
indx_H = 28
indx_O = 29

#------------------------------------------------------------------------------#
# Auxiliary functions
#------------------------------------------------------------------------------#

# Step 10: Return a string for the title of a figure in the format Time: #h##m
def titlestr(t):
    t = t / minute   # Convert from seconds to minutes
    h = int(t) / 60  # The number of hours
    m = int(t) % 60  # The number of remaining minutes
    return 'Time: {:>3}h{:>2}m'.format(h, str(m).zfill(2))

# Step 6: Creating folders for the results' output
def make_results_folders():
    os.system('mkdir -p ' + results_folder)
    os.system('mkdir -p ' + figures_folder + '/ph')
    os.system('mkdir -p ' + figures_folder + '/Eh')
    os.system('mkdir -p ' + figures_folder + '/SI')
    os.system('mkdir -p ' + figures_folder + '/elements')
    os.system('mkdir -p ' + figures_folder + '/aqueous_species')
    os.system('mkdir -p ' + figures_folder + '/solids')
    os.system('mkdir -p ' + videos_folder)

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

#------------------------------------------------------------------------------#
# Reactive transport simulations
#------------------------------------------------------------------------------#
def simulate():
    # Load the database
    database_path = 'databases/thermofun/cemdata18-thermofun.json'
    database = thermofun.Database(database_path)

    # Load species from the date files
    aqueous_species = ['Al(SO4)+', 'Al(SO4)2-', 'Al+3', 'AlO+', 'AlO2-', 'AlO2H@', 'AlOH+2', 'AlHSiO3+2', 'AlSiO5-3',
                       'Ca(CO3)@', 'Ca(HCO3)+', 'Ca(SO4)@', 'Ca+2', 'CaOH+', 'Ca(HSiO3)+', 'CaSiO3@', 'Fe(CO3)@',
                       'Fe(HCO3)+', 'Fe(HSO4)+', 'Fe(SO4)@', 'Fe+2', 'FeCl+', 'FeOH+', 'Fe(HSO4)+2', 'Fe(SO4)+',
                       'Fe(SO4)2-', 'Fe+3', 'Fe2(OH)2+4', 'Fe3(OH)4+5', 'FeCl+2', 'FeCl2+', 'FeCl3@', 'FeO+', 'FeO2-',
                       'FeO2H@', 'FeOH+2', 'FeHSiO3+2', 'K(SO4)-', 'K+', 'KOH@', 'Mg(CO3)@', 'Mg(HCO3)+', 'Mg+2',
                       'MgOH+',
                       'MgSO4@', 'Mg(HSiO3)+', 'MgSiO3@', 'Na(CO3)-', 'Na(HCO3)@', 'Na(SO4)-', 'Na+', 'NaOH@', 'HSiO3-',
                       'Si4O10-4', 'SiO2@', 'SiO3-2', 'CO2@', 'CO3-2', 'HCO3-', 'CH4@', 'ClO4-', 'Cl-', 'H2@', 'N2@',
                       'O2@',
                       'S2O3-2', 'HSO3-', 'SO3-2', 'HSO4-', 'SO4-2', 'H2S@', 'HS-', 'S-2', 'OH-', 'H+', 'H2O@']
    gaseous_species = ['CO2(g)', 'CH4(g)', 'H2(g)', 'N2(g)', 'O2(g)', 'H2S(g)', 'H2O(g)']

    # Here, natrolite and Mag are removed from the minerals
    minerals = ['monocarbonate', 'C2AClH5', 'C3AFS0.84H4.32', 'C3FS0.84H4.32', 'CSHQ-JenD', 'CSHQ-JenH', 'CSHQ-TobD',
                'CSHQ-TobH', 'KSiOH', 'NaSiOH', 'straetlingite', 'straetlingite7', 'C4AH13', 'monosulphate12',
                'tricarboalu03', 'ettringite03_ss', 'M075SH', 'M15SH', 'AlOHmic', 'Kln', 'C12A7', 'C3A', 'C3S', 'C4AF',
                'CA', 'CA2', 'C2AH7.5', 'C3AH6', 'C4AH11', 'C4AH13', 'C4AH19', 'CAH10', 'monosulphate10.5',
                'monosulphate12',
                'monosulphate14', 'monosulphate16', 'monosulphate9', 'chabazite', 'zeoliteP_Ca', 'straetlingite5.5',
                'monocarbonate9', 'hemicarbonat10.5', 'hemicarbonate', 'hemicarbonate9', 'monocarbonate', 'C4AsClH12',
                'ettringite13', 'ettringite9', 'Arg', 'Cal', 'C3FH6', 'C4FH13', 'Fe-hemicarbonate', 'Femonocarbonate',
                'Lim',
                'Portlandite', 'Anh', 'Gp', 'hemihydrate', 'FeCO3(pr)', 'Sd', 'FeOOHmic', 'Py', 'Tro', 'Melanterite',
                'K2SO4',
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
    phase_indices = np.array(
        [2, 2, 6, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         1, 1])

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
    # Cement / Initial boundary condition
    # -------------------------------------------------------------------------------------------------------------------- #
    # Bulk composition of reactive subsystem, moles of ICs
    b_cement = np.array([0.00216091995294773, 0.00131648456336122, 0.0265187935777118, 3.45606512042692e-08,
                         0.00106005014319273, 0.158129994752999, 0.000679052483286759, 0.000733976294881496,
                         0, 0.00044118313690989, 0.353701411160683, 0.00102638869260499,
                         0.118139969440399, 0])
    # Define cement problem
    problem_cement = EquilibriumInverseProblem(system)
    problem_cement.setTemperature(T, "kelvin")
    problem_cement.setPressure(P, "pascal")
    problem_cement.setElementInitialAmounts(b_cement)
    problem_cement.fixSpeciesAmount("Qtz", 0.107827, "mol")

    state_cement = equilibrate(problem_cement)
    state_cement.output("state_cement.txt")

    # -------------------------------------------------------------------------------------------------------------------- #
    # NaCl-brine / Left boundary condition
    # -------------------------------------------------------------------------------------------------------------------- #
    # Bulk composition of reactive subsystem, moles of BCs
    b_nacl = np.array([1.43733259290492e-12, 1.21929100962448e-11, 1.43733259290492e-12, 0.000238186684843354,
                       1.43733259290492e-12, 0.0480657296087116, 1.43733259290492e-12, 1.43733259290492e-12,
                       5.34889415162815e-08, 0.000238186684843354, 0.239687078800178, 1.43733259290492e-12,
                       0.107827100001437, 0])

    problem_nacl = EquilibriumInverseProblem(system)
    problem_nacl.setTemperature(T, "kelvin")
    problem_nacl.setPressure(P, "pascal")
    problem_nacl.setElementInitialAmounts(b_nacl)
    problem_nacl.add("H20", 0.001, "g")
    problem_nacl.fixSpeciesAmount("Qtz", 0.107827, "mol")

    state_nacl = equilibrate(problem_nacl)
    state_nacl.output("state_nacl.txt")

    # Step 9: Scale the volumes of the phases in the initial condition
    state_cement.scalePhaseVolume('Aqueous', 0.1, 'm3') # corresponds to the initial porosity of 10%.
    state_cement.scaleVolume(1.0, 'm3')
    #state_cement.scalePhaseVolume('Quartz', 0.882, 'm3')
    #state_cement.scalePhaseVolume('Calcite', 0.018, 'm3')

    # Step 10: Scale the boundary condition state
    state_nacl.scaleVolume(1.0, 'm3')

    # Step 11: Create the mesh for the column
    mesh = Mesh(ncells, xl, xr)

    # Step 12: Create a chemical field object with every cell having state given by state_ic
    field = ChemicalField(mesh.numCells(), state_cement)

    # Step 13: Initialize the reactive transport solver
    rt = ReactiveTransportSolver(system)
    rt.setMesh(mesh)
    rt.setVelocity(v)
    rt.setDiffusionCoeff(D)
    rt.setBoundaryState(state_nacl)
    rt.setTimeStep(dt)
    rt.initialize()

    # Step 14: Set the output of the reactive transport simulation
    output = rt.output()
    for item in output_quantities:
        output.add(item)
    output.filename(results_folder + '/state.txt')  # Set the name of the output files

    make_results_folders()

    # Step 15: Perform given number of reactive tranport steps
    t = 0.0  # current time variable
    step = 0  # current number of steps

    while step <= nsteps:  # step until the number of steps are achieved
        # Print the progress of the simulation
        print("Progress: {}/{} steps, {} min".format(step, nsteps, t/minute))

        # Perform one reactive transport time step
        rt.step(field)

        # Increment time step and number of time steps
        t += dt
        step += 1

def plot():


    print("Collecting files...")
    # Collect files with results corresponding to smart or reference (classical) solver
    files = [file for file in natsorted(os.listdir(results_folder))]

    params = {"plot_at_selected_steps": plot_at_selected_steps,
              "dt": dt,
              "results_folder": results_folder,
              "figures_folder": figures_folder,
              "videos_folder": videos_folder,
              "files": files,
              "indx_ph": indx_pH,
              "indx_Eh": indx_Eh,
              "indx_SI": indx_SI,
              "indx_C": indx_C,
              "indx_Cl": indx_Cl,
              "indx_Ca": indx_Ca,
              "indx_Fe": indx_Fe,
              "indx_K": indx_K,
              "indx_Mg": indx_Mg,
              "indx_Na": indx_Na,
              "indx_S": indx_S,
              "indx_Si": indx_Si,
              "indx_Al": indx_Al,
              "indx_Hcation": indx_Hcation,
              "indx_Clanoion": indx_Clanoion,
              "indx_Cacation": indx_Cacation,
              "indx_Mgcation": indx_Mgcation,
              "indx_Nacation": indx_Nacation,
              "indx_HCO3anoion": indx_HCO3anoion,
              "indx_CO2aq": indx_CO2aq,
              "indx_Cal": indx_Cal,
              "indx_hydrotalcite": indx_hydrotalcite,
              "indx_Portlandite": indx_Portlandite,
              "indx_C4AH11": indx_C4AH11,
              "indx_CSHQ": indx_CSHQ,
              "indx_C3AFS": indx_C3AFS,
              "indx_Brc": indx_Brc,
              "indx_ettringite": indx_ettringite,
              "indx_H": indx_H,
              "indx_O": indx_O,
              "xcells": xcells}

    plot_figures_ph(params)
    plot_figures_Eh(params)
    plot_figures_SI(params)
    plot_figures_elements(params)
    plot_figures_aqueous_species(params)
    plot_figures_solids(params)

    # Animation options
    params["animation_starts_at_frame"] = 0  # the first frame index to be considered
    # params["animation_ends_at_frame"] = 600 # the last frame index to be considered
    # params["animation_num_frames_to_jump"] = 1 # the number of frames to jump between current and next
    params["animation_ends_at_frame"] = 600  # the last frame index to be considered
    params["animation_num_frames_to_jump"] = 1  # the number of frames to jump between current and next
    params["animation_fps"] = 30  # the number of frames per second
    params["animation_interval_wait"] = 200  # the time (in milliseconds) to wait between each frame
    '''
    plot_animation_ph(params)
    plot_animation_Eh(params)
    plot_animation_SI(params)
    plot_animation_elements(params)
    plot_animation_aqueous_species(params)
    plot_animation_solids(params)
    '''
# Step 5: Define the main function
if __name__ == '__main__':

    # Create folders for the results
    make_results_folders()

    # Run the reactive transport simulations
    simulate()

    # Plotting of the figures and videos
    plot()
