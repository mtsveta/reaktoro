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
import numpy as np
import os
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
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
ncells = 50      # the number of cells in the discretization
nsteps = 600      # the number of steps in the reactive transport simulation
#nsteps = 3000      # the number of steps in the reactive transport simulation
D  = 2.3e-9       # the diffusion coefficient (in units of m2/s)
v  = 1.0/day           # the fluid pore velocity (in units of m/s)
dt = 10*minute    # the time step (30 minutes in units of s)
T = 293.15        # in kevin
P = 100000.0      # in pascal

# Premiability 1 Darcy
# Dispersivity 1e-3
# Tortuosity

dx = (xr - xl)/ncells
alpha = v*dt/dx
ndigits = len(str(nsteps))
xcells = np.linspace(xl, xr, ncells)  # the x-coordinates of the plots
folder = 'results'

CFL = dt * v  / dx
print("CFL = ", CFL)

# Options for the figure plotting
plot_at_selected_steps = [1, 10, 60, 120, 240, 480, 600]  # the time steps at which the results are plotted

# Step 4: The list of quantities to be output for each mesh cell, each time step
output_quantities = """
    pH
    Eh
    ionicStrength
    elementAmount(C)
    elementAmount(Cl)
    elementAmount(Ca)
    elementAmount(Fe)
    elementAmount(K)
    elementAmount(Mg)    
    elementAmount(Na) 
    elementAmount(S) 
    elementAmount(Si) 
    elementAmount(Al)    
    speciesMolality(H+)
    speciesMolality(Cl-)
    speciesMolality(Ca+2)
    speciesMolality(Mg+2)
    speciesMolality(Na+)
    speciesMolality(HCO3-)
    speciesMolality(CO2@)
    speciesMolality(Cal)
    speciesMolality(hydrotalcite)
    speciesMolality(Portlandite)
    speciesMolality(C4AH11)
    speciesMolality(CSHQ-JenD)
    speciesMolality(CSHQ-JenH)
    speciesMolality(CSHQ-TobD)
    speciesMolality(CSHQ-TobH)
    speciesMolality(C3AFS0.84H4.32)
    speciesMolality(Brc)
    speciesMolality(ettringite03_ss)
    speciesMolality(ettringite13)
    speciesMolality(ettringite9)
""".split()

indx_pH = 0
indx_Eh = 1
indx_IS = 2
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
indx_Cl = 14
indx_Cacation = 15
indx_Mgcation = 16
indx_Nacation = 17
indx_HCO3anoion = 18
indx_CO2aq = 19
indx_Cal = 20
indx_hydrotalcite = 21
indx_Portlandite = 22
indx_C4AH11 = 23
indx_CSHQJenD = 24
indx_CSHQJenH = 25
indx_CSHQTobD = 26
indx_CSHQTobH = 27
indx_C3AFS = 28
indx_Brc = 29
indx_ettringite03_ss = 30
indx_ettringite13 = 31
indx_ettringite9 = 32
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
    os.system('mkdir -p results')
    os.system('mkdir -p figures/ph')
    os.system('mkdir -p videos')

def plotfile(file):

    step = int((file.split('.')[0]).split('-')[1])

    print('Plotting figure', step, '...')

    t = step * dt
    filearray = np.loadtxt('results/' + file, skiprows=1)
    data = filearray.T

    ndigits = len(str(nsteps))

    x = np.linspace(xl, xr, ncells)

    plt.figure()
    plt.xlim(left=-0.02, right=0.52)
    plt.ylim(bottom=2.5, top=10.5)
    plt.title(titlestr(t))
    plt.xlabel('Distance [m]')
    plt.ylabel('pH')
    plt.plot(x, data[0])
    plt.tight_layout()
    plt.savefig('figures/ph/{}.png'.format(str(step).zfill(ndigits)))

    plt.figure()
    plt.xlim(left=-0.02, right=0.52)
    plt.ylim(bottom=-0.1, top=2.1)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.title(titlestr(t))
    plt.xlabel('Distance [m]')
    plt.ylabel('Mineral Volume [%$_{\mathsf{vol}}$]')
    plt.plot(x, data[6] * 100, label='Calcite')
    plt.plot(x, data[7] * 100, label='Dolomite')
    plt.legend(loc='center right')
    plt.tight_layout()
    plt.savefig('figures/calcite-dolomite/{}.png'.format(str(step).zfill(ndigits)))

    plt.figure()
    plt.yscale('log')
    plt.xlim(left=-0.02, right=0.52)
    plt.ylim(bottom=0.5e-5, top=2)
    plt.title(titlestr(t))
    plt.xlabel('Distance [m]')
    plt.ylabel('Concentration [molal]')
    plt.plot(x, data[2], label='Ca++')
    plt.plot(x, data[3], label='Mg++')
    plt.plot(x, data[4], label='HCO3-')
    plt.plot(x, data[5], label='CO2(aq)')
    plt.plot(x, data[1], label='H+')
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig('figures/aqueous-species/{}.png'.format(str(step).zfill(ndigits)))

    plt.close('all')

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

# Load the database
database_path = 'databases/thermofun/cemdata18-thermofun.json'
database = thermofun.Database(database_path)

# Load species from the date files
aqueous_species, gaseous_species, minerals = load_data()
print("Aqueous species: ", aqueous_species)
print("\nGaseous species: ", gaseous_species)
print("\nMinerals species: ", minerals)

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
# Sea water (SW) / Left boundary condition
# -------------------------------------------------------------------------------------------------------------------- #
# Bulk composition of reactive subsystem, moles of BCs
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
state_sw.output("state_sw.txt")

# Step 9: Scale the volumes of the phases in the initial condition
state_cement.scalePhaseVolume('Aqueous', 0.1, 'm3') # corresponds to the initial porosity of 10%.
state_cement.scaleVolume(1.0, 'm3')
#state_cement.scalePhaseVolume('Quartz', 0.882, 'm3')
#state_cement.scalePhaseVolume('Calcite', 0.018, 'm3')

# Step 10: Scale the boundary condition state
state_sw.scaleVolume(1.0, 'm3')

# Step 11: Create the mesh for the column
mesh = Mesh(ncells, xl, xr)

# Step 12: Create a chemical field object with every cell having state given by state_ic
field = ChemicalField(mesh.numCells(), state_cement)

# Step 13: Initialize the reactive transport solver
rt = ReactiveTransportSolver(system)
rt.setMesh(mesh)
rt.setVelocity(v)
rt.setDiffusionCoeff(D)
rt.setBoundaryState(state_sw)
rt.setTimeStep(dt)
rt.initialize()

# Step 14: Set the output of the reactive transport simulation
output = rt.output()
for item in output_quantities:
    output.add(item)
output.filename('results/state.txt')  # Set the name of the output files

make_results_folders()

# Step 15: Perform given number of reactive tranport steps
t = 0.0  # current time variable
step = 0  # current number of steps
"""
while step <= nsteps:  # step until the number of steps are achieved
    # Print the progress of the simulation
    print("Progress: {}/{} steps, {} min".format(step, nsteps, t/minute))

    # Perform one reactive transport time step
    rt.step(field)

    # Increment time step and number of time steps
    t += dt
    step += 1
"""
print("Collecting files...")
# Collect files with results corresponding to smart or reference (classical) solver
files = [file for file in natsorted( os.listdir(folder) )]

params = {"plot_at_selected_steps": plot_at_selected_steps,
          "dt":  dt,
          "folder": folder,
          "files": files,
          "indx_ph": indx_pH,
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
          "xcells": xcells}

plot_figures_ph(params)
plot_figures_elements(params)

# Animation options
params["animation_starts_at_frame"] = 0 # the first frame index to be considered
#params["animation_ends_at_frame"] = 600 # the last frame index to be considered
#params["animation_num_frames_to_jump"] = 1 # the number of frames to jump between current and next
params["animation_ends_at_frame"] = 10 * 30 # the last frame index to be considered
params["animation_num_frames_to_jump"] = 1 # the number of frames to jump between current and next
params["animation_fps"] = 30  # the number of frames per second
params["animation_interval_wait"] = 200  # the time (in milliseconds) to wait between each frame

#plot_animation_ph(params)
