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
from numpy import *
import os
from joblib import Parallel, delayed
import matplotlib.pyplot as plt

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
folder = 'results'

CFL = dt * v  / dx
print("CFL = ", CFL)

# Step 4: The list of quantities to be output for each mesh cell, each time step
output_quantities = """
    pH
    speciesMolality(H+)
    speciesMolality(Ca++)
    speciesMolality(Mg++)
    speciesMolality(HCO3-)
    speciesMolality(CO2(aq))
    speciesMolality(Calcite)
    speciesMolality(Dolomite)
    speciesMolality(Quartz)
    
    speciesMolality(hydrotalcite)
    speciesMolality(Portlandite)
    speciesMolality(C4AH11)
    CSHQ-JenD
    CSHQ-JenH
    CSHQ-TobD
    CSHQ-TobH  
""".split()

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
    os.system('mkdir -p figures/aqueous-species')
    os.system('mkdir -p figures/calcite-dolomite')
    os.system('mkdir -p videos')

# Step 8: Plot all result files and generate a video
def plot():
    # Plot all result files
    files = sorted(os.listdir('results'))
    Parallel(n_jobs=16)(delayed(plotfile)(file) for file in files)
    # Create videos for the figures
    ffmpegstr = 'ffmpeg -y -r 30 -i figures/{0}/%0' + str(ndigits) + 'd.png -codec:v mpeg4 -flags:v +qscale -global_quality:v 0 videos/{0}.mp4'
    os.system(ffmpegstr.format('calcite-dolomite'))
    os.system(ffmpegstr.format('aqueous-species'))
    os.system(ffmpegstr.format('ph'))

def plotfile(file):

    step = int((file.split('.')[0]).split('-')[1])

    print('Plotting figure', step, '...')

    t = step * dt
    filearray = loadtxt('results/' + file, skiprows=1)
    data = filearray.T

    ndigits = len(str(nsteps))

    x = linspace(xl, xr, ncells)

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

#------------------------------------------------------------------------------#
# Reactive transport simulations
#------------------------------------------------------------------------------#

# Step 4: Construct the chemical system with its phases and species
database_path = 'databases/thermofun/cemdata18-thermofun.json'
db = Database(database_path)

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


# Step 5: Create the ChemicalSystem object using the configured editor
system = ChemicalSystem(editor)
# Create EquilibriumSolver for the chemical system
solver = EquilibriumSolver(system)

# Step 6: Define the initial condition of the reactive transport modeling problem
# With swapped values 0.00044118313690989 for Na and 0 for N
b_cement = np.array([0.00216091995294773, 0.00131648456336122, 0.0265187935777118, 3.45606512042692e-08,
                     0.00106005014319273, 0.158129994752999, 0.000679052483286759, 0.000733976294881496,
                     0, 0.00044118313690989, 0.353701411160683, 0.00102638869260499,
                     0.118139969440399, 0])
state_cement = ChemicalState(system)
solver.solve(state_cement, T, P, b_cement)
state_cement.output("state_cement.txt")

# Step 7: Define the boundary condition of the reactive transport modeling problem
# Replace by the array with swapped 8th and 9th elements
b_sw = np.array([1e-12, 8.39047799118904e-07, 4.60274559011631e-06, 0.00024470481263518,
                 1e-12, 0.0480862507085869, 4.5682055764741e-06, 2.37779781110928e-05,
                 5.4730517594269e-08, 0.000209584727160819, 0.239750374257753, 1.26338032622899e-05,
                 0.107827168643251, 0])
state_sw = ChemicalState(system)
solver.solve(state_sw, T, P, b_sw)
state_sw.output("state_sw.txt")

input()

# Step 9: Scale the volumes of the phases in the initial condition
state_cement.scalePhaseVolume('Aqueous', 0.1, 'm3') # corresponds to the initial porosity of 10%.
state_cement.scaleVolume(1.0, 'm3')
#state_cement.scaleMin
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
rt.setBoundaryState(state_bc)
rt.setTimeStep(dt)
rt.initialize()

# Step 14: Set the output of the reactive transport simulation
output = rt.output()
output.add("pH")
output.add("speciesMolality(H+)")
output.add("speciesMolality(Ca++)")
output.add("speciesMolality(Mg++)")
output.add("speciesMolality(HCO3-)")
output.add("speciesMolality(CO2(aq))")
output.add("phaseVolume(Calcite)")
output.add("phaseVolume(Dolomite)")
output.filename('results/rtsolver.txt')  # Set the name of the output files

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


plot()
