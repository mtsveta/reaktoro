# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: notebooks//ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Reactive transport in granite simulation
#
# In this tutorial, we show how Reaktoro can be used for one-dimensional reactive transport calculations for modeling
# the geochemical reactions that occur along a porous rock column with granite as an acidic brine is continuously
# injected on its left side.
#
# ## Import the reaktoro Python package (and other packages)
# First, we import the **reaktoro** Python package so that we can use its classes
# and methods for performing the chemical reaction calculations.

from reaktoro import *
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from natsort import natsorted
import time
from progress.bar import IncrementalBar

mpl.set_loglevel("critical")
mpl.rc_file("matplotlibrc")
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'TeX Gyre Adventor'
mpl.rcParams['font.style'] = 'normal'
mpl.rcParams['font.size'] = 12
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['axes.unicode_minus'] = False

# ## Defining auxiliary time-related constants
# In this step, we initialize auxiliary time-related constants from seconds to years. This is only done for
# convenience, so that we can specify later, for example, fluid velocity as 1 m/week.

second = 1
minute = 60
hour = 60 * minute
day = 24 * hour
week = 7 * day
year = 365 * day

# ## Defining parameters for the reactive transport simulation
# Next, we define reactive transport and numerical discretization parameters.

# +
# Discretization parameters
xl = 0.0                # the x-coordinate of the left boundary
xr = 1.0                # the x-coordinate of the right boundary
ncells = 100            # the number of cells in the discretization
nsteps = 1000           # the number of steps in the reactive transport simulation
dt = 30*minute          # the time step (30 minutes in units of s)
dx = (xr - xl)/ncells   # length of the mesh cells (in units of m)

# Physical and chemical parameters
D  = 1.0e-9             # the diffusion coefficient (in units of m2/s)
v  = 1.0/week           # the fluid pore velocity (1 m/week in units of m/s)
T = 300.0               # the temperature (in units of degC)
P = 85.88               # the pressure (in units of bar)
phi = 0.1               # the porosity
# -

# Activity model for the aqueous species:

# +
#activity_model = "hkf-full"
activity_model = "hkf-selected-species"
#activity_model = "pitzer-full"
#activity_model = "pitzer-selected-species"
#activity_model = "dk-full"
#activity_model = "dk-selected-species"
# -

tag = "-dt-" + "{:d}".format(dt) + \
      "-ncells-" + str(ncells) + \
      "-nsteps-" + str(nsteps) + \
      "-" + activity_model
folder_results = 'results-rtsolver-granite'
os.system('mkdir -p ' + folder_results)
folder_result_plots = "plots-" + folder_results
os.system('mkdir -p ' + folder_result_plots)


# The seconds spent on equilibrium and transport calculations per time step
time_steps = np.linspace(0, nsteps, nsteps)

# Next, we generate the coordinates of the mesh nodes (array `x`) by equally dividing the interval *[xr, xl]* with
# the number of cells `ncells`. The length between each consecutive mesh nodes is computed and stored in `dx` (the
# length of the mesh cells).

xcells = np.linspace(xl, xr, ncells)    # interval [xl, xr] split into `ncells`

# To make sure that the applied finite-volume scheme is stable, we need to keep track of Courant–Friedrichs–Lewy (CFL)
# number, which should be less than 1.0.

CFL = v*dt/dx
assert CFL <= 1.0, f"Make sure that CFL = {CFL} is less that 1.0"

# ## Reactive transport simulations
#
# ### Defining the chemical system
#
# We need to define a chemical system that can represent both our fluid and rock. We use class
# [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html) below to define a system with an
# aqueous phase and three mineral phases: quartz, calcite, and dolomite. Initially, our rock has no dolomite
# (CaMg(CO<sub>3</sub>)<sub>2</sub>), but since this is a mineral that could potentially precipitate given the fluid
# composition injected ( containing CaCl<sub>2</sub> and MgCl<sub>2</sub> dissolved), we add it here in the
# chemical system to ensure that the calculations are able to model dolomite precipitation.

editor = ChemicalEditor()

selected_elements = "Al Cl H K Na O Si"
selected_species = "H2O(l) H+ OH- Cl- HCl(aq) Na+ NaOH(aq) NaHSiO3(aq) NaCl(aq) " \
                    "K+ KOH(aq) KCl(aq) Al+++ AlOH++"

if activity_model == "pitzer-full":
    editor.addAqueousPhaseWithElements(selected_elements) \
        .setChemicalModelPitzerHMW() \
        .setActivityModelDrummondCO2()
elif activity_model == "hkf-full":
    editor.addAqueousPhaseWithElements(selected_elements)
elif activity_model == "dk-full":
    editor.addAqueousPhaseWithElements(selected_elements) \
        .setChemicalModelDebyeHuckel()
elif activity_model == "pitzer-selected-species":
    editor.addAqueousPhase(selected_species) \
        .setChemicalModelPitzerHMW() \
        .setActivityModelDrummondCO2()
elif activity_model == "hkf-selected-species":
    editor.addAqueousPhase(selected_species)
elif activity_model == "dk-selected-species":
    editor.addAqueousPhase(selected_species) \
        .setChemicalModelDebyeHuckel()

editor.addMineralPhase("Quartz") # SiO2
editor.addMineralPhase("Diaspore") # AlO(OH)
editor.addMineralPhase("Gibbsite") # Al(OH)3
editor.addMineralPhase("Andalusite") # Al2SiO5
editor.addMineralPhase("Kyanite") # Al2SiO5
editor.addMineralPhase("Sillimanite") # Al2SiO5
editor.addMineralPhase("Muscovite") # KAl2(AlSi3)O10(OH)2
editor.addMineralPhase("Paragonite") # NaAl2(AlSi3)O10(OH)2
editor.addMineralPhase("Pyrophyllite") # Al2Si4O10(OH)2
editor.addMineralPhase("Kaolinite") # Al2Si2O5(OH)4
editor.addMineralPhase("Albite") # Na(AlSi3)O8
editor.addMineralPhase("K-Feldspar") # K(AlSi3)O8

# > **Note**: The aqueous phase is defined above by using a list of compounds, which is then broken automatically by
# > Reaktoro into a list of element names. These element names are then used to find in the database all the aqueous
# > species that could be formed out of them.
#
# ### Constructing the chemical system
#
# This step is where we create an object of class
# [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html) using the
# chemical system definition details stored in the object ``editor``.

system = ChemicalSystem(editor)

# ### Initial condition for the fluid composition
#
# Below, we define the **initial condition** for the fluid composition in the rock. We want an aqueous fluid that is
# 0.7 molal of NaCl and in equilibrium with calcite and quartz. To achieve this, we mix 1 kg of water, 0.7 mol of
# NaCl, and plenty of calcite and quartz (10 mol each) to ensure that the aqueous solution is saturated with respect
# to these minerals.

problem_ic = EquilibriumProblem(system)
problem_ic.setTemperature(T, 'celsius')
problem_ic.setPressure(P, 'bar')
problem_ic.add("H2O", 55.51, "mol") # H2O 55.51 M
problem_ic.add("NaCl", 0.27, "mol") # NaCl (aq) 0.27 M
problem_ic.add("KCl", 0.03, "mol") # KCl (aq)  0.03 M

# ### Boundary condition for the fluid composition
#
# Next, we define the **boundary condition** for the fluid composition on the left side of the rock, which should be
# the one that represents the fluid being continuously injected: 0.9 molal NaCl, 0.05 molal MgCl<sub>2</sub>,
# 0.01 molal CaCl<sub>2</sub> and almost CO<sub>2</sub>-saturated, with 0.75 molal of CO<sub>2</sub> dissolved. To
# provide that, we mix 1 kg of HO<sub>2</sub> with 0.9 mol of NaCl, 0.05 mol of MgCl<sub>2</sub>, 0.01 mol
# of CaCl<sub>2</sub>, and 0.75 mol of CO<sub>2</sub>.

problem_bc = EquilibriumProblem(system)
problem_bc.setTemperature(T, 'celsius')
problem_bc.setPressure(P, 'bar')
problem_bc.add("H2O", 55.51, "mol") # H2O 55.51 M
problem_bc.add("HCl", 0.1, "mol") # HCl (aq) 0.1 M
problem_bc.add("NaCl", 0.17, "mol") # NaCl (aq) 0.17 M
problem_bc.add("KCl", 0.03, "mol") # KCl (aq)  0.03 M

# ### Calculating the IC and BC fluid compositions
# In this step, we use the [equilibrate](https://reaktoro.org/cpp/namespaceReaktoro.html#af2d3b39d3e0b8f9cb5a4d9bbb06b697e)
# function to calculate the chemical equilibrium state of the system with the given initial and boundary equilibrium
# conditions stored in the object `problem_ic` and `problem_bc`, respectively. The result is stored in the
# corresponding instances of the class [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html),
# i.e., `state_ic` and `state_bc`.

state_ic = equilibrate(problem_ic)
state_bc = equilibrate(problem_bc)

# ### Scaling the phases in the initial condition
#
# The initial chemical state `state_ic` computed before has, at this point, phases with volumes that do not
# correspond to our desired porosity of 10 % and rock mineral composition of 98 %<sub>vol</sub> of quartz and
# 2 %<sub>vol</sub> of calcite.
#
# To obtain this, we scale the volumes of the aqueous and mineral phases as shown
# below:

# Scale the volumes of the phases in the initial condition
state_ic.scalePhaseVolume("Aqueous", 0.1, "m3") # 10% if the 1.0m3
state_ic.scalePhaseVolume("Quartz", 0.3 * 0.9, "m3") # 30% of 90% of remaining volume
state_ic.scalePhaseVolume("Muscovite", 0.05 * 0.9, "m3") # 5% of 90% of remaining volume
state_ic.scalePhaseVolume("Albite", 0.33 * 0.9, "m3") # 33% of 90% of remaining volume
state_ic.scalePhaseVolume("K-Feldspar", 0.32 * 0.9, "m3") # 32% of 90% of remaining volume


# > **Note**: After this scaling step, the sum of the phase volumes in ``state_ic`` is 1 m<sup>3</sup>. This also
# > ensures that the amounts of the species in the chemical system are normalized by m<sup>3</sup>, and thus they can
# > be regarded as concentrations in a unit of mol/m<sup>3</sup> (*bulk volume, not fluid volume!*).
#
# ### Scaling the boundary condition state
#
# Next, we scale the boundary condition state to 1 m<sup>3</sup>, so that we have the amounts of fluid species in
# `state_bc` also normalized by m<sup>3</sup>.
#
# > **Note**: The chemical state represented by `state_bc` has no other stable phase than the aqueous phase (i.e.,
# > all mineral phases have zero or negligible amounts such as 10<sup>-21</sup> mol).

state_bc.scaleVolume(1.0, 'm3')

# ### Creating the mesh
#
# We define the mesh with the class [Mesh](https://reaktoro.org/cpp/classReaktoro_1_1Mesh.html) to use in
# the initialization of class [ReactiveTransportSolver](
# https://reaktoro.org/cpp/classReaktoro_1_1ReactiveTransportSolver.html) later. Here, we specify the number of cells
# in the mesh and the x-coordinates of the left and right boundaries (in m).

mesh = Mesh(ncells, xl, xr)

# ### Creating a chemical field object
#
# We have been using class [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html) to represent
# an individual chemical state. We will now use class [ChemicalField](
# https://reaktoro.org/cpp/classReaktoro_1_1ChemicalField.html) to represent a collection of chemical states: one for
# each mesh cell.
#
# > **Note**: Different initial conditions across the mesh cells are possible by assigning different chemical states to
# > each mesh cell. Here, the same chemical state in `state_ic` is used for all cells.

field = ChemicalField(mesh.numCells(), state_ic)

# ### Initializing the reactive transport solver
#
# At last, we define the object responsible for solving the reactive transport problem, which is handled by the
# class [ReactiveTransportSolver](https://reaktoro.org/cpp/classReaktoro_1_1ReactiveTransportSolver.html).
# Here, we set the mesh and problem parameters such as velocity, diffusion coefficient, the chemical state
# representing the boundary condition, and the time step. We also initialize the reactive solver object with the
# chemical field object specified on the previous step, at this point containing the initial condition for theuse_smart_eqilibirum_solver
# chemical state of each mesh cell.

equilibrium_options = EquilibriumOptions()

reactive_transport_options = ReactiveTransportOptions()
reactive_transport_options.equilibrium = equilibrium_options

rtsolver = ReactiveTransportSolver(system)
rtsolver.setOptions(reactive_transport_options)
rtsolver.setMesh(mesh)
rtsolver.setVelocity(v)
rtsolver.setDiffusionCoeff(D)
rtsolver.setBoundaryState(state_bc)
rtsolver.setTimeStep(dt)
rtsolver.initialize()

profiler = ReactiveTransportProfiler()


# ### Defining the output quantities
#
# Before starting the reactive transport calculations, we define the quantities that will be output for every mesh
# cell, at every time step. For this, we use an object of the class
# [ChemicalOutput](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalOutput.html).
# The name of the output file is to `reactive-transport.txt`. We specify the parameters that we are interested in
# outputting. In this case, it is pH, molality of `H+`, `Ca++`, `Mg++`, `HCO3-`, `CO2(aq)`, as well as a phase volume
# of calcite and dolomite.

# +

# Create output class
output = rtsolver.output()
output.add("pH")
output.add("speciesMolality(Cl-)")
output.add("speciesMolality(HCl(aq))")
output.add("speciesMolality(Na+)")
output.add("speciesMolality(NaCl(aq))")
output.add("speciesMolality(OH-)")
output.add("speciesMolality(NaOH(aq))")
output.add("speciesMolality(NaHSiO3(aq))")
output.add("speciesMolality(K+)")
output.add("speciesMolality(KOH(aq))")
output.add("speciesMolality(KCl(aq))")
output.add("speciesMolality(Al+++)")
output.add("speciesMolality(AlOH++)")
output.add("speciesMolality(Quartz)")
output.add("speciesMolality(Diaspore)")
output.add("speciesMolality(Gibbsite)")
output.add("speciesMolality(Andalusite)")
output.add("speciesMolality(Kyanite)")
output.add("speciesMolality(Sillimanite)")
output.add("speciesMolality(Muscovite)")
output.add("speciesMolality(Paragonite)")
output.add("speciesMolality(Pyrophyllite)")
output.add("speciesMolality(Kaolinite)")
output.add("speciesMolality(Albite)")
output.add("speciesMolality(K-Feldspar)")
output.filename(folder_results + '/state.txt')  # Set the name of the output files
# -

# ### Running the reactive transport simulation
#
# As shown below, we perform a sequence of reactive transport calculations, one for each time step, during which the
# chemical state of each mesh cell is updated. The iterations continue until the maximum number of steps is
# achieved.
#
# Using **tqdm** we track the progress of simulations using the progress bar. For that, we wrap the while-loop with
# the function 'tqdm()'. We then use the method `step` of class
# [ReactiveTransportSolver](https://reaktoro.org/cpp/classReaktoro_1_1ReactiveTransportSolver.html) to perform a
# single  reactive transport time-stepping. This method also produces a new output file containing the requested
# output properties for every mesh cell. In each such file, rows correspond to cells, whereas the columns correspond
# to the requested output properties, i.e., pH, molality of `H+`, `Ca++`, `Mg++`, `HCO3-`, `CO2(aq)`, as
# well as the phase volume of calcite and dolomite.

# Step 15: Perform given number of reactive tranport steps
t = 0.0  # current time variable
step = 0  # current number of steps

#from tqdm.notebook import tqdm
#with tqdm(total=nsteps, desc="Reactive transport simulations") as pbar:
bar = IncrementalBar('Run reactive transport:', max = nsteps)
while step < nsteps:  # step until the number of steps are achieved

    # Perform one reactive transport time step
    rtsolver.step(field)

    # Update the profiler after every call to step method
    profiler.update(rtsolver.result())

    # Increment time step and number of time steps
    t += dt
    step += 1

    # Update a progress bar
    #pbar.update(1)
    bar.next()

bar.finish()

files = [file for file in natsorted( os.listdir(folder_results) ) ]

indx_ph         = 0
indx_Clanion    = 1
indx_HClaq      = 2
indx_Nacation   = 3
indx_NaClaq     = 4
indx_OHanion    = 5
indx_NaOHaq     = 6
indx_NaHSiO3aq  = 7
indx_Kcation    = 8
indx_KOHaq      = 9
indx_KClaq      = 10
indx_Alcation   = 11
indx_AlOHcation = 12

indx_Quartz       = 13
indx_Diaspore     = 14
indx_Gibbsite     = 15
indx_Andalusite   = 16
indx_Kyanite      = 17
indx_Sillimanite  = 18
indx_Muscovite    = 19
indx_Paragonite   = 20
indx_Pyrophyllite = 21
indx_Kaolinite    = 22
indx_Albite       = 23
indx_KFeldspar    = 24

#plot_at_selected_steps = [1, 10, 60, 100]
plot_at_selected_steps = [1, 10, 60, 120, 240, 480, 860, 960, 1000]

def titlestr(t):
    d = int(t / day)                 # The number of days
    h = int(int(t % day) / hour)     # The number of remaining hours
    m = int(int(t % hour) / minute)  # The number of remaining minutes
    return '{:>3d}d {:>2}h {:>2}m'.format(int(d), str(int(h)).zfill(2), str(int(m)).zfill(2))

def line(color):
    return {'linestyle': '-', 'color': color, 'zorder': 1, 'linewidth': 2}

def plot_figures_quartz():

    for i in plot_at_selected_steps:
        t = i * dt
        filearray = np.loadtxt(folder_results + '/' + files[i-1], skiprows=1)
        data = filearray.T
        data_quartz = data[indx_Quartz]

        plt.axes(xlim=(xl - 0.01, xr + 0.01), ylim=(130.0, 190.0))
        plt.ylabel('Concentration [mol/m3]')
        plt.xlabel('Distance [m]')
        plt.title(titlestr(t))
        plt.plot(xcells, data_quartz, label='Quartz', **line('darkviolet'))
        plt.legend(loc='center right')
        plt.savefig(folder_result_plots + '/quarzt-{}.png'.format(i))
        plt.tight_layout()
        plt.close()

def plot_figures_kfeldspar():

    for i in plot_at_selected_steps:
        t = i * dt
        filearray = np.loadtxt(folder_results + '/' + files[i-1], skiprows=1)
        data = filearray.T
        data_kfeldspar = data[indx_KFeldspar]

        plt.axes(xlim=(xl - 0.01, xr + 0.01), ylim=(22.0, 38.0))
        plt.ylabel('Concentration [mol/m3]')
        plt.xlabel('Distance [m]')
        plt.title(titlestr(t))
        plt.plot(xcells, data_kfeldspar, label='K-Feldspar', **line('C3'))
        plt.legend(loc='center right')
        plt.savefig(folder_result_plots + '/kfeldspar-{}.png'.format(i))
        plt.tight_layout()
        plt.close()

def plot_figures_muscovite():

    for i in plot_at_selected_steps:
        t = i * dt
        filearray = np.loadtxt(folder_results + '/' + files[i-1], skiprows=1)
        data = filearray.T
        data_muscovite = data[indx_Muscovite]

        plt.axes(xlim=(xl - 0.01, xr + 0.01), ylim=(3.0, 14.0))
        plt.ylabel('Concentration [mol/m3]')
        plt.xlabel('Distance [m]')
        plt.title(titlestr(t))
        plt.plot(xcells, data_muscovite, label='Muscovite', **line('C1'))
        plt.legend(loc='center right')
        plt.savefig(folder_result_plots + '/muscovite-{}.png'.format(i))
        plt.tight_layout()
        plt.close()

def plot_figures_albite():

    for i in plot_at_selected_steps:
        t = i * dt
        filearray = np.loadtxt(folder_results + '/' + files[i-1], skiprows=1)
        data = filearray.T
        data_albite = data[indx_Albite]

        plt.axes(xlim=(xl - 0.01, xr + 0.01), ylim=(10.0, 42.0))
        plt.ylabel('Concentration [mol/m3]')
        plt.xlabel('Distance [m]')
        plt.title(titlestr(t))
        plt.plot(xcells, data_albite, label='Albite', **line('C2'))
        plt.legend(loc='center right')
        plt.savefig(folder_result_plots + '/albite-{}.png'.format(i))
        plt.tight_layout()
        plt.close()

plot_figures_albite()
plot_figures_muscovite()
plot_figures_kfeldspar()
plot_figures_quartz()