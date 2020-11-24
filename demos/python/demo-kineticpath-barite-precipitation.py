# Import reaktoro and math packages functionality:

from reaktoro import *
from math import *

# Define the time interval of the kinetic simulation and corresponding to it file-name (for the future results to be
# stored):

# +
# Construct the chemical system with its phases and species
db = Database('supcrt07.xml')

# Fetch Debye-Huckel activity model parameters
dhModel = DebyeHuckelParams()
dhModel.setPHREEQC()

# Define the content of the chemical system
editor = ChemicalEditor(db)
editor.addAqueousPhaseWithElements("H Cl S O Ba Ca Sr Na K Mg C Si"). \
    setChemicalModelDebyeHuckel(dhModel)
editor.addMineralPhase('Barite')    # BaSiO4

# Create chemical system
system = ChemicalSystem(editor)

# Ionic strength function
I = ChemicalProperty.ionicStrength(system)
# -

# We set the reaction equations, and particular parameters for calcite, dolomite, halite, k-feldspar, quartz, and
# kaolinite using chemical editor:

# +
eq_str_barite = "Barite = SO4-- + Ba++"
min_reaction_barite = editor.addMineralReaction("Barite") \
    .setEquation(eq_str_barite) \
    .addMechanism("logk = -8.6615 mol/(m2*s); Ea = 22 kJ/mol") \
    .setSpecificSurfaceArea(0.006, "m2/g")
# Create reaction for barite
reaction_barite = createReaction(min_reaction_barite, system)
reaction_barite.setName("Barite kinetics")
# Barite kinetic rate function
def rate_func_barite_shell(properties):


    # The number of chemical species in the system
    num_species = system.numSpecies()

    # Final mineral reaction rate using specified surface area
    res = ChemicalScalar(num_species)
    # Auxiliary variable for calculating mineral reaction rate
    f = ChemicalScalar(num_species, 1.0)

    # The universal gas constant (in units of kJ/(mol*K))
    R = 8.3144621e-3

    # The temperature and pressure of the system
    T = properties.temperature().val

    # Calculate the saturation index of the mineral
    lnK = reaction_barite.lnEquilibriumConstant(properties)
    lnQ = reaction_barite.lnReactionQuotient(properties)
    lnOmega = lnQ.val - lnK.val

    # Calculate the saturation index
    Omega = exp(lnOmega)

    # The composition of the chemical system
    n = properties.composition()

    # The initial composition of the chemical system
    n0 = reaction_barite.initialAmounts()

    # The index of the mineral
    imineral = system.indexSpeciesWithError(min_reaction_barite.mineral())

    # The number of moles of the mineral
    nm = n[imineral].val
    nm0 = n0[imineral]

    # Prevent negative mole numbers here for the solution of the ODEs
    nm = max(nm, 0.0)
    nm0 = max(nm0, 0.0)

    # Get H+ activity and ionic strength
    lna = properties.lnActivities().val
    i_h = system.indexSpeciesWithError("H+")

    activity_h = exp(lna[i_h])
    ionic_strength = I(properties).val

    # The molar mass of the mineral (in units of kg/mol)
    molar_mass = system.species(imineral).molarMass()

    # If (m <= 0) and (SRmin < 1) Then GoTo 240
    # If (SRmin = 1) Then GoTo 240
    if ((nm <= 0 and Omega < 1) or (Omega == 1)): # the is no way to precipitate further
        res = 0
    ssa = 0.006

    # If (SRmin > 1) Then GoTo 130
    if Omega > 1: # precipitation kinetics
        # If (m <= 1e-5) then GoTo 170
        if nm > 1e-5:
            kappa_pre = -2.18e-9 * exp(- 22 / R * (1.0/T - 1.0/298.15)) \
                        * pow(10, 0.6 * pow(ionic_strength, 0.5))
            res += f * ssa * nm * molar_mass * kappa_pre * (Omega - 1)
        else:
            # Set nucleation rate
            res = -1e-12 * f

    else: # dissolution kinetics
        kappa = 2.75e-8 * exp(- 25 / R * (1.0/T - 1.0/298.15)) \
                * pow(activity_h, 0.03) \
                * pow(10, 0.6 * pow(ionic_strength, 0.5))
        res += f * ssa * nm * molar_mass * pow(nm / nm0, 2 / 3) * kappa * (1 - pow(Omega, 0.2))

    return res

reaction_barite.setRate(rate_func_barite_shell)
# -

reactions = ReactionSystem(system, [reaction_barite])

# Specifying the partition including the kinetic species:
partition = Partition(system)
partition.setKineticPhases(["Barite"])

T = 60.0    # temperature (in units of celcius)
P = 200.0   # pressure (in units of atm)
water_kg = 1.00

# +
problem_ic = EquilibriumInverseProblem(system)
problem_ic.setTemperature(T, "celsius")
problem_ic.setPressure(P, "atm")
problem_ic.add('H2O', water_kg, 'kg')
problem_ic.add("SO4", 10 * water_kg, "ug")
problem_ic.add("Ca", 995 * water_kg, "mg")
problem_ic.add("Ba", 995 * water_kg, "mg")
problem_ic.add("Sr", 105 * water_kg, "mg")
problem_ic.add("Na", 27250 * water_kg, "mg")
problem_ic.add("K", 1730 * water_kg, "mg")
problem_ic.add("Mg", 110 * water_kg, "mg")
problem_ic.add("Cl", 45150 * water_kg, "mg")
problem_ic.add("HCO3", 1980 * water_kg, "mg")
problem_ic.pH(7.0, "HCl", "NaOH")

# Calculate the equilibrium states for the initial conditions
state_ic = equilibrate(problem_ic)
state_ic.setSpeciesAmount("Barite", 0.1, "mcmol")
state_ic.scaleVolume(1.0, "m3")

reaction_barite.setInitialAmounts(state_ic.speciesAmounts())
# -

# +
problem_bc = EquilibriumInverseProblem(system)
problem_bc.setTemperature(T, "celsius")
problem_bc.setPressure(P, "atm")
problem_bc.add('H2O', water_kg, 'kg')
problem_bc.add("SO4--", 2710 * water_kg, "mg")
problem_bc.add("Ca++", 411 * water_kg, "mg")
problem_bc.add("Ba++", 0.01 * water_kg, "mg")
problem_bc.add("Sr++", 8 * water_kg, "mg")
problem_bc.add("Na+", 10760 * water_kg, "mg")
problem_bc.add("K+", 399 * water_kg, "mg")
problem_bc.add("Mg++", 1290 * water_kg, "mg")
problem_bc.add("Cl-", 19350 * water_kg, "mg")
problem_bc.add("HCO3-", 142 * water_kg, "mg")
problem_bc.pH(8.1, "HCl", "NaOH")

# Calculate the equilibrium states for the initial conditions
state_bc = equilibrate(problem_bc)
state_bc.scaleVolume(1.0, "m3")
# -

# +
# Performing the kinetic path calculation:
path = KineticPath(reactions, partition)

# To analyse the result of kinetic simulations, we save the evolution of different properties of the chemical system
# into file `result_file_name`:

result_file_name = "kineticpath-barite-kinetics-n-steps.txt"
output = path.output()
output.filename(result_file_name)
output.add("time(units=s)")
output.add("speciesAmount(Barite units=mol)", "Barite")
output.add("speciesAmount(Ba++ units=mol)", "Ba++")
output.add("speciesAmount(SO4-- units=mol)", "SO4--")
output.filename(result_file_name)
# -

# Solving the chemical kinetics problem:

t0 = 0
dt = 5.0
tfinal = 50.0
n = int(tfinal / dt)
state_ic = state_ic + state_bc
#path.solve(state_ic, t0, t1, "minute")
path.solve(state_ic, t0, dt, n, "day")

# For plotting of the results of equilibrium path calculation, we load the results into the `data` array:

filearray = numpy.loadtxt(result_file_name, skiprows=1) # load data from the file skipping the one row
data = filearray.T  # transpose the matrix with data
[time_indx, barite_indx, ba_indx, so4_indx] = numpy.arange(0, 4)

import os
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm, rcParams
import matplotlib as mpl

fpath = os.path.join(rcParams["datapath"], "texgyreadventor-regular.otf")
prop = fm.FontProperties(fname=fpath)
prop.set_size(14)

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'TeX Gyre Adventor'
mpl.rcParams['font.style'] = 'normal'
mpl.rcParams['font.size'] = 14
mpl.set_loglevel("critical")

# +
time = data[time_indx, :]  # fetch time from the data matrix

plt.xlabel(r'Time [s]', fontproperties=prop)
plt.ylabel('Amount [mol]', fontproperties=prop)
plt.title("Minerals concentration w.r.t. time", fontproperties=prop)
plt.plot(time, data[barite_indx], label='Barite')
plt.legend(loc='center right', prop=prop)
plt.show()
# -