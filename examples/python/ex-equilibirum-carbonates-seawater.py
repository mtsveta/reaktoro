# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright (C) 2014-2020 Allan Leal
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


from reaktoro import *

import numpy as np
import matplotlib.pyplot as plt

def carbonates_in_seawater(system, T, P):

    n0Calcite = 10.0
    n0Dolomite = 10.0
    water_kg = 1.0

    state_sw = ChemicalState(system)
    state_sw.setTemperature(T, "celsius")
    state_sw.setPressure(P, "atm")
    state_sw.setSpeciesMass("H2O", 1.0, "kg")
    state_sw.setSpeciesMass("Ca+2", 412.3 * water_kg, "mg")
    state_sw.setSpeciesMass("Mg+2", 1290 * water_kg, "mg")
    state_sw.setSpeciesMass("Na+", 10768.0 * water_kg, "mg")
    state_sw.setSpeciesMass("K+", 399.1 * water_kg, "mg")
    state_sw.setSpeciesMass("Cl-", 19353.0 * water_kg, "mg")
    state_sw.setSpeciesMass("HCO3-", 141.682 * water_kg, "mg")
    state_sw.setSpeciesMass("SO4-2", 2712.0 * water_kg, "mg")

    # Calculate chemical state corresponding to the carbonate water
    solver = EquilibriumSolver(system)

    res = solver.solve(state_sw)
    if not res.optima.succeeded:
        raise RuntimeError("Equilibrium calculation did not succeed!")

    state_sw.setSpeciesAmount("Dolomite", n0Dolomite, "mol")
    state_sw.setSpeciesAmount("Calcite", n0Calcite, "mol")

    # Equilibrate the seawater with carbonates
    solver.solve(state_sw)

    nDolomite = state_sw.speciesAmount("Dolomite")
    nCalcite = state_sw.speciesAmount("Calcite")

    nCa2 = state_sw.speciesAmount("Ca+2")
    nMg2 = state_sw.speciesAmount("Mg+2")
    nH = state_sw.speciesAmount("H+")
    nHCO3 = state_sw.speciesAmount("HCO3-")

    return (nCalcite[0], nDolomite[0], nCa2[0], nMg2[0], nH[0], nHCO3[0])

#database = "phreeqc"
database = "pitzer"

db = PhreeqcDatabase(database + ".dat")

aqueousphase = AqueousPhase(speciate("H O C Ca Cl Na K Mg S Si"))

if database == "phreeqc":
    aqueousphase.setActivityModel(chain(
        ActivityModelHKF(),
        ActivityModelDrummond("CO2"),
    ))
else:
    aqueousphase.setActivityModel(chain(
        ActivityModelPitzerHMW()
    ))

calcitephase = MineralPhase("Calcite")
dolomitephase = MineralPhase("Dolomite")

phases = Phases(db)
phases.add(aqueousphase)
phases.add(calcitephase)
phases.add(dolomitephase)

system = ChemicalSystem(phases)

T = np.arange(25.0, 91.0, 5.0)
P = 1.0
#P = 100.0

mCalcite_Dolomite = [carbonates_in_seawater(system, x, P) for x in T]  # [0] is needed to get the value of autodiff.real
mCalcite = [molals[0] for molals in mCalcite_Dolomite]
mDolomite = [molals[1] for molals in mCalcite_Dolomite]

mCa2 = [molals[2] for molals in mCalcite_Dolomite]
mMg2 = [molals[3] for molals in mCalcite_Dolomite]
mH = [molals[4] for molals in mCalcite_Dolomite]
mHCO3 = [molals[5] for molals in mCalcite_Dolomite]

# print(mCa2)
# print(mMg2)
# print(mH)
# print(mHCO3)

from numpy import savetxt
savetxt('reaktoro-seawater-' + database + '-Calcite-p-' + str(P) + '.txt', mCalcite)
savetxt('reaktoro-seawater-' + database + '-Dolomite-p-' + str(P) + '.txt', mDolomite)
savetxt('reaktoro-seawater-' + database + '-Ca2-p-' + str(P) + '.txt', mCa2)
savetxt('reaktoro-seawater-' + database + '-Mg2-p-' + str(P) + '.txt', mMg2)
savetxt('reaktoro-seawater-' + database + '-H-p-' + str(P) + '.txt', mH)
savetxt('reaktoro-seawater-' + database + '-HCO3-p-' + str(P) + '.txt', mHCO3)

fig, ax = plt.subplots()

ax.plot(T, mCalcite, label=f"Calcite molal")
ax.plot(T, mDolomite, label=f"Dolomite molal")

ax.legend(loc="upper right")

ax.set(xlabel='Temperature [degC]', ylabel='Solubility [mol/kgw]',
       title='Carbonates equilibration with seawater')
ax.grid()

fig.savefig('carbonates-in-seawater-' + database + '-p-' + str(P) + '.png')


fig, ax = plt.subplots()

plt.yscale("log")
ax.plot(T, mCa2, label=f"Ca2+ molal")
ax.plot(T, mMg2, label=f"Mg2+ molal")
ax.plot(T, mHCO3, label=f"HCO3- molal")
ax.plot(T, mH, label=f"H+ molal")

ax.legend(loc="center right")

ax.set(xlabel='Temperature [degC]', ylabel='Solubility [mol/kgw]',
       title='Aqueous species in carbonates equilibration with seawater')
ax.grid()

fig.savefig('aqueous-species-carbonates-in-seawater-' + database + '-p-' + str(P) + '.png')