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


def carbonates_dissolution_in_co2_saturated_water(system, T, P):

    n0CO2g = 10.0
    n0Calcite = 10.0
    n0Dolomite = 10.0
    water_kg = 1.0

    # Define solution corresponding to the carbonate water
    state_cw = ChemicalState(system)
    state_cw.setTemperature(T, "celsius")
    state_cw.setPressure(P, "atm")
    state_cw.setSpeciesMass("H2O", water_kg, "kg")
    state_cw.setSpeciesAmount("CO2(g)", n0CO2g, "mol")
    state_cw.setSpeciesAmount("Calcite", n0Calcite, "mol")
    state_cw.setSpeciesAmount("Dolomite", n0Dolomite, "mol")

    options = EquilibriumOptions()
    options.optima.output.active = False
    options.optima.max_iterations = 100
    options.optima.kkt.method = optima.SaddlePointMethod.Nullspace
    options.optima.linesearch.trigger_when_current_error_is_greater_than_initial_error_by_factor = 1000000000000000000000.0  # zero to disable line search
    options.optima.linesearch.trigger_when_current_error_is_greater_than_previous_error_by_factor = 1000000000000000000000.0
    options.epsilon = 1e-40

    # Calculate chemical state corresponding to the carbonate water
    solver = EquilibriumSolver(system)
    solver.setOptions(options)

    res = solver.solve(state_cw)

    if not res.optima.succeeded:
        raise RuntimeError("Equilibrium calculation did not succeed!")

    nCO2g = state_cw.speciesAmount("CO2(g)")
    nCalcite = state_cw.speciesAmount("Calcite")
    nDolomite = state_cw.speciesAmount("Dolomite")

    # print("delta CO2g = ", n0CO2g - nCO2g)
    # print("delta Calcite = ", n0Calcite - nCalcite)
    # print("delta Dolomite = ", n0Dolomite - nDolomite)
    # print("delta Dolomite = ", ((n0CO2g - nCO2g)[0], (n0Calcite - nCalcite)[0], (n0Dolomite - nDolomite)[0]))
    # input()

    # n = state_cw.speciesAmounts()
    #
    # for i in range(len(n)):
    #     print(f"{system.species(i).name()} {n[i]}")
    #
    # input()

    return ((n0CO2g - nCO2g)[0], (n0Calcite - nCalcite)[0], (n0Dolomite - nDolomite)[0])

database = "phreeqc"
#database = "pitzer"

db = PhreeqcDatabase(database + ".dat")

aqueousphase = AqueousPhase(speciate("H O C Ca Mg"))

if database == "phreeqc":
    aqueousphase.setActivityModel(chain(
        ActivityModelHKF(),
        ActivityModelDrummond("CO2"),
    ))
else:
    aqueousphase.setActivityModel(chain(
        ActivityModelPitzerHMW()
    ))

gaseousphase = GaseousPhase("CO2(g)")
gaseousphase.setActivityModel(ActivityModelPengRobinson())

calcitephase = MineralPhase("Calcite")
dolomitephase = MineralPhase("Dolomite")

phases = Phases(db)
phases.add(aqueousphase)
phases.add(gaseousphase)
phases.add(calcitephase)
phases.add(dolomitephase)

system = ChemicalSystem(phases)

T = np.arange(25.0, 91.0, 5.0)
P = 1.0
#P = 100.0

deltaCO2_Calcite_Dolomite = [carbonates_dissolution_in_co2_saturated_water(system, x, P) for x in T]  # [0] is needed to get the value of autodiff.real
deltaCO2 = [molals[0] for molals in deltaCO2_Calcite_Dolomite]
deltaCalcite = [molals[1] for molals in deltaCO2_Calcite_Dolomite]
deltaDolomite = [molals[2] for molals in deltaCO2_Calcite_Dolomite]

from numpy import savetxt
savetxt('reaktoro-' + database + '-delta-CO2g-p-' + str(P) + '.txt', deltaCO2)
savetxt('reaktoro-' + database + '-delta-Calcite-p-' + str(P) + '.txt', deltaCalcite)
savetxt('reaktoro-' + database + '-delta-Dolomite-p-' + str(P) + '.txt', deltaDolomite)

fig, ax = plt.subplots()

ax.plot(T, deltaCO2, label=f"CO2(g) molal")
ax.plot(T, deltaCalcite, label=f"Calcite molal")
ax.plot(T, deltaDolomite, label=f"Dolomite molal")

ax.legend(loc="upper right")

ax.set(xlabel='Temperature [degC]', ylabel='Solubility [mol/kgw]',
       title='Carbonates dissolution in CO2 saturated water')
ax.grid()

fig.savefig('carbonates-dissolution-' + database + '-p-' + str(P) + '.png')