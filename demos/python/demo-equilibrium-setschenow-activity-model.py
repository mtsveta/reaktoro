from reaktoro import *

db = Database("supcrt98.xml")

editor = ChemicalEditor(db)
aqueousphase = editor.addAqueousPhaseWithElements("H O Na Cl")
aqueousphase.setActivityModelSetschenow("NaCl(aq)", 1.0)
editor.addMineralPhase("Halite")

system = ChemicalSystem(editor)

T_left = 30
T_right = 90
temperature_points = 7

import numpy as np
temperatures = np.linspace(T_left, T_right, temperature_points)  # the x-coordinates of the plots
molalities = [0.5, 1.0, 2.0]
species_amounts = np.zeros((len(molalities), len(temperatures)))
for molality, i in zip(molalities, list(range(len(molalities)))):
    for T, j in zip(temperatures, list(range(len(temperatures)))):

        problem = EquilibriumProblem(system)
        problem.setTemperature(T, "celsius")
        problem.setPressure(1.0, "bar")
        problem.add("H2O", 1.0, "kg")
        problem.add("NaCl", molality, "mol")

        state = equilibrate(problem)
        species_amounts[i, j] = state.speciesAmount("Na+")

    print(*species_amounts[i], sep=', ')