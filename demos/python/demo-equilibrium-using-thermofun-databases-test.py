from reaktoro import *
import thermofun.PyThermoFun as thermofun

T = 500.0   # in celsius
P = 3000.0  # in bar

#database = thermofun.Database("databases/thermofun/aq17-thermofun.json") # OK
#database = thermofun.Database("databases/thermofun/heracles-thermofun.json") # Error
'''
*******************************************************************************
*** Error: Could not parse the species `Al2F6(l) in the database.
*** Reason: The type of the species is unknown.
*** Location:  This error was encountered in Reaktoro/Thermodynamics/Core/Database.cpp:489.
*******************************************************************************
'''
database = thermofun.Database("databases/thermofun/cemdata18-thermofun.json") # Error
'''
************************************************************************************************************************************
*** Error: Unknown JSON type  
*** Reason: The JSON object needs to be an element, a substance or a reaction file databases/thermofun/cemdata18-thermofun.json.
*** Location: :311
************************************************************************************************************************************
'''
#database = thermofun.Database("databases/thermofun/mines16-thermofun.json") # OK
#database = thermofun.Database("databases/thermofun/psinagra07-thermofun.json") # Error
'''
*************************************************************************************************************************************
*** Error: Unknown JSON type  
*** Reason: The JSON object needs to be an element, a substance or a reaction file databases/thermofun/psinagra07-thermofun.json.
*** Location: :311
*************************************************************************************************************************************
'''

editor = ChemicalEditor(database)
editor.setTemperatures([T], "celsius")
editor.setPressures([P], "bar")

editor.addAqueousPhase(["OH-", "H+", "H2O@", "Na+", "Cl-", "Ca+2", "Mg+2", "HCO3-", "CO2@", "CO3-2"])
editor.addMineralPhase("Calcite")
#editor.addMineralPhase("Cal")
editor.addMineralPhase("Dolomite")
editor.addMineralPhase("Quartz")

system = ChemicalSystem(editor)
print(system)

problem = EquilibriumProblem(system)
problem.add("H2O",   1.0, "kg");
problem.add("NaCl",  0.7, "mol");
problem.add("CaCO3", 10,  "mol");
problem.add("SiO2",  10,  "mol");
problem.setTemperature(T, "celsius")
problem.setPressure(P, "bar")

state = equilibrate(problem)
print("Species name : n (in mol)")
for species in system.species():
    name = species.name()
    amount = state.speciesAmount(name)
    # Output according to the threshold
    print(f"{name:>13} = {amount}")
state.output("result.txt")