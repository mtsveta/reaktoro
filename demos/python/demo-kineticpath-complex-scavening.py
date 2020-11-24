
# In this tutorial, we expand the complexity of the souring example with kinetically controlled mineral phases, such
# that, chlorite, kaolinite, quartz, siderite, and calcite. We keep precipitating pyrrhotite as the equilibrium phase.
# Calcite can be represented as an equilibrium phase, but for the testing purposes, we represent it using kinetics.

# First, we import **reaktoro** and **math** packages functionality:

from reaktoro import *
from math import *

# Second, we construct the chemical system with its phases and species.

# +
# Define the database to be used
db = Database('supcrt07.xml')

# Fetch Debye-Huckel activity model parameters
dhModel = DebyeHuckelParams()
dhModel.setPHREEQC()

# Define the content of the chemical system
editor = ChemicalEditor(db)
editor.addAqueousPhaseWithElements("Al C Ca Cl Fe H K Mg Na O S Si"). \
    setChemicalModelDebyeHuckel(dhModel)
editor.addMineralPhase("Pyrrhotite")
editor.addMineralPhase("Calcite")
editor.addMineralPhase("Daphnite,14A") # Chlorite
editor.addMineralPhase("Siderite")
editor.addMineralPhase("Kaolinite")
editor.addMineralPhase("Quartz")

# Create chemical system
system = ChemicalSystem(editor)
print(system)
# -

# Next, we define the parameters of the kinetically controlled minerals.
# First, the kinetic parameters of the calcite are defined:

# +
# Define reaction equation
eq_str_calcite = "Calcite + H+  =  Ca++ + HCO3-"
# Define mineral reaction
min_reaction_calcite = editor.addMineralReaction("Calcite") \
    .setEquation(eq_str_calcite) \
    .addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol") \
    .addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0") \
    .setSpecificSurfaceArea(0.7, "m2/g")
# Create reaction
reaction_calcite = createReaction(min_reaction_calcite, system)
reaction_calcite.setName("Calcite kinetics")
# Kinetic rate function
def rate_func_calcite_shell(properties):

    # The number of chemical species in the system
    num_species = system.numSpecies()

    # Reaction rate variables
    res = ChemicalScalar(num_species)
    res_growth = ChemicalScalar(num_species)
    res_nuc = ChemicalScalar(num_species)

    # Auxiliary variable for calculating mineral reaction rate
    f = ChemicalScalar(num_species, 1.0)

    # The universal gas constant (in units of kJ/(mol*K))
    R = 8.3144621e-3

    # The temperature and pressure of the system
    T = properties.temperature().val

    # Calculate the saturation index of the mineral
    lnK = reaction_calcite.lnEquilibriumConstant(properties)
    lnQ = reaction_calcite.lnReactionQuotient(properties)
    lnOmega = lnQ.val - lnK.val

    # Calculate the saturation index
    Omega = exp(lnOmega)

    # The composition of the chemical system
    n = properties.composition()

    # The initial composition of the chemical system
    n0 = reaction_calcite.initialAmounts()

    print(min_reaction_calcite)

    # The index of the mineral
    imineral = system.indexSpeciesWithError("Calcite")

    print(imineral)
    input()

    # The number of moles of the mineral
    nm = n[imineral].val
    nm0 = n0[imineral]

    # Prevent negative mole numbers here for the solution of the ODEs
    nm = max(nm, 0.0)
    nm0 = max(nm0, 0.0)

    # Get H+ activity
    lna = properties.lnActivities().val
    i_h = system.indexSpeciesWithError("H+")
    activity_h = exp(lna[i_h])

    # The molar mass of the mineral (in units of kg/mol)
    molar_mass = system.species(imineral).molarMass()

    # If (m <= 0) and (SRmin < 1) Then GoTo 350
    # If (SRmin = 1) Then GoTo 350
    if ((nm <= 0 and Omega < 1) or (Omega == 1)): # the is no way to precipitate further
        res = 0

    # Average BET value in m2/g => 1.7 * 1e3 m2/kg
    ssa = 0.7 * 1e3
    # Molar volume
    molar_volume = 3.693e-05

    # If (SRmin > 1) Then GoTo 130
    if Omega > 1: # precipitation kinetics


        # Interfacial energy in J m-2
        surf_energy = 0.047

        # Molecular volume in m3
        molecular_volume = 6.13e-29

        u = 16 * 3.14159 * pow(surf_energy, 3.0) * pow(molecular_volume, 2.0) / (3 * pow(1.38e-23 * T, 3.0))

        # Nucleation rate in nuclei kg-1 sec-1 (after Fritz et al 2009)
        j0 = 1e20

        # Critical saturation threshold
        Omega_crit = exp(sqrt(u / log10(j0)))


        if Omega < Omega_crit:

            # Nucleation rate in nuclie sec-1 multiplied by time to account for nucleation over time
            nuclie = j0 * exp(- u / pow(log10(Omega), 2.0))

            # Number of growth units in critical nuclie radius
            nj = 2 * u / pow(log10(Omega), 3.0)

            # Critical nuclie volume
            vol = nj * molecular_volume

            # Moles of nuclie formed
            res_nuc += - f * nuclie * vol / molar_volume

        # Define rate constant of kinetic growth
        kappa_neu = -1.6E-11 * exp(- 108 / R * (1.0 / T - 1.0 / 298.15))
        kappa_acid = 5.9E-06 * exp(- 56.0 / R * (1.0/T - 1.0/298.15)) * pow(activity_h, 0.6)
        kappa_growth = - (kappa_neu + kappa_acid)

        # Calculate kinetic rate of the mineral growth
        res_growth += f * ssa * nm * molar_mass * kappa_growth * pow(pow(Omega, 0.5) - 1, 1.12)

        # Sum the rate of growth and nucleation
        res += res_growth + res_nuc

    else: # dissolution kinetics

        # Define rate constants of mineral dissolution
        kappa_neu = 1.6E-6 * exp(- 24.0 / R * (1.0/T - 1.0/298.15))
        kappa_acid = 5E-1 * exp(- 14.0 / R * (1.0/T - 1.0/298.15)) * activity_h
        kappa = kappa_neu + kappa_acid

        # Calculate kinetic rate of mineral dissolution
        res += f * ssa * nm * molar_mass * pow(nm / nm0, 2 / 3) * kappa * (1 - Omega)

    return res
# Set kinetic rate function
reaction_calcite.setRate(rate_func_calcite_shell)
# -

# Define kinetic parameters of the daphnite:

# +
# Define reaction equation
eq_str_daphnite = "Daphnite,14A + 16*H+ = 2*Al+++ + 5*Fe++ + 3*SiO2(aq) + 12*H2O(l)"
# Define mineral reaction
min_reaction_daphnite = editor.addMineralReaction("Daphnite,14A") \
    .setEquation(eq_str_daphnite) \
    .setSpecificSurfaceArea(0.0027, "m2/g")
# Create reaction
reaction_daphnite = createReaction(min_reaction_daphnite, system)
reaction_daphnite.setName("Daphnite kinetics")
# Kinetic rate function
def rate_func_daphnite_shell(properties):

    # The number of chemical species in the system
    num_species = system.numSpecies()

    # Reaction rate variables
    res = ChemicalScalar(num_species)

    # Auxiliary variable for calculating mineral reaction rate
    f = ChemicalScalar(num_species, 1.0)

    # The universal gas constant (in units of kJ/(mol*K))
    R = 8.3144621e-3

    # The temperature and pressure of the system
    T = properties.temperature().val

    # Calculate the saturation index of the mineral
    lnK = reaction_daphnite.lnEquilibriumConstant(properties)
    lnQ = reaction_daphnite.lnReactionQuotient(properties)
    lnOmega = lnQ.val - lnK.val

    # Calculate the saturation index
    Omega = exp(lnOmega)

    # The composition of the chemical system
    n = properties.composition()

    # The initial composition of the chemical system
    n0 = reaction_daphnite.initialAmounts()

    # The index of the mineral
    imineral = system.indexSpeciesWithError(min_reaction_daphnite.mineral())

    # The number of moles of the mineral
    nm = n[imineral].val
    nm0 = n0[imineral]

    # Prevent negative mole numbers here for the solution of the ODEs
    nm = max(nm, 0.0)
    nm0 = max(nm0, 0.0)

    # Get H+ activity and OH- activity
    lna = properties.lnActivities().val
    i_h = system.indexSpeciesWithError("H+")
    activity_h = exp(lna[i_h])
    i_oh = system.indexSpeciesWithError("OH-")
    activity_oh = exp(lna[i_oh])

    # The molar mass of the mineral (in units of kg/mol)
    molar_mass = system.species(imineral).molarMass()

    # If (m <= 0) and (SRmin < 1) Then GoTo 350
    # If (SRmin = 1) Then GoTo 350
    if ((nm <= 0 and Omega < 1) or (Omega == 1)): # the is no way to precipitate further
        res = 0

    # Specific surface area: 0.2% BET (03bra/bos); suggested value in 0.0027 m2/g => 0.0027 m2/kg
    ssa = 0.0027 * 1e3

    # If (SRmin > 1) Then GoTo 130
    if Omega < 1: # dissolution kinetics

        # Define rate constants of mineral dissolution
        kappa_neu = 6.4E-17 * exp(- 16.0 / R * (1.0/T - 1.0/298.15))
        kappa_acid = 8.2E-09 * exp(- 17.0 / R * (1.0/T - 1.0/298.15)) * pow(activity_h, 0.28)
        kappa_oh = 6.9E-09 * exp(- 16.0 / R * (1.0/T - 1.0/298.15)) * pow(activity_oh, 0.34)
        kappa = kappa_neu + kappa_acid + kappa_oh

        # Calculate kinetic rate of mineral dissolution
        res += f * ssa * nm * molar_mass * pow(nm / nm0, 2 / 3) * kappa * (1 - Omega)

    return res
# Set kinetic rate function
reaction_daphnite.setRate(rate_func_daphnite_shell)
# -

# Define kinetic parameters of the siderite:

# +
# Define reaction equation
eq_str_siderite = "Siderite + H+ = HCO3- + Fe++"
# Define mineral reaction
min_reaction_siderite = editor.addMineralReaction("Siderite") \
    .setEquation(eq_str_siderite) \
    .setSpecificSurfaceArea(0.13, "m2/g")
# Create reaction
reaction_siderite = createReaction(min_reaction_siderite, system)
reaction_siderite.setName("Siderite kinetics")
# Kinetic rate function
def rate_func_siderite_shell(properties):

    # The number of chemical species in the system
    num_species = system.numSpecies()

    # Reaction rate variables
    res = ChemicalScalar(num_species)
    res_growth = ChemicalScalar(num_species)
    res_nuc = ChemicalScalar(num_species)

    # Auxiliary variable for calculating mineral reaction rate
    f = ChemicalScalar(num_species, 1.0)

    # The universal gas constant (in units of kJ/(mol*K))
    R = 8.3144621e-3

    # The temperature and pressure of the system
    T = properties.temperature().val

    # Calculate the saturation index of the mineral
    lnK = reaction_siderite.lnEquilibriumConstant(properties)
    lnQ = reaction_siderite.lnReactionQuotient(properties)
    lnOmega = lnQ.val - lnK.val

    # Calculate the saturation index
    Omega = exp(lnOmega)

    # The composition of the chemical system
    n = properties.composition()

    # The initial composition of the chemical system
    n0 = reaction_siderite.initialAmounts()

    # The index of the mineral
    imineral = system.indexSpeciesWithError(min_reaction_siderite.mineral())

    # The number of moles of the mineral
    nm = n[imineral].val
    nm0 = n0[imineral]

    # Prevent negative mole numbers here for the solution of the ODEs
    nm = max(nm, 0.0)
    nm0 = max(nm0, 0.0)

    # Get H+ activity
    lna = properties.lnActivities().val
    i_h = system.indexSpeciesWithError("H+")
    activity_h = exp(lna[i_h])

    # The molar mass of the mineral (in units of kg/mol)
    molar_mass = system.species(imineral).molarMass()

    # If (m <= 0) and (SRmin < 1) Then GoTo 350
    # If (SRmin = 1) Then GoTo 350
    if ((nm <= 0 and Omega < 1) or (Omega == 1)): # the is no way to precipitate further
        res = 0

    # Average BET value in 0.13 m2/g => 0.13 * 1e3 m2/kg
    ssa = 0.13 * 1e3
    # Molar volume
    molar_volume = 2.049E-05

    # If (SRmin > 1) Then GoTo 130
    if Omega > 1: # precipitation kinetics


        # Interfacial energy in J m-2
        surf_energy = 0.06

        # Molecular volume in m3
        molecular_volume = 3.4e-29

        u = 16 * 3.14159 * pow(surf_energy, 3.0) * pow(molecular_volume, 2.0) / (3 * pow(1.38E-23 * T, 3.0))

        # Nucleation rate in nuclei kg-1 sec-1 (after Fritz et al 2009)
        j0 = 1e20

        # Critical saturation threshold
        Omega_crit = exp(sqrt(u / log10(j0)))


        if Omega < Omega_crit:

            # Nucleation rate in nuclie sec-1 multiplied by time to account for nucleation over time
            nuclie = j0 * exp(- u / pow(log10(Omega), 2.0))

            # Number of growth units in critical nuclie radius
            nj = 2 * u / pow(log10(Omega), 3.0)

            # Critical nuclie volume
            vol = nj * molecular_volume

            # Moles of nuclie formed
            res_nuc += - nuclie * vol / molar_volume

        # Define rate constant of kinetic growth
        kappa_growth = - 1.6E-11 * exp(- 108 / R * (1.0 / T - 1.0 / 298.15))

        # Calculate kinetic rate of the mineral growth
        res_growth += f * ssa * nm * molar_mass * kappa_growth * abs(Omega - 1)

        # Sum the rate of growth and nucleation
        res += res_growth + res_nuc

    else: # dissolution kinetics

        # Define rate constants of mineral dissolution
        kappa_neu = 2.1E-09 * exp(- 56.0 / R * (1.0/T - 1.0/298.15))
        kappa_acid = 5.9E-06 * exp(- 56.0 / R * (1.0/T - 1.0/298.15)) * pow(activity_h, 0.6)
        kappa = kappa_neu + kappa_acid

        # Calculate kinetic rate of mineral dissolution
        res += f * ssa * nm * molar_mass * pow(nm / nm0, 2 / 3) * kappa * (1 - Omega)

    return res
# Set kinetic rate function
reaction_siderite.setRate(rate_func_siderite_shell)
# -

# Define kinetic parameters of the kaolinite:

# +
# Define reaction equation
eq_str_kaolinite = "Kaolinite + 6*H+ = 2*Al+++ + 2*SiO2(aq) + 5*H2O(l)"
# Define mineral reaction
min_reaction_kaolinite = editor.addMineralReaction("Kaolinite") \
    .setEquation(eq_str_kaolinite) \
    .setSpecificSurfaceArea(11.8, "m2/g")
# Create reaction
reaction_kaolinite = createReaction(min_reaction_kaolinite, system)
reaction_kaolinite.setName("Kaolinite kinetics")
# Kinetic rate function
def rate_func_kaolinite_shell(properties):

    # The number of chemical species in the system
    num_species = system.numSpecies()

    # Reaction rate variables
    res = ChemicalScalar(num_species)
    res_growth = ChemicalScalar(num_species)
    res_nuc = ChemicalScalar(num_species)

    # Auxiliary variable for calculating mineral reaction rate
    f = ChemicalScalar(num_species, 1.0)

    # The universal gas constant (in units of kJ/(mol*K))
    R = 8.3144621e-3

    # The temperature and pressure of the system
    T = properties.temperature().val

    # Calculate the saturation index of the mineral
    lnK = reaction_kaolinite.lnEquilibriumConstant(properties)
    lnQ = reaction_kaolinite.lnReactionQuotient(properties)
    lnOmega = lnQ.val - lnK.val

    # Calculate the saturation index
    Omega = exp(lnOmega)

    # The composition of the chemical system
    n = properties.composition()

    # The initial composition of the chemical system
    n0 = reaction_kaolinite.initialAmounts()

    # The index of the mineral
    imineral = system.indexSpeciesWithError(min_reaction_kaolinite.mineral())

    # The number of moles of the mineral
    nm = n[imineral].val
    nm0 = n0[imineral]

    # Prevent negative mole numbers here for the solution of the ODEs
    nm = max(nm, 0.0)
    nm0 = max(nm0, 0.0)

    # Get H+ activity
    lna = properties.lnActivities().val
    i_h = system.indexSpeciesWithError("H+")
    activity_h = exp(lna[i_h])
    i_oh = system.indexSpeciesWithError("OH-")
    activity_oh = exp(lna[i_oh])

    # The molar mass of the mineral (in units of kg/mol)
    molar_mass = system.species(imineral).molarMass()

    # If (m <= 0) and (SRmin < 1) Then GoTo 350
    # If (SRmin = 1) Then GoTo 350
    if ((nm <= 0 and Omega < 1) or (Omega == 1)): # the is no way to precipitate further
        res = 0

    # average BET; suggested value in 11.8 m2/g
    ssa = 11.8 * 1e3
    # Molar volume
    molar_volume = 9.876E-05

    # If (SRmin > 1) Then GoTo 130
    if Omega > 1: # precipitation kinetics

        # Interfacial energy in J m-2
        surf_energy = 0.1

        # Thickness of one mineral layer in m
        sheet_thick = 7e-10

        # Molecular volume in m3
        molecular_volume = 1.64E-28

        u = 2 * sqrt(3) * sheet_thick * pow(surf_energy, 2.0) * molecular_volume / pow(1.38E-23 * T, 2.0)

        # Nucleation rate in nuclei kg-1 sec-1 (after Fritz et al 2009)
        j0 = 1e20

        # Critical saturation threshold
        Omega_crit = exp(u / log10(j0))

        if Omega < Omega_crit:

            # Nucleation rate in nuclie sec-1 multiplied by time to account for nucleation over time
            nuclie = j0 * exp(- u / log10(Omega))

            # Number of growth units in critical nuclie radius
            nj = 2 * u / pow(log10(Omega), 3.0)

            # Critical nuclie volume
            vol = nj * molecular_volume

            # Moles of nuclie formed
            res_nuc += - nuclie * vol / molar_volume

        # Define rate constant of kinetic growth
        kappa_growth = -5.5E-13 * exp(- 66.0 / R * (1.0 / T - 1.0 / 298.15))

        # Calculate kinetic rate of the mineral growth
        res_growth += f * ssa * nm * molar_mass * kappa_growth * pow(abs(pow(Omega, 0.06) - 1), 1.68)

        # Sum the rate of growth and nucleation
        res += res_growth + res_nuc

    else: # dissolution kinetics

        # Define rate constants of mineral dissolution
        kappa_neu = 2.1E-09 * exp(- 56.0 / R * (1.0/T - 1.0/298.15))
        kappa_acid = 5.9E-06 * exp(- 56.0 / R * (1.0/T - 1.0/298.15)) * pow(activity_h, 0.6)
        kappa_oh = 2.5E-11 * exp(- 46.0 / R * (1.0/T - 1.0/298.15)) * pow(activity_oh, 0.58)
        kappa = kappa_neu + kappa_acid + kappa_oh

        # Calculate kinetic rate of mineral dissolution
        res += f * ssa * nm * molar_mass * kappa * (1 - Omega)

    return res
# Set kinetic rate function
reaction_kaolinite.setRate(rate_func_kaolinite_shell)
# -

# +
# Define reaction equation
eq_str_quartz = "Quartz + H2O(l) = H+ + HSiO3-"
# Define mineral reaction
min_reaction_quartz = editor.addMineralReaction("Quartz") \
    .setEquation(eq_str_quartz) \
    .setSpecificSurfaceArea(0.1, "m2/g")
# Create reaction
reaction_quartz = createReaction(min_reaction_quartz, system)
reaction_quartz.setName("Quartz kinetics")
# Kinetic rate function
def rate_func_quartz_shell(properties):

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
    lnK = reaction_quartz.lnEquilibriumConstant(properties)
    lnQ = reaction_quartz.lnReactionQuotient(properties)
    lnOmega = lnQ.val - lnK.val

    # Calculate the saturation index
    Omega = exp(lnOmega)

    # The composition of the chemical system
    n = properties.composition()

    # The initial composition of the chemical system
    n0 = reaction_quartz.initialAmounts()

    # The index of the mineral
    imineral = system.indexSpeciesWithError(min_reaction_quartz.mineral())

    # The number of moles of the mineral
    nm = n[imineral].val
    nm0 = n0[imineral]

    # Prevent negative mole numbers here for the solution of the ODEs
    nm = max(nm, 0.0)
    nm0 = max(nm0, 0.0)

    # Get OH- activity
    lna = properties.lnActivities().val
    i_oh = system.indexSpeciesWithError("OH-")
    activity_oh = exp(lna[i_oh])

    # The molar mass of the mineral (in units of kg/mol)
    molar_mass = system.species(imineral).molarMass()

    # If (m <= 0) and (SRmin < 1) Then GoTo 240
    # If (SRmin = 1) Then GoTo 240
    if ((nm <= 0 and Omega < 1) or (Omega == 1)): # the is no way to precipitate further
        res = 0
    ssa = 0.006

    # If (SRmin > 1) Then GoTo 130
    if Omega > 1: # precipitation kinetics

        if nm > 1e-8:
            kappa_pre = -3.24e-12
            res += f * ssa * nm * molar_mass * kappa_pre * pow(abs(pow(Omega, 4.58) - 1), 0.54)
        else:
            # Set nucleation rate
            res = -1e-10 * f

    else: # dissolution kinetics

        # Define rate constants of mineral dissolution
        kappa_neu = 6.42E-14 * exp(- 76.7 / R * (1.0/T - 1.0/298.15))
        kappa_oh = 0.000000000192 * exp(- 80.0 / R * (1.0/T - 1.0/298.15)) * pow(activity_oh, 0.339)

        # Sum the neutral and oh-promoted rates
        kappa = kappa_neu + kappa_oh

        # Calculate kinetic rate of mineral dissolution
        res += f * ssa * nm * molar_mass * pow(nm / nm0, 2 / 3) * kappa * (1 - Omega)

    return res
# Set kinetic rate function
reaction_quartz.setRate(rate_func_quartz_shell)
# -

# Define the reaction system
reactions = ReactionSystem(system, [reaction_calcite, reaction_siderite, reaction_daphnite, reaction_kaolinite, reaction_quartz])

# Specifying the partition including the kinetic species:
partition = Partition(system)
partition.setKineticPhases(["Calcite", "Siderite", "Daphnite,14A", "Kaolinite", "Quartz"])

# Define thermodynamic parameters
T = 25.0    # temperature (in units of celcius)
P = 1.0   # pressure (in units of atm)
water_kg = 1.00 # amount of water in brine

# Define the initial state corresponding to the formation water:

# +
problem_ic = EquilibriumInverseProblem(system)
problem_ic.setTemperature(T, "celsius")
problem_ic.setPressure(P, "atm")
problem_ic.add("H2O", water_kg, "kg")
problem_ic.add("HCO3-", 1441.11, "mg") # 1441.11 mg/l
problem_ic.add("Ca++", 310, "mg") # 310 mg/l
problem_ic.add("Cl-", 19100, "mg") # 19100 charge
problem_ic.add("Fe++", 0.0034, "mg") # 0.0034 mg/l
problem_ic.add("Mg++", 1240, "mg") # 1240 mg/l
problem_ic.add("Na+", 11400, "mg") # 11400 mg/l
problem_ic.add("SO4--", 2420, "mg") # 2420 mg/l
# Sulfide spike
problem_ic.add("S", 100, "mg") # 100 mg/l
# Equilibrium phases
problem_ic.add("Calcite", 0.0, "mol") # Calcite 0 0
problem_ic.add("Pyrrhotite", 0.0, "mol") # Pyrrhotite 0 0
problem_ic.pH(7.1015)
problem_ic.pE(-3.6887)

# Calculate the equilibrium states for the initial conditions
state_ic = equilibrate(problem_ic)
state_ic.setSpeciesAmount("Calcite", 0.1, "mol") # MM(100.09) = 713.5 g / mol
state_ic.setSpeciesAmount("Daphnite,14A", 0.1, "mol") # MM(Daphnite) = 713.5 g / mol
state_ic.setSpeciesAmount("Siderite", 0.1, "mol") # MM(Siderite) = 115.86 g / mol
state_ic.setSpeciesAmount("Kaolinite", 0.1, "mol") # MM(Kaolinite) = 258.071 g / mol
state_ic.setSpeciesAmount("Quartz", 0.1, "mol") # MM(Quartz) = 60.083 g / mol

state_ic.scaleVolume(1.0, "m3")

reaction_calcite.setInitialAmounts(state_ic.speciesAmounts())
reaction_daphnite.setInitialAmounts(state_ic.speciesAmounts())
reaction_siderite.setInitialAmounts(state_ic.speciesAmounts())
reaction_kaolinite.setInitialAmounts(state_ic.speciesAmounts())
reaction_quartz.setInitialAmounts(state_ic.speciesAmounts())
# -

# Performing the kinetic path calculation:

# Define KineticPath instance to perform kinetic simulations
path = KineticPath(reactions, partition)

# To analyse the result of kinetic simulations, we save the evolution of different properties of the chemical system
# into file `result_file_name`:

result_file_name = "kineticpath-barite-kinetics-n-steps.txt"
output = path.output()
output.add("time(units=s)")
output.add("pH")
output.add("speciesAmount(Ca++)")
output.add("speciesAmount(Mg++)")
output.add("speciesAmount(Na+)")
output.add("speciesAmount(Fe++)")
output.add("speciesAmount(Fe+++)")
output.add("speciesAmount(SO4--)")
output.add("speciesAmount(HS-)")
output.add("speciesAmount(HCO3-)")
output.add("speciesAmount(H+")
output.add("speciesAmount(Al+++")
output.add("speciesAmount(SiO2(aq))")
output.add("speciesAmount(HSiO3-")
output.add("speciesAmount(Pyrrhotite)")
output.add("speciesAmount(Calcite)")
output.add("speciesAmount(Quartz)")
output.add("speciesAmount(Kaolinite)")
output.add("speciesAmount(Daphnite,14A)")
output.add("speciesAmount(Siderite)")

# Solving the chemical kinetics problem on n uniform steps:

t0 = 0
dt = 3600.0
n = 1440

result_file_name = "kineticpath-complex-scavenging-" + str(n) + "-steps.txt"
output.filename(result_file_name)

path.solve(state_ic, t0, dt, n, "second")

# For plotting of the results of equilibrium path calculation, we load the results into the `data` array:

# +
import matplotlib.pyplot as plt
import numpy as np
results_file = result_file_name

# Load the Reaktoro values from the file
data_reaktoro_pyrrhotite = np.loadtxt(results_file, skiprows=1, usecols=(-6))
data_reaktoro_calcite = np.loadtxt(results_file, skiprows=1, usecols=(-5))
data_reaktoro_quartz = np.loadtxt(results_file, skiprows=1, usecols=(-4))
data_reaktoro_kaolinite = np.loadtxt(results_file, skiprows=1, usecols=(-3))
data_reaktoro_daphnite = np.loadtxt(results_file, skiprows=1, usecols=(-2))
data_reaktoro_siderite = np.loadtxt(results_file, skiprows=1, usecols=(-1))
data_reaktoro_time = np.loadtxt(results_file, skiprows=1, usecols=(0))

fig, (ax1, ax2) = plt.subplots(2, 1)
ax1.set_xlim([data_reaktoro_time[0] - 2, data_reaktoro_time[-1] + 2])
ax1.set_ylim([0.0875, 0.101])
ax1.set_ylabel('Amount [mol]', fontproperties=prop)
ax1.plot(data_reaktoro_time, data_reaktoro_siderite, label='Siderite', color="C2")
ax1.legend(loc='center right', prop=prop)
ax2.set_xlim([data_reaktoro_time[0] - 2, data_reaktoro_time[-1] + 2])
ax2.set_ylim([-0.001, 0.0125])
ax2.set_ylabel('Amount [mol]', fontproperties=prop)
ax2.plot(data_reaktoro_time, data_reaktoro_pyrrhotite, label='Pyrrhotite', color="C1")
ax2.legend(loc='center right', prop=prop)
ax2.set_title("", fontproperties=prop, fontsize=20)
plt.xlabel('Time [s]', fontproperties=prop)
plt.savefig('reaktoro-pyrrhotite-siderite.png')
plt.close()


fig = plt.figure()
plt.axes(xlim=(data_reaktoro_time[0] - 2, data_reaktoro_time[-1] + 2))
plt.xlabel('Time [s]', fontproperties=prop)
plt.ylabel('Amount [mol]', fontproperties=prop)
plt.plot(data_reaktoro_time, data_reaktoro_daphnite, label='Daphnite', color="C3")
plt.plot(data_reaktoro_time, data_reaktoro_kaolinite, label='Kaolinite', color="C4")
plt.plot(data_reaktoro_time, data_reaktoro_quartz, label='Quartz', color="C5")
plt.plot(data_reaktoro_time, data_reaktoro_calcite, label='Calcite', color="C6")
plt.legend(loc='upper left', prop=prop)
plt.savefig('reaktoro-daphnite-kaolinite-quartz-calcite.png')
# -