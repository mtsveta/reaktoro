from reaktoro import *
import thermofun.PyThermoFun as thermofun

database = thermofun.Database("databases/thermofun/aq17-thermofun.json")

T = 293.15
P = 100000.0

editor = ChemicalEditor(database)
editor.setTemperatures(T, "Kelvin")
editor.setPressures(P, "Pascal")

editor.addAqueousPhase([
    "Al(OH)2+", "Al(OH)3@", "Al(OH)4-", "Al+3", "AlH3SiO4+2", "AlOH+2",
    "Ca+2", "CaCO3@", "CaCl+", "CaCl2@", "CaHCO3+", "CaHSiO3+", "CaOH+", "CaSiO3@", "K+", "KAlO2@",
    "KCl@", "KOH@", "KCO3-", "KHCO3@", "Mg+2", "MgCO3@", "MgCl+", "MgCl2@", "MgHCO3+", "MgHSiO3+",
    "MgOH+", "MgSiO3@", "Na+", "NaAl(OH)4@", "NaCO3-", "NaCl@", "NaHCO3@", "NaHSiO3@",
    "NaOH@", "HSiO3-", "SiO2@", "CO@", "CO2@", "CO3-2", "HCO3-", "CH4@", "Cl-",
    "HCl@", "H2@", "O2@", "OH-", "H+", "H2O@"])

'''
## (5) Data for Dependent Components
'Al(SO4)+' 'Al(SO4)2-' 'Al+3' 'AlO+' 'AlO2-' 'AlO2H@' 'AlOH+2' 'AlHSiO3+2' 'AlSiO5-3' 'Ca(CO3)@' 'Ca(HCO3)+' 'Ca(SO4)@' 'Ca+2' 'CaOH+' 'Ca(HSiO3)+' 'CaSiO3@' 'Fe(CO3)@' 'Fe(HCO3)+' 'Fe(HSO4)+' 'Fe(SO4)@' 'Fe+2' 'FeCl+' 'FeOH+' 'Fe(HSO4)+2' 'Fe(SO4)+' 'Fe(SO4)2-' 'Fe+3' 'Fe2(OH)2+4' 'Fe3(OH)4+5' 'FeCl+2' 'FeCl2+' 'FeCl3@' 'FeO+' 'FeO2-' 'FeO2H@' 'FeOH+2' 'FeHSiO3+2' 'K(SO4)-' 'K+' 'KOH@'
'Mg(CO3)@' 'Mg(HCO3)+' 'Mg+2' 'MgOH+' 'MgSO4@' 'Mg(HSiO3)+' 'MgSiO3@' 'Na(CO3)-' 'Na(HCO3)@' 'Na(SO4)-' 'Na+' 'NaOH@' 'HSiO3-' 'Si4O10-4' 'SiO2@' 'SiO3-2' 'CO2@' 'CO3-2' 'HCO3-' 'CH4@' 'ClO4-' 'Cl-' 'H2@' 'N2@' 'O2@' 'S2O3-2' 'HSO3-' 'SO3-2' 'HSO4-' 'SO4-2' 'H2S@' 'HS-' 'S-2' 'OH-' 'H+' 'H2O@' 
'CO2' 'CH4' 'H2' 'N2' 'O2' 'H2S' 'H2O' 'monocarbonate' 'C2AClH5' 'C3AFS0.84H4.32' 'C3FS0.84H4.32' 'CSHQ-JenD' 'CSHQ-JenH' 'CSHQ-TobD' 'CSHQ-TobH' 'KSiOH' 'NaSiOH' 'straetlingite' 'straetlingite7' 'C4AH13' 'monosulphate12' 'tricarboalu03' 'ettringite03_ss' 'M075SH' 'M15SH' 'AlOHmic' 'Kln' 'Gr' 'C12A7' 'C2S' 'C3A' 'C3S' 'C4AF' 'CA' 'CA2' 'C2AH7.5' 'C3AH6' 'C4AH11' 'C4AH13' 'C4AH19' 'CAH10' 'monosulphate10.5' 'monosulphate12' 'monosulphate14'
'monosulphate16' 'monosulphate9' 'chabazite' 'zeoliteP_Ca' 'straetlingite5.5' 'monocarbonate9' 'hemicarbonat10.5' 'hemicarbonate' 'hemicarbonate9' 'monocarbonate' 'C4AClH10' 'C4AsClH12' 'ettringite13' 'ettringite9' 'Arg' 'Cal' 'C3FH6' 'C4FH13' 'C3FS0.84H4.32' 'C3FS1.34H3.32' 'Fe-hemicarbonate' 'Femonocarbonate' 'Dis-Dol' 'Ord-Dol' 'Lim' 'Portlandite' 'Anh' 'Gp' 'hemihydrate' 'thaumasite' 'Fe' 'FeCO3(pr)' 'Sd' 'Mag' 'Fe(OH)3(am)' 'FeOOHmic' 'Py' 'Tro' 'Melanterite' 'K2SO4'
'syngenite' 'K2O' 'hydrotalcite' 'Mgs' 'Brc' 'Na2SO4' 'natrolite' 'zeoliteX' 'zeoliteY' 'Na2O' 'Sulfur' 'Qtz' 'Amor-Sl'

# xDC: Speciation - amounts of DCs in equilibrium (primal solution), moles (GEM output/input) [nDCb]
<xDC>
7.08390205407382e-23 3.05792140115698e-30 9.46429572449386e-16 3.38222265993579e-14 1.16939763075074e-12 2.25272922927293e-13 7.89337370586244e-15 9.27837865538291e-21 2.92542199589053e-28 2.4941732333023e-21 7.44386960014572e-20 1.84774555939163e-20 1.4373321272283e-12 3.70236949679406e-19 2.93246636232349e-23 2.03906048121015e-26 1.97089343351131e-28 3.23458112235902e-27 7.82524797623901e-35 8.13304052568418e-29 7.28469372559797e-21 1.03067270957778e-21 4.24726410837069e-24 3.59392121966752e-32 2.72617255684425e-26 2.50890769361049e-34 2.62393470512793e-19 6.15265151833959e-24 2.06779162516089e-28 3.18992129323225e-19 1.53978234556082e-19 4.30494795195437e-21 9.15847575258294e-13 6.11763364772575e-15 5.1403007957578e-13 1.33655574185217e-15 6.8053150064198e-22 4.12700321827989e-21 1.43733256763097e-12 2.11469419530701e-20
1.49375478413927e-21 6.79050879060711e-20 1.43732347098944e-12 9.03112256091319e-18 2.13334759017737e-20 6.03825917427575e-23 2.14274044876232e-25 3.75831596989159e-14 1.35902158952379e-12 4.83202429846381e-13 0.000238186675912991 7.05055407480896e-12 2.98958603678772e-15 0 1.43433453497021e-12 4.29224068516673e-21 1.38112208616446e-12 1.66762916903346e-14 9.39850682283554e-12 0 9.89464075658034e-30 0.000238186684843353 0 2.67444707581408e-08 6.97817229238136e-09 0 0 0 1.92970529701553e-18 9.54128189344437e-13 0 0 0 4.1403441039378e-11 5.52155873495934e-11 0.0240328647467665 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0.1078271 0

'''
'''
## (6) Data for Phases
#  'aq_gen' 'gas_gen' 'Friedels-ss' 'C3(AF)S0.84H' 'CSHQ' 'straetlingite' 'SO4_OH_AFm' 'SO4_CO3_AFt' 'MSH' 'Al(OH)3mic' 'Kaolinite' 'Graphite' 'Mayenite' 'Belite' 'Aluminate' 'Alite' 'Ferrite' 'CA' 'CA2' 'C2AH75' 'C3AH6' 'C4AH11' 'C4AH13' 'C4AH19' 'CAH10' 'C4AsH105' 'C4AsH12' 'C4AsH14' 'C4AsH16' 'C4AsH9' 'Chabazite' 'ZeoliteP' 'C2ASH55' 'C4AcH9' 'C4Ac0.5H105' 'C4Ac0.5H12' 'C4Ac0.5H9' 'C4AcH11' 'Friedels' 'Kuzels'
#  'C6AsH13' 'C6AsH9' 'Aragonite' 'Calcite' 'C3FH6' 'C4FH13' 'C3FS0.84H4.32' 'C3FS1.34H3.32' 'C4Fc05H10' 'C4FcH12' 'Dolomite-dis' 'Dolomite-ord' 'lime' 'Portlandite' 'Anhydrite' 'Gypsum' 'hemihydrate' 'thaumasite' 'Iron' 'Fe-carbonate' 'Siderite' 'Magnetite' 'Ferrihydrite-am' 'Ferrihydrite-mc' 'Pyrite' 'Troilite' 'Melanterite' 'arcanite' 'syngenite' 'K-oxide' 'OH-hydrotalcite' 'Magnesite' 'Brucite' 'thenardite' 'Natrolite' 'ZeoliteX' 'ZeoliteY' 'Na-oxide' 'Sulphur' 'Quartz'
#  'Silica-amorph'
# aPH: Specific surface areas of phases, m2/kg (GEM input) [nPHb]
'''
editor.addMineralPhase("Albite")
editor.addMineralPhase("Andalusite")
editor.addMineralPhase("Calcite")
editor.addMineralPhase("Corundum")
editor.addMineralPhase("Diopside")
editor.addMineralPhase("Dolomite")
editor.addMineralPhase("Enstatite")
editor.addMineralPhase("Grossular")
editor.addMineralPhase("Margarite")
editor.addMineralPhase("Microcline")
editor.addMineralPhase("Muscovite")
editor.addMineralPhase("Pargasite-Mg")
editor.addMineralPhase("Phlogopite")
editor.addMineralPhase("Quartz")
editor.addMineralPhase("Sanidine")
editor.addMineralPhase("Sillimanite")
editor.addMineralPhase("Zoisite")

system = ChemicalSystem(editor)

problem = EquilibriumProblem(system)
problem.add("H2O",             	1000,	"g")
problem.add("CO2",             	0.001,  "g")
problem.add("CaCO3",           	1,	    "g")
problem.add("MgSiO3",          	1,	    "g")
problem.add("NaCl",            	5,      "g")
problem.add("NaAlSi3O8",        37,	    "g")
problem.add("KAl3Si3O10(OH)2",  13,	    "g")
problem.add("SiO2",          	30,	    "g")
problem.add("KAlSi3O8",        	20,	    "g")
problem.setTemperature(T, "Kelvin")
problem.setPressure(P, "Pascal")

# Step 7.2: Create the ChemicalSystem object using the configured editor
system = ChemicalSystem(editor)

# Step 7.3: Define the initial condition of the reactive transport modeling problem
problem_ic = EquilibriumProblem(system)
problem_ic.setTemperature(T)
problem_ic.setPressure(P)
problem_ic.add('H2O', 1.0, 'kg')
problem_ic.add('NaCl', 0.7, 'mol')
problem_ic.add('CaCO3', 10, 'mol')
problem_ic.add('SiO2', 10, 'mol')

# Step 7.4: Define the boundary condition of the reactive transport modeling problem
problem_bc = EquilibriumProblem(system)
problem_bc.setTemperature(T)
problem_bc.setPressure(P)
problem_bc.add('H2O', 1.0, 'kg')
problem_bc.add('NaCl', 0.90, 'mol')
problem_bc.add('MgCl2', 0.05, 'mol')
problem_bc.add('CaCl2', 0.01, 'mol')
problem_bc.add('CO2', 0.75, 'mol')

options = EquilibriumOptions()
options.smart.reltol = 0.1
options.smart.abstol = 0.1

# Step 7.5: Calculate the equilibrium states for the initial and boundary conditions
state_ic = equilibrate(problem_ic)
state_bc = equilibrate(problem_bc)

# Step 7.6: Scale the volumes of the phases in the initial condition
state_ic.scalePhaseVolume('Aqueous', 0.1, 'm3')
state_ic.scalePhaseVolume('Quartz', 0.882, 'm3')
state_ic.scalePhaseVolume('Calcite', 0.018, 'm3')

# Scale the boundary condition state to 1 m3
state_bc.scaleVolume(1.0, 'm3')