import numpy as np
import matplotlib.pyplot as plt
import os

from matplotlib import font_manager as fm, rcParams
fpath = os.path.join(rcParams["datapath"], "texgyreadventor-regular.otf")
prop = fm.FontProperties(fname=fpath)
prop.set_size(14)

import matplotlib as mpl
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'TeX Gyre Adventor'
mpl.rcParams['font.style'] = 'normal'
mpl.rcParams['font.size'] = 14
mpl.set_loglevel("critical")

# Plotting params
circ_area = 6 ** 2
custom_font = { }

C0 = '#107ab0'
C1 = '#fc5a50'

T_left = 25
T_right = 90
temperature_points = 14
temperatures = np.linspace(T_left, T_right, temperature_points)  # the x-coordinates of the plots

def empty_marker(color):
    return {'facecolor': 'white', 'edgecolor': color, 's': circ_area, 'zorder': 2, 'linewidths': 1.5 }

def filled_marker(color):
    return {'color': color, 's': circ_area, 'zorder': 2, 'linewidths': 1.5 }

def line_empty_marker(color):
    return {'marker': 'd', 'markerfacecolor': 'white', 'markeredgecolor':color, 'markersize': 6, 'markeredgewidth': 1.5 }

def line_filled_marker(color):
    return {'color': color, 'markersize': 6, 'markeredgewidth': 1.5 }

def line(color):
    return {'linestyle': '-', 'color': color, 'zorder': 1, 'linewidth': 2}

def line_error(color):
    return {'linestyle': ':', 'marker': 'D', 'color': color, 'zorder': 1, 'linewidth': 0.5, 'markersize': 4}

def plot_concentrations(data_reaktoro, data_phreeqc, tag, title, folder):

    colors = ['C2', 'C3']
    plt.axes(xlim=(temperatures[0] - 2, temperatures[-1] + 2))

    plt.xlabel(r'Temperatures [C$\degree$]', fontproperties=prop)
    plt.ylabel('Concentration [mol/kgw]', fontproperties=prop)
    plt.title(title, fontproperties=prop)

    plt.plot(temperatures, data_reaktoro[0], label='Calcite, Reaktoro', **line(colors[0]))
    plt.plot(temperatures, data_phreeqc[0], 's', **line_filled_marker(colors[0]))[0],
    plt.plot([], [], 's', label='Calcite, PHREEQC', **line_filled_marker(colors[0]))

    plt.plot(temperatures, data_reaktoro[1], label='Dolomite, Reaktoro', **line(colors[1]))
    plt.plot(temperatures, data_phreeqc[1], 'o', **line_filled_marker(colors[1]))[0],
    plt.plot([], [], 'o', label='Dolomite, PHREEQC', **line_filled_marker(colors[1]))

    plt.legend(loc='center right', prop=prop)
    plt.savefig(folder + '/reaktoro-phreeqc-carbonates-seawater' + tag +'.pdf')
    plt.close()

def plot_aqueous_concentrations(data_reaktoro, data_phreeqc, tag, title, folder):

    colors = ['C1', 'C4', 'C5', 'C6']
    plt.axes(xlim=(temperatures[0] - 2, temperatures[-1] + 2))

    plt.xlabel(r'Temperature [$\degree$C]', fontproperties=prop)
    plt.ylabel('Concentration [mol/kgw]', fontproperties=prop)
    plt.title(title, fontproperties=prop)
    plt.yscale("log")

    plt.plot(temperatures, data_reaktoro[0], label=r'Ca$^{2+}$, Reaktoro', **line(colors[0]))
    plt.plot(temperatures, data_phreeqc[0], 's', **line_filled_marker(colors[0]))[0],
    plt.plot([], [], 's', label='Ca$^{2+}$, PHREEQC', **line_filled_marker(colors[0]))

    plt.plot(temperatures, data_reaktoro[1], label='Mg$^{2+}$, Reaktoro', **line(colors[1]))
    plt.plot(temperatures, data_phreeqc[1], 'o', **line_filled_marker(colors[1]))[0],
    plt.plot([], [], 'o', label='Mg$^{2+}$, PHREEQC', **line_filled_marker(colors[1]))

    plt.plot(temperatures, data_reaktoro[2], label='H$^{+}$, Reaktoro', **line(colors[2]))
    plt.plot(temperatures, data_phreeqc[2], 's', **line_filled_marker(colors[2]))[0],
    plt.plot([], [], 's', label='H$^{+}$, PHREEQC', **line_filled_marker(colors[2]))

    plt.plot(temperatures, data_reaktoro[3], label='HCO$_3^{-}$, Reaktoro', **line(colors[3]))
    plt.plot(temperatures, data_phreeqc[3], 'o', **line_filled_marker(colors[3]))[0],
    plt.plot([], [], 'o', label='HCO$_3^{-}$, PHREEQC', **line_filled_marker(colors[3]))

    plt.legend(loc='center right', prop=prop)
    plt.savefig(folder + '/reaktoro-phreeqc-carbonates-seawater-aqueous' + tag +'.pdf')
    plt.close()

if __name__ == '__main__':

    results_folder = 'results-carbonates-seawater'

    # ---------------------------------------------------------------------------------------------------------------- #
    # P = 1.0 atm
    # ---------------------------------------------------------------------------------------------------------------- #

    data_reaktoro_calcite = np.loadtxt(results_folder + '/reaktoro-seawater-phreeqc-Calcite-p-1.0.txt')
    data_reaktoro_dolomite = np.loadtxt(results_folder + '/reaktoro-seawater-phreeqc-Dolomite-p-1.0.txt')
    data_reaktoro = np.vstack((data_reaktoro_calcite.T, data_reaktoro_dolomite.T))

    # Load from the file -d_CO2(g) concentrations
    data_phreeqc_calcite = np.loadtxt(results_folder + '/phreeqc-seawater-with-carbonates-p-1.txt', skiprows=2, usecols=(-4))
    data_phreeqc_dolomite = np.loadtxt(results_folder + '/phreeqc-seawater-with-carbonates-p-1.txt', skiprows=2, usecols=(-2))
    data_phreeqc = np.vstack((data_phreeqc_calcite, data_phreeqc_dolomite))

    title = r'Minerals concentrations using phreeqc.dat, P = 1 atm'
    plot_concentrations(data_reaktoro, data_phreeqc, '-phreeqc-dat-P-1-atm', title, results_folder)

    data_reaktoro_ca2 = np.loadtxt(results_folder + '/reaktoro-seawater-phreeqc-Ca2-p-1.0.txt')
    data_reaktoro_mg2 = np.loadtxt(results_folder + '/reaktoro-seawater-phreeqc-Mg2-p-1.0.txt')
    data_reaktoro_h = np.loadtxt(results_folder + '/reaktoro-seawater-phreeqc-H-p-1.0.txt')
    data_reaktoro_hco3 = np.loadtxt(results_folder + '/reaktoro-seawater-phreeqc-HCO3-p-1.0.txt')
    data_reaktoro = np.vstack((data_reaktoro_ca2.T, data_reaktoro_mg2.T, data_reaktoro_h.T, data_reaktoro_hco3.T))

    # Load from the file -d_CO2(g) concentrations
    data_phreeqc_ca2 = np.loadtxt(results_folder + '/phreeqc-seawater-with-carbonates-p-1.txt', skiprows=2, usecols=(-8))
    data_phreeqc_mg2 = np.loadtxt(results_folder + '/phreeqc-seawater-with-carbonates-p-1.txt', skiprows=2, usecols=(-7))
    data_phreeqc_h = np.loadtxt(results_folder + '/phreeqc-seawater-with-carbonates-p-1.txt', skiprows=2, usecols=(-6))
    data_phreeqc_hco3 = np.loadtxt(results_folder + '/phreeqc-seawater-with-carbonates-p-1.txt', skiprows=2, usecols=(-5))
    data_phreeqc = np.vstack((data_phreeqc_ca2, data_phreeqc_mg2, data_phreeqc_h, data_phreeqc_hco3))

    title = r'Aqueous species concentrations using phreeqc.dat, P = 1 atm'
    plot_aqueous_concentrations(data_reaktoro, data_phreeqc, '-phreeqc-dat-P-1-atm', title, results_folder)

    # ---------------------------------------------------------------------------------------------------------------- #
    # P = 100.0 atm
    # ---------------------------------------------------------------------------------------------------------------- #

    data_reaktoro_calcite = np.loadtxt(results_folder + '/reaktoro-seawater-phreeqc-Calcite-p-100.0.txt')
    data_reaktoro_dolomite = np.loadtxt(results_folder + '/reaktoro-seawater-phreeqc-Dolomite-p-100.0.txt')
    data_reaktoro = np.vstack((data_reaktoro_calcite.T, data_reaktoro_dolomite.T))

    # Load from the file -d_CO2(g) concentrations
    data_phreeqc_calcite = np.loadtxt(results_folder + '/phreeqc-seawater-with-carbonates-p-100.txt', skiprows=2, usecols=(-4))
    data_phreeqc_dolomite = np.loadtxt(results_folder + '/phreeqc-seawater-with-carbonates-p-100.txt', skiprows=2, usecols=(-2))
    data_phreeqc = np.vstack((data_phreeqc_calcite, data_phreeqc_dolomite))

    title = r'Minerals concentrations using phreeqc.dat, P = 1 atm'
    plot_concentrations(data_reaktoro, data_phreeqc, '-phreeqc-dat-P-100-atm', title, results_folder)

    data_reaktoro_ca2 = np.loadtxt(results_folder + '/reaktoro-seawater-phreeqc-Ca2-p-100.0.txt')
    data_reaktoro_mg2 = np.loadtxt(results_folder + '/reaktoro-seawater-phreeqc-Mg2-p-100.0.txt')
    data_reaktoro_h = np.loadtxt(results_folder + '/reaktoro-seawater-phreeqc-H-p-100.0.txt')
    data_reaktoro_hco3 = np.loadtxt(results_folder + '/reaktoro-seawater-phreeqc-HCO3-p-100.0.txt')
    data_reaktoro = np.vstack((data_reaktoro_ca2.T, data_reaktoro_mg2.T, data_reaktoro_h.T, data_reaktoro_hco3.T))

    # Load from the file -d_CO2(g) concentrations
    data_phreeqc_ca2 = np.loadtxt(results_folder + '/phreeqc-seawater-with-carbonates-p-100.txt', skiprows=2, usecols=(-8))
    data_phreeqc_mg2 = np.loadtxt(results_folder + '/phreeqc-seawater-with-carbonates-p-100.txt', skiprows=2, usecols=(-7))
    data_phreeqc_h = np.loadtxt(results_folder + '/phreeqc-seawater-with-carbonates-p-100.txt', skiprows=2, usecols=(-6))
    data_phreeqc_hco3 = np.loadtxt(results_folder + '/phreeqc-seawater-with-carbonates-p-100.txt', skiprows=2, usecols=(-5))
    data_phreeqc = np.vstack((data_phreeqc_ca2, data_phreeqc_mg2, data_phreeqc_h, data_phreeqc_hco3))

    title = r'Aqueous species concentrations using phreeqc.dat, P = 1 atm'
    plot_aqueous_concentrations(data_reaktoro, data_phreeqc, '-phreeqc-dat-P-100-atm', title, results_folder)

    # ---------------------------------------------------------------------------------------------------------------- #
    # Pitzer
    # ---------------------------------------------------------------------------------------------------------------- #
    # P = 1.0 atm
    # ---------------------------------------------------------------------------------------------------------------- #

    data_reaktoro_calcite = np.loadtxt(results_folder + '/reaktoro-seawater-pitzer-Calcite-p-1.0.txt')
    data_reaktoro_dolomite = np.loadtxt(results_folder + '/reaktoro-seawater-pitzer-Dolomite-p-1.0.txt')
    data_reaktoro = np.vstack((data_reaktoro_calcite.T, data_reaktoro_dolomite.T))

    # Load from the file -d_CO2(g) concentrations
    data_phreeqc_calcite = np.loadtxt(results_folder + '/phreeqc-pitzer-seawater-with-carbonates-p-1.txt', skiprows=2, usecols=(-4))
    data_phreeqc_dolomite = np.loadtxt(results_folder + '/phreeqc-pitzer-seawater-with-carbonates-p-1.txt', skiprows=2, usecols=(-2))
    data_phreeqc = np.vstack((data_phreeqc_calcite, data_phreeqc_dolomite))

    title = r'Minerals concentrations using pitzer.dat, P = 1 atm'
    plot_concentrations(data_reaktoro, data_phreeqc, '-pitzer-dat-P-1-atm', title, results_folder)

    data_reaktoro_ca2 = np.loadtxt(results_folder + '/reaktoro-seawater-pitzer-Ca2-p-1.0.txt')
    data_reaktoro_mg2 = np.loadtxt(results_folder + '/reaktoro-seawater-pitzer-Mg2-p-1.0.txt')
    data_reaktoro_h = np.loadtxt(results_folder + '/reaktoro-seawater-pitzer-H-p-1.0.txt')
    data_reaktoro_hco3 = np.loadtxt(results_folder + '/reaktoro-seawater-pitzer-HCO3-p-1.0.txt')
    data_reaktoro = np.vstack((data_reaktoro_ca2.T, data_reaktoro_mg2.T, data_reaktoro_h.T, data_reaktoro_hco3.T))

    # Load from the file -d_CO2(g) concentrations
    data_phreeqc_ca2 = np.loadtxt(results_folder + '/phreeqc-pitzer-seawater-with-carbonates-p-1.txt', skiprows=2, usecols=(-8))
    data_phreeqc_mg2 = np.loadtxt(results_folder + '/phreeqc-pitzer-seawater-with-carbonates-p-1.txt', skiprows=2, usecols=(-7))
    data_phreeqc_h = np.loadtxt(results_folder + '/phreeqc-pitzer-seawater-with-carbonates-p-1.txt', skiprows=2, usecols=(-6))
    data_phreeqc_hco3 = np.loadtxt(results_folder + '/phreeqc-pitzer-seawater-with-carbonates-p-1.txt', skiprows=2, usecols=(-5))
    data_phreeqc = np.vstack((data_phreeqc_ca2, data_phreeqc_mg2, data_phreeqc_h, data_phreeqc_hco3))
    #print(data_phreeqc)
    #input()
    title = r'Aqueous species concentrations using pitzer.dat, P = 1 atm'
    plot_aqueous_concentrations(data_reaktoro, data_phreeqc, '-pitzer-dat-P-1-atm', title, results_folder)

    # ---------------------------------------------------------------------------------------------------------------- #
    # P = 100.0 atm
    # ---------------------------------------------------------------------------------------------------------------- #

    data_reaktoro_calcite = np.loadtxt(results_folder + '/reaktoro-seawater-pitzer-Calcite-p-100.0.txt')
    data_reaktoro_dolomite = np.loadtxt(results_folder + '/reaktoro-seawater-pitzer-Dolomite-p-100.0.txt')
    data_reaktoro = np.vstack((data_reaktoro_calcite.T, data_reaktoro_dolomite.T))

    # Load from the file -d_CO2(g) concentrations
    data_phreeqc_calcite = np.loadtxt(results_folder + '/phreeqc-pitzer-seawater-with-carbonates-p-100.txt', skiprows=2, usecols=(-4))
    data_phreeqc_dolomite = np.loadtxt(results_folder + '/phreeqc-pitzer-seawater-with-carbonates-p-100.txt', skiprows=2, usecols=(-2))
    data_phreeqc = np.vstack((data_phreeqc_calcite, data_phreeqc_dolomite))

    title = r'Minerals concentrations using pitzer.dat, P = 1 atm'
    plot_concentrations(data_reaktoro, data_phreeqc, '-pitzer-dat-P-100-atm', title, results_folder)

    data_reaktoro_ca2 = np.loadtxt(results_folder + '/reaktoro-seawater-pitzer-Ca2-p-100.0.txt')
    data_reaktoro_mg2 = np.loadtxt(results_folder + '/reaktoro-seawater-pitzer-Mg2-p-100.0.txt')
    data_reaktoro_h = np.loadtxt(results_folder + '/reaktoro-seawater-pitzer-H-p-100.0.txt')
    data_reaktoro_hco3 = np.loadtxt(results_folder + '/reaktoro-seawater-pitzer-HCO3-p-100.0.txt')
    data_reaktoro = np.vstack((data_reaktoro_ca2.T, data_reaktoro_mg2.T, data_reaktoro_h.T, data_reaktoro_hco3.T))

    # Load from the file -d_CO2(g) concentrations
    data_phreeqc_ca2 = np.loadtxt(results_folder + '/phreeqc-pitzer-seawater-with-carbonates-p-100.txt', skiprows=2, usecols=(-8))
    data_phreeqc_mg2 = np.loadtxt(results_folder + '/phreeqc-pitzer-seawater-with-carbonates-p-100.txt', skiprows=2, usecols=(-7))
    data_phreeqc_h = np.loadtxt(results_folder + '/phreeqc-pitzer-seawater-with-carbonates-p-100.txt', skiprows=2, usecols=(-6))
    data_phreeqc_hco3 = np.loadtxt(results_folder + '/phreeqc-pitzer-seawater-with-carbonates-p-100.txt', skiprows=2, usecols=(-5))
    data_phreeqc = np.vstack((data_phreeqc_ca2, data_phreeqc_mg2, data_phreeqc_h, data_phreeqc_hco3))

    title = r'Aqueous species concentrations using pitzer.dat, P = 1 atm'
    plot_aqueous_concentrations(data_reaktoro, data_phreeqc, '-pitzer-dat-P-100-atm', title, results_folder)
