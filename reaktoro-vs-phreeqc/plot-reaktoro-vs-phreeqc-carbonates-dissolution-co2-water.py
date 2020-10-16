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

    colors = ['C1', 'C2', 'C3']
    plt.axes(xlim=(temperatures[0] - 2, temperatures[-1] + 2))

    plt.xlabel(r'Temperature [$\degree$C]', fontproperties=prop)
    plt.ylabel('Solubility [mol/kgw]', fontproperties=prop)
    plt.title(title, fontproperties=prop)

    plt.plot(temperatures, data_reaktoro[0], label='CO2(g), Reaktoro', **line(colors[0]))
    plt.plot(temperatures, data_phreeqc[0], 'D', **line_filled_marker(colors[0]))[0],
    plt.plot([], [], 'D', label='CO2(g), PHREEQC', **line_filled_marker(colors[0]))
    """
    plt.plot(temperatures, data_reaktoro[1], label='Calcite, Reaktoro', **line(colors[1]))
    plt.plot(temperatures, data_phreeqc[1], 's', **line_filled_marker(colors[1]))[0],
    plt.plot([], [], 's', label='Calcite, PHREEQC', **line_filled_marker(colors[1]))

    plt.plot(temperatures, data_reaktoro[2], label='Dolomite, Reaktoro', **line(colors[2]))
    plt.plot(temperatures, data_phreeqc[2], 'o', **line_filled_marker(colors[2]))[0],
    plt.plot([], [], 'o', label='Dolomite, PHREEQC', **line_filled_marker(colors[2]))
    """
    plt.legend(loc='upper right', prop=prop)
    plt.savefig(folder + '/reaktoro-phreeqc-comparison-co2' + tag +'.png')
    plt.close()

def plot_concentrations_minerals(data_reaktoro, data_phreeqc, tag, title, folder):

    colors = ['C2', 'C3']
    plt.axes(xlim=(temperatures[0] - 2, temperatures[-1] + 2))

    plt.xlabel(r'Temperature [$\degree$C]', fontproperties=prop)
    plt.ylabel('Solubility [mol/kgw]', fontproperties=prop)
    plt.title(title, fontproperties=prop)

    plt.plot(temperatures, data_reaktoro[1], label='Calcite, Reaktoro', **line(colors[0]))
    plt.plot(temperatures, data_phreeqc[1], 's', **line_filled_marker(colors[0]))[0],
    plt.plot([], [], 's', label='Calcite, PHREEQC', **line_filled_marker(colors[0]))

    plt.plot(temperatures, data_reaktoro[2], label='Dolomite, Reaktoro', **line(colors[1]))
    plt.plot(temperatures, data_phreeqc[2], 'o', **line_filled_marker(colors[1]))[0],
    plt.plot([], [], 'o', label='Dolomite, PHREEQC', **line_filled_marker(colors[1]))

    plt.legend(loc='upper right', prop=prop)
    plt.savefig(folder + '/reaktoro-phreeqc-comparison-calcite-dolomite' + tag +'.png')
    plt.close()

if __name__ == '__main__':

    results_folder = 'results-carbonates-dissolution-co2-saturated-water'

    # ---------------------------------------------------------------------------------------------------------------- #
    # P = 1 atm
    # ---------------------------------------------------------------------------------------------------------------- #

    # Load the Reaktoro values from the file
    data_reaktoro_co2g = np.loadtxt(results_folder + '/reaktoro-phreeqc-delta-CO2g-p-1.0.txt')
    data_reaktoro_calcite = np.loadtxt(results_folder + '/reaktoro-phreeqc-delta-Calcite-p-1.0.txt')
    data_reaktoro_dolomite = np.loadtxt(results_folder + '/reaktoro-phreeqc-delta-Dolomite-p-1.0.txt')
    data_reaktoro = np.vstack((data_reaktoro_co2g.T, data_reaktoro_calcite.T, data_reaktoro_dolomite.T))

    # Load the PHREEQC values from the file, -d_CO2(g) concentrations
    data_phreeqc_co2g = -np.loadtxt(results_folder + '/phreeqc-calcite-dolomite-dissolution-p-1.txt', skiprows=2, usecols=(-5))
    data_phreeqc_calcite = -np.loadtxt(results_folder + '/phreeqc-calcite-dolomite-dissolution-p-1.txt', skiprows=2, usecols=(-3))
    data_phreeqc_dolomite = -np.loadtxt(results_folder + '/phreeqc-calcite-dolomite-dissolution-p-1.txt', skiprows=2, usecols=(-1))
    data_phreeqc = np.vstack((data_phreeqc_co2g, data_phreeqc_calcite, data_phreeqc_dolomite))

    title = r'CO$_2$(g) solubility, P = 1 atm'
    plot_concentrations(data_reaktoro, data_phreeqc, '-P-1-atm', title, results_folder)
    title = r'Carbonates solubility, P = 1 atm'
    plot_concentrations_minerals(data_reaktoro, data_phreeqc, '-P-1-atm', title, results_folder)

    # ---------------------------------------------------------------------------------------------------------------- #
    # P = 100 atm
    # ---------------------------------------------------------------------------------------------------------------- #

    # Load the Reaktoro values from the file
    data_reaktoro_co2g = np.loadtxt(results_folder + '/reaktoro-phreeqc-delta-CO2g-p-100.0.txt')
    data_reaktoro_calcite = np.loadtxt(results_folder + '/reaktoro-phreeqc-delta-Calcite-p-100.0.txt')
    data_reaktoro_dolomite = np.loadtxt(results_folder + '/reaktoro-phreeqc-delta-Dolomite-p-100.0.txt')
    data_reaktoro = np.vstack((data_reaktoro_co2g.T, data_reaktoro_calcite.T, data_reaktoro_dolomite.T))

    # Load the PHREEQC values from the file, -d_CO2(g) concentrations
    data_phreeqc_co2g = -np.loadtxt(results_folder + '/phreeqc-calcite-dolomite-dissolution-p-100.txt', skiprows=2, usecols=(-5))
    data_phreeqc_calcite = -np.loadtxt(results_folder + '/phreeqc-calcite-dolomite-dissolution-p-100.txt', skiprows=2, usecols=(-3))
    data_phreeqc_dolomite = -np.loadtxt(results_folder + '/phreeqc-calcite-dolomite-dissolution-p-100.txt', skiprows=2, usecols=(-1))
    data_phreeqc = np.vstack((data_phreeqc_co2g, data_phreeqc_calcite, data_phreeqc_dolomite))

    title = r'CO$_2$(g) solubility, P = 100 atm'
    plot_concentrations(data_reaktoro, data_phreeqc, '-P-100-atm', title, results_folder)
    title = r'Carbonates solubility, P = 100 atm'
    plot_concentrations_minerals(data_reaktoro, data_phreeqc, '-P-100-atm', title, results_folder)

