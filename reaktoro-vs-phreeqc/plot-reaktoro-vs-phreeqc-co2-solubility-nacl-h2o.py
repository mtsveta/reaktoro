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

def plot_concentrations(data_reaktoro, data_phreeqc, molalities, tag, title, folder):

    colors = ['C1', 'C2', 'C3']
    plt.axes(xlim=(temperatures[0] - 2, temperatures[-1] + 2))

    for molality, i in zip(molalities, list(range(len(temperatures)))):

        plt.xlabel(r'Temperature [°C]', fontproperties=prop)
        plt.ylabel('Solubility [mol/kgw]', fontproperties=prop)
        plt.title(title, fontproperties=prop)

        plt.plot(temperatures, data_reaktoro[i], label='Reaktoro, ' + str(molality) + ' molal', **line(colors[i]))
        plt.plot(temperatures, data_phreeqc[i], 'D', **line_filled_marker(colors[i]))[0],
        plt.plot([], [], 'D', label='PHREEQC, ' + str(molality) + ' molal', **line_filled_marker(colors[i]))

    plt.legend(loc='upper right', prop=prop)
    plt.savefig(folder + '/reaktoro-phreeqc-comparison' + tag +'.pdf')
    plt.close()

if __name__ == '__main__':

    molalities = [1.0, 2.0, 4.0]

    results_folder = 'results-co2-solubility-nacl-h2o'

    # Load the Reaktoro values from the file
    data_reaktoro_nacl_1 = np.loadtxt(results_folder + '/reaktoro-phreeqc-nacl-1-p-1.0-atm.txt')
    data_reaktoro_nacl_2 = np.loadtxt(results_folder + '/reaktoro-phreeqc-nacl-2-p-1.0-atm.txt')
    data_reaktoro_nacl_4 = np.loadtxt(results_folder + '/reaktoro-phreeqc-nacl-4-p-1.0-atm.txt')
    data_reaktoro = np.vstack((data_reaktoro_nacl_1.T,data_reaktoro_nacl_2.T, data_reaktoro_nacl_4.T))
    #print(data_reaktoro)
    #input()

    # Load the PHREEQC values from the file -d_CO2(g) concentrations
    data_phreeqc_nacl_1 = np.loadtxt(results_folder + '/phreeqc-nacl-1-p-1.txt', skiprows=2, usecols=(-3))
    data_phreeqc_nacl_2 = np.loadtxt(results_folder + '/phreeqc-nacl-2-p-1.txt', skiprows=2, usecols=(-3))
    data_phreeqc_nacl_4 = np.loadtxt(results_folder + '/phreeqc-nacl-4-p-1.txt', skiprows=2, usecols=(-3))
    data_phreeqc = np.vstack((data_phreeqc_nacl_1, data_phreeqc_nacl_2, data_phreeqc_nacl_4))
    #print(data_phreeqc)
    #input()

    title = r'Solubility of CO$_2$ in NaCl brine phreeqc.dat, P = 1 atm'
    plot_concentrations(data_reaktoro, data_phreeqc, molalities, '-nacl-P-1-atm', title, results_folder)

    # Load the Reaktoro values from the file
    data_reaktoro_nacl_1 = np.loadtxt(results_folder + '/reaktoro-pitzer-nacl-1-p-1.0-atm.txt')
    data_reaktoro_nacl_2 = np.loadtxt(results_folder + '/reaktoro-pitzer-nacl-2-p-1.0-atm.txt')
    data_reaktoro_nacl_4 = np.loadtxt(results_folder + '/reaktoro-pitzer-nacl-4-p-1.0-atm.txt')
    data_reaktoro_pitzer = np.vstack((data_reaktoro_nacl_1.T,data_reaktoro_nacl_2.T, data_reaktoro_nacl_4.T))

    # Load the PHREEQC values from the file -d_CO2(g) concentrations
    data_phreeqc_nacl_1 = np.loadtxt(results_folder + '/phreeqc-pitzer-nacl-1-p-1.txt', skiprows=2, usecols=(-3))
    data_phreeqc_nacl_2 = np.loadtxt(results_folder + '/phreeqc-pitzer-nacl-2-p-1.txt', skiprows=2, usecols=(-3))
    data_phreeqc_nacl_4 = np.loadtxt(results_folder + '/phreeqc-pitzer-nacl-4-p-1.txt', skiprows=2, usecols=(-3))
    data_phreeqc_pitzer = np.vstack((data_phreeqc_nacl_1, data_phreeqc_nacl_2, data_phreeqc_nacl_4))

    title = r'Solubility of CO$_2$ in NaCl brine using pitzer.dat, P = 1 atm'
    plot_concentrations(data_reaktoro_pitzer, data_phreeqc_pitzer, molalities, '-pitzer-nacl-P-1-atm', title, results_folder)

    # P = 100 atm
    # Load the Reaktoro values from the file
    data_reaktoro_nacl_1 = np.loadtxt(results_folder + '/reaktoro-phreeqc-nacl-1-p-100.0-atm.txt')
    data_reaktoro_nacl_2 = np.loadtxt(results_folder + '/reaktoro-phreeqc-nacl-2-p-100.0-atm.txt')
    data_reaktoro_nacl_4 = np.loadtxt(results_folder + '/reaktoro-phreeqc-nacl-4-p-100.0-atm.txt')
    data_reaktoro = np.vstack((data_reaktoro_nacl_1.T,data_reaktoro_nacl_2.T, data_reaktoro_nacl_4.T))

    # Load the PHREEQC values from the file -d_CO2(g) concentrations
    data_phreeqc_nacl_1 = np.loadtxt(results_folder + '/phreeqc-nacl-1-p-100.txt', skiprows=2, usecols=(-3))
    data_phreeqc_nacl_2 = np.loadtxt(results_folder + '/phreeqc-nacl-2-p-100.txt', skiprows=2, usecols=(-3))
    data_phreeqc_nacl_4 = np.loadtxt(results_folder + '/phreeqc-nacl-4-p-100.txt', skiprows=2, usecols=(-3))
    data_phreeqc = np.vstack((data_phreeqc_nacl_1, data_phreeqc_nacl_2, data_phreeqc_nacl_4))

    title = r'Solubility of CO$_2$ in NaCl brine using phreeqc.dat, P = 100 atm'
    plot_concentrations(data_reaktoro, data_phreeqc, molalities, '-nacl-P-100-atm', title, results_folder)

    # P = 100 atm
    # Load the Reaktoro values from the file
    data_reaktoro_nacl_1 = np.loadtxt(results_folder + '/reaktoro-pitzer-nacl-1-p-100.0-atm.txt')
    data_reaktoro_nacl_2 = np.loadtxt(results_folder + '/reaktoro-pitzer-nacl-2-p-100.0-atm.txt')
    data_reaktoro_nacl_4 = np.loadtxt(results_folder + '/reaktoro-pitzer-nacl-4-p-100.0-atm.txt')
    data_reaktoro_pitzer = np.vstack((data_reaktoro_nacl_1.T,data_reaktoro_nacl_2.T, data_reaktoro_nacl_4.T))

    # Load the PHREEQC from the file -d_CO2(g) concentrations
    data_phreeqc_nacl_1 = np.loadtxt(results_folder + '/phreeqc-pitzer-nacl-1-p-100.txt', skiprows=2, usecols=(-3))
    data_phreeqc_nacl_2 = np.loadtxt(results_folder + '/phreeqc-pitzer-nacl-2-p-100.txt', skiprows=2, usecols=(-3))
    data_phreeqc_nacl_4 = np.loadtxt(results_folder + '/phreeqc-pitzer-nacl-4-p-100.txt', skiprows=2, usecols=(-3))
    data_phreeqc_pitzer = np.vstack((data_phreeqc_nacl_1, data_phreeqc_nacl_2, data_phreeqc_nacl_4))

    title = r'Solubility of CO$_2$ in NaCl brine using pitzer.dat, P = 100 atm'
    plot_concentrations(data_reaktoro_pitzer, data_phreeqc_pitzer, molalities, '-pitzer-nacl-P-100-atm', title, results_folder)



    # # Load the Reaktoro values from the file
    # data_reaktoro_nacl_1 = np.loadtxt(results_folder + '/reaktoro-phreeqc-nacl-1-p-1.0-atm.txt')
    # data_reaktoro_nacl_2 = np.loadtxt(results_folder + '/reaktoro-phreeqc-nacl-2-p-1.0-atm.txt')
    # data_reaktoro_nacl_4 = np.loadtxt(results_folder + '/reaktoro-phreeqc-nacl-4-p-1.0-atm.txt')
    # data_reaktoro = np.vstack((data_reaktoro_nacl_1.T,data_reaktoro_nacl_2.T, data_reaktoro_nacl_4.T))
    # print(data_reaktoro)
    # input()
    #
    # # Load the PHREEQC values from the file -d_CO2(g) concentrations
    # data_phreeqc_nacl_1 = -np.loadtxt(results_folder + '/phreeqc-nacl-1-p-1.txt', skiprows=2, usecols=(-1))
    # data_phreeqc_nacl_2 = -np.loadtxt(results_folder + '/phreeqc-nacl-2-p-1.txt', skiprows=2, usecols=(-1))
    # data_phreeqc_nacl_4 = -np.loadtxt(results_folder + '/phreeqc-nacl-4-p-1.txt', skiprows=2, usecols=(-1))
    # data_phreeqc = np.vstack((data_phreeqc_nacl_1, data_phreeqc_nacl_2, data_phreeqc_nacl_4))
    # print(data_phreeqc)
    # input()
    #
    # title = r'Solubility of CO$_2$ in NaCl brine, P = 1 atm'
    # plot_concentrations(data_reaktoro, data_phreeqc, molalities, '-nacl-P-1-atm', title, results_folder)
    #
    # # Load the Reaktoro values from the file
    # data_reaktoro_nacl_1 = np.loadtxt(results_folder + '/reaktoro-pitzer-nacl-1-p-1.0-atm.txt')
    # data_reaktoro_nacl_2 = np.loadtxt(results_folder + '/reaktoro-pitzer-nacl-2-p-1.0-atm.txt')
    # data_reaktoro_nacl_4 = np.loadtxt(results_folder + '/reaktoro-pitzer-nacl-4-p-1.0-atm.txt')
    # data_reaktoro_pitzer = np.vstack((data_reaktoro_nacl_1.T,data_reaktoro_nacl_2.T, data_reaktoro_nacl_4.T))
    #
    # # Load the PHREEQC values from the file -d_CO2(g) concentrations
    # data_phreeqc_nacl_1 = -np.loadtxt(results_folder + '/phreeqc-pitzer-nacl-1-p-1.txt', skiprows=2, usecols=(-1))
    # data_phreeqc_nacl_2 = -np.loadtxt(results_folder + '/phreeqc-pitzer-nacl-2-p-1.txt', skiprows=2, usecols=(-1))
    # data_phreeqc_nacl_4 = -np.loadtxt(results_folder + '/phreeqc-pitzer-nacl-4-p-1.txt', skiprows=2, usecols=(-1))
    # data_phreeqc_pitzer = np.vstack((data_phreeqc_nacl_1, data_phreeqc_nacl_2, data_phreeqc_nacl_4))
    #
    # title = r'Solubility of CO$_2$ in NaCl brine (Pitzer model), P = 1 atm'
    # plot_concentrations(data_reaktoro_pitzer, data_phreeqc_pitzer, molalities, '-pitzer-nacl-P-1-atm', title, results_folder)
    #
    # # P = 100 atm
    # # Load the Reaktoro values from the file
    # data_reaktoro_nacl_1 = np.loadtxt(results_folder + '/reaktoro-phreeqc-nacl-1-p-100.0-atm.txt')
    # data_reaktoro_nacl_2 = np.loadtxt(results_folder + '/reaktoro-phreeqc-nacl-2-p-100.0-atm.txt')
    # data_reaktoro_nacl_4 = np.loadtxt(results_folder + '/reaktoro-phreeqc-nacl-4-p-100.0-atm.txt')
    # data_reaktoro = np.vstack((data_reaktoro_nacl_1.T,data_reaktoro_nacl_2.T, data_reaktoro_nacl_4.T))
    #
    # # Load the PHREEQC values from the file -d_CO2(g) concentrations
    # data_phreeqc_nacl_1 = -np.loadtxt(results_folder + '/phreeqc-nacl-1-p-100.txt', skiprows=2, usecols=(-1))
    # data_phreeqc_nacl_2 = -np.loadtxt(results_folder + '/phreeqc-nacl-2-p-100.txt', skiprows=2, usecols=(-1))
    # data_phreeqc_nacl_4 = -np.loadtxt(results_folder + '/phreeqc-nacl-4-p-100.txt', skiprows=2, usecols=(-1))
    # data_phreeqc = np.vstack((data_phreeqc_nacl_1, data_phreeqc_nacl_2, data_phreeqc_nacl_4))
    #
    # title = r'Solubility of CO$_2$ in NaCl brine, P = 100 atm'
    # plot_concentrations(data_reaktoro, data_phreeqc, molalities, '-nacl-P-100-atm', title, results_folder)
    #
    # # P = 100 atm
    # # Load the Reaktoro values from the file
    # data_reaktoro_nacl_1 = np.loadtxt(results_folder + '/reaktoro-pitzer-nacl-1-p-100.0-atm.txt')
    # data_reaktoro_nacl_2 = np.loadtxt(results_folder + '/reaktoro-pitzer-nacl-2-p-100.0-atm.txt')
    # data_reaktoro_nacl_4 = np.loadtxt(results_folder + '/reaktoro-pitzer-nacl-4-p-100.0-atm.txt')
    # data_reaktoro_pitzer = np.vstack((data_reaktoro_nacl_1.T,data_reaktoro_nacl_2.T, data_reaktoro_nacl_4.T))
    #
    # # Load the PHREEQC from the file -d_CO2(g) concentrations
    # data_phreeqc_nacl_1 = -np.loadtxt(results_folder + '/phreeqc-pitzer-nacl-1-p-100.txt', skiprows=2, usecols=(-1))
    # data_phreeqc_nacl_2 = -np.loadtxt(results_folder + '/phreeqc-pitzer-nacl-2-p-100.txt', skiprows=2, usecols=(-1))
    # data_phreeqc_nacl_4 = -np.loadtxt(results_folder + '/phreeqc-pitzer-nacl-4-p-100.txt', skiprows=2, usecols=(-1))
    # data_phreeqc_pitzer = np.vstack((data_phreeqc_nacl_1, data_phreeqc_nacl_2, data_phreeqc_nacl_4))
    #
    # title = r'Solubility of CO$_2$ in NaCl brine (Pitzer model), P = 100 atm'
    # plot_concentrations(data_reaktoro_pitzer, data_phreeqc_pitzer, molalities, '-pitzer-nacl-P-100-atm', title, results_folder)
