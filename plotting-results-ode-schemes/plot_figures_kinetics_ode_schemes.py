import numpy as np
import matplotlib.pyplot as plt
import os
from natsort import natsorted
import matplotlib as mpl
import json

# Options for the figure plotting
#plot_at_selected_steps = [1, 10, 60, 120, 480, 960, 1000, 1200, 2400, 3600, 5800, 7200]  # the time steps at which the results are plotted
plot_at_selected_steps = [1, 10, 60, 120, 240, 480, 960, 1000]  # the time steps at which the results are plotted
#plot_at_selected_steps = [1, 10, 20, 40, 80, 100]  # the time steps at which the results are plotted

# Auxiliary time related constants
second = 1
minute = 60
hour = 60 * minute
day = 24 * hour
year = 365 * day

# Discretisation parameters
xl = 0.0          # the x-coordinate of the left boundary
xr = 1.0          # the x-coordinate of the right boundary
nsteps = 1000    # the number of steps in the reactive transport simulation
ncells = 100      # the number of cells in the discretization
eqreltol = 1e-1   # relative tolerance in equilibrium
eqabstol = 1e-12  # absolute tolerance in equilibrium
kinreltol = 1e-1  # relative tolerance in kinetics
kinabstol = 1e-5  # absolute tolerance in kinetics
cutoff = -1e-5    # absolute tolerance in kinetics

D  = 1.0e-9       # the diffusion coefficient (in units of m2/s)
v  = 1.0/day      # the fluid pore velocity (in units of m/s)
dt = 30 * minute  # the time step (in units of s)
T = 60.0          # the temperature (in units of K)
P = 100           # the pressure (in units of Pa)
phi = 0.1           # the porosity

dirichlet = False # the parameter that defines whether Dirichlet BC must be used
smrt_solv = True  # the parameter that defines whether classic or smart
                  # EquilibriumSolver must be used

tag_smart = "-dt-" + "{:d}".format(dt) + \
      "-ncells-" + str(ncells) + \
      "-nsteps-" + str(nsteps) + \
      "-eqrel-" + "{:.{}e}".format(eqreltol, 1) + \
      "-eqabs-" + "{:.{}e}".format(eqabstol, 1) + \
      "-kinrel-" + "{:.{}e}".format(kinreltol, 1) + \
      "-kinabs-" + "{:.{}e}".format(kinabstol, 1) + \
      "-cutoff-" + "{:.{}e}".format(cutoff, 1)

tag_class = "-dt-" + "{:d}".format(dt) + \
      "-ncells-" + str(ncells) + \
      "-nsteps-" + str(nsteps)

test_tag_smart = tag_smart + "-smart-kin-conv-eq"
test_tag_class = tag_class + "-conv-kin-conv-eq"

folder = 'rt-sa-50'
folder_class = 'rt-sa-50'
folder_smart   = "../" + folder + test_tag_smart
folder_class   = "../" + folder_class + test_tag_class
folder_general = "results-smart-kin-" + folder + tag_smart

os.system('mkdir -p ' + folder_general)

fillz = len(str(123))

# Indices of the loaded data to plot
indx_ph        = 0
indx_Hcation   = 1
indx_Cacation  = 2
indx_Mgcation  = 3
indx_HCO3anion = 4
indx_CO2aq     = 5
indx_calcite   = 6
indx_dolomite  = 7
indx_mol_calcite   = 8
indx_mol_dolomite  = 9

# Plotting params
circ_area = 6 ** 2
custom_font = { }
time_steps = np.linspace(0, nsteps, nsteps)

C0 = '#107ab0'
C1 = '#fc5a50'

xcells = np.linspace(xl, xr, ncells)  # the x-coordinates of the plots

#font
#mpl.rcParams['font.family'] = 'sans-serif'
#mpl.rcParams['font.sans-serif'] = 'Fira Sans'
# mpl.rcParams['font.family'] = 'sans-serif'
# mpl.rcParams['font.serif'] ='TeXGyreSchola'
# mpl.rcParams['font.sans-serif'] ='Fira Code'
# mpl.rcParams['font.sans-serif'] = 'Century Gothic'
# mpl.rcParams['font.size'] = 12
# mpl.rcParams['legend.fontsize'] = 'medium'
#tick label spacing and tick width
# mpl.rcParams['xtick.major.pad'] = 4
# mpl.rcParams['ytick.major.pad'] = 5
# mpl.rcParams['xtick.major.width'] = 1
# mpl.rcParams['ytick.major.width'] = 1
#legend style
# mpl.rcParams['legend.frameon'] = True
# mpl.rcParams['legend.numpoints'] = 3
#mpl.rcParams['backend'] = 'PDF'
#mpl.rcParams['savefig.dpi'] = 200

def empty_marker(color):
    return {'facecolor': 'white', 'edgecolor': color, 's': circ_area, 'zorder': 2, 'linewidths': 1.5 }

def filled_marker(color):
    return {'color': color, 's': circ_area, 'zorder': 2, 'linewidths': 1.5 }

def line_empty_marker(color):
    return {'markerfacecolor': 'white', 'markeredgecolor':color, 'markersize': 6, 'markeredgewidth': 1.5 }

def line_filled_marker(color):
    return {'color': color, 'markersize': 6, 'markeredgewidth': 1.5 }

def line(color):
    return {'linestyle': '-', 'color': color, 'zorder': 1, 'linewidth': 2}

def titlestr(t):
    d = int(t / day)                 # The number of days
    h = int(int(t % day) / hour)     # The number of remaining hours
    m = int(int(t % hour) / minute)  # The number of remaining minutes
    return '{:>3d}d {:>2}h {:>2}m'.format(int(d), str(int(h)).zfill(2), str(int(m)).zfill(2))


def plot_figures_ph():

    for i in plot_at_selected_steps:
        print("On pH figure at time step: {}".format(i))
        t = i * dt
        filearray_class = np.loadtxt(folder_class + '/' + files_class[i-1], skiprows=1)
        filearray_smart = np.loadtxt(folder_smart + '/' + files_smart[i-1], skiprows=1)
        data_class = filearray_class.T
        data_smart = filearray_smart.T
        data_class_ph = data_class[indx_ph]
        data_smart_ph = data_smart[indx_ph]
        plt.axes(xlim=(-0.01, 0.501), ylim=(2.5, 12.0))
        plt.xlabel('Distance [m]')
        plt.ylabel('pH')
        plt.title(titlestr(t))
        plt.plot(xcells, data_class_ph, label='pH', **line('teal'))
        plt.plot(xcells[status[i-1]==0], data_smart_ph[status[i-1]==0], 'o', **line_empty_marker('teal'))
        plt.plot(xcells[status[i-1]==1], data_smart_ph[status[i-1]==1], 'o', **line_filled_marker('teal'))
        plt.plot([], [], 'o', label='Smart Prediction', **line_filled_marker('black'))
        plt.plot([], [], 'o', label='Learning', **line_empty_marker('black'))
        plt.legend(loc='lower right')
        plt.savefig(folder_general + '/pH-{}.png'.format(i))
        plt.tight_layout()
        plt.close()


def plot_figures_calcite_dolomite():

    for i in plot_at_selected_steps:
        print("On calcite-dolomite (mineral volume) figure at time step: {}".format(i))
        t = i * dt
        filearray_class = np.loadtxt(folder_class + '/' + files_class[i-1], skiprows=1)
        filearray_smart = np.loadtxt(folder_smart + '/' + files_smart[i-1], skiprows=1)
        data_class = filearray_class.T
        data_smart = filearray_smart.T
        data_class_calcite, data_class_dolomite = data_class[indx_calcite], data_class[indx_dolomite]
        data_smart_calcite, data_smart_dolomite = data_smart[indx_calcite], data_smart[indx_dolomite]
        plt.axes(xlim=(-0.01, 0.501), ylim=(-0.1, 2.1))
        plt.xlabel('Distance [m]')
        plt.ylabel('Mineral Volume [%$_{\mathsf{vol}}$]')
        plt.title(titlestr(t))
        plt.plot(xcells, data_class_calcite * 100/(1 - phi), label='Calcite', **line('C0'))
        plt.plot(xcells, data_class_dolomite * 100/(1 - phi), label='Dolomite', **line('C1'))
        plt.plot(xcells[status[i-1]==0], data_smart_calcite[status[i-1]==0] * 100/(1 - phi), 'o', **line_empty_marker('C0'))
        plt.plot(xcells[status[i-1]==1], data_smart_calcite[status[i-1]==1] * 100/(1 - phi), 'o', **line_filled_marker('C0'))
        plt.plot(xcells[status[i-1]==0], data_smart_dolomite[status[i-1]==0] * 100/(1 - phi), 'o', **line_empty_marker('C1'))
        plt.plot(xcells[status[i-1]==1], data_smart_dolomite[status[i-1]==1] * 100/(1 - phi), 'o', **line_filled_marker('C1'))
        plt.plot([], [], 'o', label='Smart Prediction', **line_filled_marker('black'))
        plt.plot([], [], 'o', label='Learning', **line_empty_marker('black'))
        plt.legend(loc='center right')
        plt.savefig(folder_general + '/calcite-dolomite-{}.png'.format(i))
        plt.tight_layout()
        plt.close()

def plot_figures_calcite_dolomite_moles():

    for i in plot_at_selected_steps:
        print("On calcite-dolomite figure at time step: {}".format(i))
        t = i * dt
        filearray_class = np.loadtxt(folder_class + '/' + files_class[i-1], skiprows=1)
        filearray_smart = np.loadtxt(folder_smart + '/' + files_smart[i-1], skiprows=1)
        data_class = filearray_class.T
        data_smart = filearray_smart.T
        data_class_calcite, data_class_dolomite = data_class[indx_mol_calcite], data_class[indx_mol_dolomite]
        data_smart_calcite, data_smart_dolomite = data_smart[indx_mol_calcite], data_smart[indx_mol_dolomite]
        plt.axes(xlim=(-0.01, 0.501))
        plt.ylabel('Concentration [mol/m3]')
        plt.xlabel('Distance [m]')
        plt.title(titlestr(t))
        plt.plot(xcells, data_class_calcite, label='Calcite', **line('C0'))
        plt.plot(xcells, data_class_dolomite, label='Dolomite', **line('C1'))
        plt.plot(xcells[status[i-1]==0], data_smart_calcite[status[i-1]==0], 'o', **line_empty_marker('C0'))
        plt.plot(xcells[status[i-1]==1], data_smart_calcite[status[i-1]==1], 'o', **line_filled_marker('C0'))
        plt.plot(xcells[status[i-1]==0], data_smart_dolomite[status[i-1]==0], 'o', **line_empty_marker('C1'))
        plt.plot(xcells[status[i-1]==1], data_smart_dolomite[status[i-1]==1], 'o', **line_filled_marker('C1'))
        plt.plot([], [], 'o', label='Smart Prediction', **line_filled_marker('black'))
        plt.plot([], [], 'o', label='Learning', **line_empty_marker('black'))
        plt.legend(loc='center right')
        plt.savefig(folder_general + '/calcite-dolomite-moles-{}.png'.format(i))
        plt.tight_layout()
        plt.close()


def plot_figures_aqueous_species():

    for i in plot_at_selected_steps:
        print("On aqueous-species figure at time step: {}".format(i))
        t = i * dt
        filearray_class = np.loadtxt(folder_class + '/' + files_class[i-1], skiprows=1)
        filearray_smart = np.loadtxt(folder_smart + '/' + files_smart[i-1], skiprows=1)
        data_class = filearray_class.T
        data_smart = filearray_smart.T
        data_class_cacation  = data_class[indx_Cacation]
        data_class_mgcation  = data_class[indx_Mgcation]
        data_class_hco3anion = data_class[indx_HCO3anion]
        data_class_co2aq     = data_class[indx_CO2aq]
        data_class_hcation   = data_class[indx_Hcation]
        data_smart_cacation  = data_smart[indx_Cacation]
        data_smart_mgcation  = data_smart[indx_Mgcation]
        data_smart_hco3anion = data_smart[indx_HCO3anion]
        data_smart_co2aq     = data_smart[indx_CO2aq]
        data_smart_hcation   = data_smart[indx_Hcation]
        plt.axes(xlim=(-0.01, 0.501), ylim=(0.5e-5, 2))
        plt.xlabel('Distance [m]')
        plt.ylabel('Concentration [molal]')
        plt.yscale('log')
        plt.title(titlestr(t))
        plt.plot(xcells, data_class_cacation, label=r'$\mathrm{Ca^{2+}}$', **line('C0'))[0],
        plt.plot(xcells, data_class_mgcation, label=r'$\mathrm{Mg^{2+}}$', **line('C1'))[0],
        plt.plot(xcells, data_class_hco3anion, label=r'$\mathrm{HCO_3^{-}}$',**line('C2'))[0],
        plt.plot(xcells, data_class_co2aq, label=r'$\mathrm{CO_2(aq)}$',**line('red'))[0],
        plt.plot(xcells, data_class_hcation, label=r'$\mathrm{H^+}$', **line('darkviolet'))[0],
        plt.plot(xcells[status[i-1]==0], data_smart_cacation[status[i-1]==0], 'o', **line_empty_marker('C0'))[0],
        plt.plot(xcells[status[i-1]==1], data_smart_cacation[status[i-1]==1], 'o', **line_filled_marker('C0'))[0],
        plt.plot(xcells[status[i-1]==0], data_smart_mgcation[status[i-1]==0], 'o', **line_empty_marker('C1'))[0],
        plt.plot(xcells[status[i-1]==1], data_smart_mgcation[status[i-1]==1], 'o', **line_filled_marker('C1'))[0],
        plt.plot(xcells[status[i-1]==0], data_smart_hco3anion[status[i-1]==0], 'o', **line_empty_marker('C2'))[0],
        plt.plot(xcells[status[i-1]==1], data_smart_hco3anion[status[i-1]==1], 'o', **line_filled_marker('C2'))[0],
        plt.plot(xcells[status[i-1]==0], data_smart_co2aq[status[i-1]==0], 'o', **line_empty_marker('red'))[0],
        plt.plot(xcells[status[i-1]==1], data_smart_co2aq[status[i-1]==1], 'o', **line_filled_marker('red'))[0],
        plt.plot(xcells[status[i-1]==0], data_smart_hcation[status[i-1]==0], 'o', **line_empty_marker('darkviolet'))[0],
        plt.plot(xcells[status[i-1]==1], data_smart_hcation[status[i-1]==1], 'o', **line_filled_marker('darkviolet'))[0],
        plt.plot([], [], 'o', label='Smart Prediction', **line_filled_marker('black'))
        plt.plot([], [], 'o', label='Learning', **line_empty_marker('black'))
        plt.legend(loc='upper right')
        plt.savefig(folder_general + '/aqueous-species-{}.png'.format(i))
        plt.tight_layout()
        plt.close()


def plot_computing_costs():


    with open(folder_class + '/analysis-conventional-kin-conventional-eq.json') as read_file:
        data_class = json.load(read_file)
    with open(folder_smart + '/analysis-smart-kin-conventional-eq.json') as read_file:
        data_smart = json.load(read_file)

    timing_class = data_class.get('computing_costs_per_time_step')
    timing_smart = data_smart.get('computing_costs_per_time_step')

    timings_transport = np.array(timing_class.get('transport')) * 1e6  # in microseconds
    timings_kinetics_class = np.array(timing_class.get('kinetics')) * 1e6  # in microseconds
    timings_kinetics_smart = np.array(timing_smart.get('smart_kinetics')) * 1e6  # in microseconds
    timings_kinetics_smart_ideal = (np.array(timing_smart.get('smart_kinetics')) \
                                    - np.array(timing_smart.get('smart_kinetics_nearest_neighbor_search')) \
                                    - np.array(timing_smart.get('smart_kinetics_chemical_properties'))) * 1e6  # in microseconds

    plt.xlabel('Time Step')
    plt.ylabel('Computing Cost [μs]')
    plt.yscale('log')
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(timings_kinetics_class, label="Chemical Kinetics (Conventional)", color='C0', linewidth=2)
    plt.plot(timings_kinetics_smart, label="Chemical Kinetics (Smart)", color='C1', linewidth=2)
    #plt.plot(timings_kinetics_smart_ideal, label="Chemical Kinetics (Ideal Smart)", color='C3', linewidth=2)
    plt.plot(timings_transport, label="Transport", color='C2', linewidth=2)
    #leg = plt.legend(loc='lower right', bbox_to_anchor=(1, 0.13))
    #for line in leg.get_lines(): line.set_linewidth(2.0)
    #plt.tight_layout()
    plt.savefig(folder_general + '/computing-costs-nolegend.png')
    #plt.savefig(folder_general + '/computing-costs-nolegend-with-ideal.png')
    plt.close()


def plot_on_demand_learning_countings():

    # Load the status data, where 0 stands for conventional learning and 1 for smart prediction
    with open(folder_smart + '/analysis-smart-kin-conventional-eq.json') as read_file:
        data = json.load(read_file)

    status_learnings = data.get('smart_kinetics_cells_where_learning_was_required_at_step')

    # Collect the number of learnings on each step
    learnings = [len(x) for x in status_learnings]

    plt.xlabel('Time Step')
    plt.ylim(bottom=0, top=np.max(learnings)+1)
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(learnings, color='C0', linewidth=2)
    plt.tight_layout()
    plt.savefig(folder_general + '/on-demand-learning-countings.png')
    plt.close()


def plot_on_demand_learnings_total():

    step = 1

    # Load the status data, where 0 stands for conventional learning and 1 for smart prediction
    with open(folder_smart + '/analysis-smart-kin-conventional-eq.json') as read_file:
        data = json.load(read_file)

    status_learnings = data.get('smart_kinetics_cells_where_learning_was_required_at_step')

    # Collect the number of learnings on each step and accumulate them w.r.t. time-step
    learnings = [len(x) for x in status_learnings]
    accum_learnings = [np.sum(learnings[0:i]) for i in range(0, nsteps)]

    plt.xlabel('Time Step')
    plt.ylim(bottom=0, top=np.max(accum_learnings)+1)
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(time_steps[0:nsteps:step], accum_learnings[0:nsteps:step], color='C0', linewidth=2)
    plt.tight_layout()
    plt.savefig(folder_general + '/on-demand-learning-total.png')
    plt.close()


def plot_computing_costs_of_estimation():

    step = 20
    with open(folder_smart + '/analysis-smart-kin-conventional-eq.json') as read_file:
        data_smart = json.load(read_file)

    timing_smart = data_smart.get('computing_costs_per_time_step')

    timings_estimate = np.array(timing_smart.get('smart_kinetics_estimate')) * 1e6       # in microseconds
    timings_search = np.array(timing_smart.get('smart_kinetics_nearest_neighbor_search')) * 1e6  # in microseconds

    plt.xlabel('Time Step')
    plt.ylabel('Computing Cost [μs]')
    #plt.xscale('log')
    plt.yscale('log')
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(time_steps[0:nsteps:step], timings_search[0:nsteps:step], label="Search", color='C0', linewidth=2)
    plt.plot(time_steps[0:nsteps:step], timings_estimate[0:nsteps:step], label="Estimation", color='C1', linewidth=2)

    leg = plt.legend(loc='lower right')
    for line in leg.get_lines(): line.set_linewidth(2.0)
    plt.tight_layout()
    plt.savefig(folder_general + '/estimate-search-costs.png')
    plt.close()


def plot_speedups():

    step = 80
    # Load the status data, where 0 stands for conventional learning and 1 for smart prediction
    with open(folder_smart + '/analysis-smart-kin-conventional-eq.json') as read_file:
        data_smart = json.load(read_file)
    with open(folder_class + '/analysis-conventional-kin-conventional-eq.json') as read_file:
        data_class = json.load(read_file)

    timing_class = data_class.get('computing_costs_per_time_step')
    timing_smart = data_smart.get('computing_costs_per_time_step')

    speedup = np.array(timing_class.get('kinetics')) / np.array(timing_smart.get('smart_kinetics'))
    speedup_ideal = np.array(timing_class.get('kinetics')) / \
                    (np.array(timing_smart.get('smart_kinetics'))
                     - np.array(timing_smart.get('smart_kinetics_nearest_neighbor_search'))
                     - np.array(timing_smart.get('smart_kinetics_chemical_properties')))

    plt.xlabel('Time Step')
    plt.ylabel('Speedup (-)')
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(time_steps[0:nsteps:step], speedup_ideal[0:nsteps:step], label="Conventional vs. Smart (Ideal search)", color='C2', linewidth=2)
    plt.plot(time_steps[0:nsteps:step], speedup[0:nsteps:step], label="Conventional vs. Smart ", color='C0', linewidth=2)

    #leg = plt.legend(loc='upper left')
    #for line in leg.get_lines():
    #    line.set_linewidth(2.0)
    #plt.tight_layout()
    #plt.savefig(folder_general + '/speedups.png')
    #plt.show()
    plt.savefig(folder_general + '/speedups-nolegend.png')
    plt.close()

def count_trainings(status):
    counter = [st for st in status if st == 0]
    return (len(counter), len(status))

if __name__ == '__main__':

    # Load the status data, where 0 stands for conventional learning and 1 for smart prediction
    with open(folder_smart + '/analysis-smart-kin-conventional-eq.json') as read_file:
        status_data = json.load(read_file)

    status_learnings = status_data.get('smart_kinetics_cells_where_learning_was_required_at_step')
    status = np.ones([nsteps, ncells])
    for i in range(0, nsteps):
        for j in range(0, len(status_learnings[i])):
            status[i][j] = 0

    # Count the percentage of the trainings needed
    training_counter = count_trainings(np.array(status).flatten())
    title = "%2.2f percent is training (%d out of %d cells)" % (
        100 * training_counter[0] / training_counter[1], training_counter[0], training_counter[1])
    text = ["Number of chemical equilibrium trainings: %d" % training_counter[0],
            "Number of smart chemical equilibrium estimations: %d" % (training_counter[1] - training_counter[0]),
            ("Percentage of smart chemical kineics predictions: %2.2f" % (
                    100 - 100 * training_counter[0] / training_counter[1])) + "$\%$"]
    print(title)

    print("Collecting files...")
    # Collect files with results corresponding to smart or reference (classical) solver
    files_smart = [file for file in natsorted( os.listdir(folder_smart) ) if ("analysis" not in file)]
    files_class = [file for file in natsorted( os.listdir(folder_class) ) if ("analysis" not in file)]

    plot_computing_costs()
    plot_on_demand_learning_countings()
    plot_on_demand_learnings_total()
    plot_speedups()
    plot_computing_costs_of_estimation()
    plot_figures_ph()
    plot_figures_calcite_dolomite()
    plot_figures_aqueous_species()
    plot_figures_calcite_dolomite_moles()