import numpy  as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import colors as mcolors

minute = 60
hour = 60 * minute
day = 24 * hour

# https://github.com/sciplot/gnuplot-palettes/blob/master/doc/overview.png

def line(color):
    return {'linestyle': '-', 'color': color, 'zorder': 1, 'linewidth': 2}

def titlestr(t):
    d = int(t / day)                 # The number of days
    h = int(int(t % day) / hour)     # The number of remaining hours
    m = int(int(t % hour) / minute)  # The number of remaining minutes
    return '{:>3d}d {:>2}h {:>2}m'.format(int(d), str(int(h)).zfill(2), str(int(m)).zfill(2))

def plot_figures_ph(params):

    plot_at_selected_steps = params["plot_at_selected_steps"]
    dt = params["dt"]
    folder = params["folder"]
    files = params["files"]
    indx_ph = params["indx_ph"]
    xcells = params["xcells"]

    for i in plot_at_selected_steps:
        print("On pH figure at time step: {}".format(i))
        t = i * dt
        filearray = np.loadtxt(folder + '/' + files[i-1], skiprows=1)
        data = filearray.T
        data_ph = data[indx_ph]
        plt.axes(xlim=(-0.01, 1.001), ylim=(8.5, 10.0))
        plt.xlabel('Distance [m]')
        plt.ylabel('pH')
        plt.title(titlestr(t))
        plt.plot(xcells, data_ph, label='pH', **line('teal'))
        plt.legend(loc='lower right')
        plt.savefig('figures/ph/pH-{}.pdf'.format(i))
        plt.tight_layout()
        plt.close()

def plot_figures_elements(params):

    plot_at_selected_steps = params["plot_at_selected_steps"]
    dt = params["dt"]
    folder = params["folder"]
    files = params["files"]
    xcells = params["xcells"]

    indices = [params["indx_C"], params["indx_Cl"], params["indx_Ca"], params["indx_Fe"], params["indx_K"],
               params["indx_Mg"], params["indx_Na"], params["indx_S"], params["indx_Si"], params["indx_Al"]]
    labels = ['C', 'Cl', 'Ca', 'Fe', 'K', 'Mg', 'Na', 'S', 'Si', 'Al']

    """
    indx_C = params["indx_C"]
    indx_Cl = params["indx_Cl"]
    indx_Ca = params["indx_Ca"]
    indx_Fe = params["indx_Fe"]
    indx_K = params["indx_K"]
    indx_Mg = params["indx_Mg"]
    indx_Na = params["indx_Na"]
    indx_S = params["indx_S"]
    indx_Si = params["indx_Si"]
    indx_Al = params["indx_Al"]
    indx_Fe = params["indx_S"]
    """
    colors = ['grey', 'tan', 'gold', 'yellowgreen', 'palevioletred',
              'steelblue', 'salmon', 'ligthseagreen', 'firebrick', 'olivedrab']
    xcells = params["xcells"]

    for i in plot_at_selected_steps:
        print("On pH figure at time step: {}".format(i))
        t = i * dt
        filearray = np.loadtxt(folder + '/' + files[i-1], skiprows=1)
        data = filearray.T
        """
        data_C = data[indx_C]
        data_Cl = data[indx_Cl]
        data_Ca = data[indx_Ca]
        data_Fe = data[indx_Fe]
        data_K = data[indx_K]
        data_Mg = data[indx_Mg]
        data_Na = data[indx_Na]
        data_S = data[indx_S]
        data_Si = data[indx_Si]
        data_Al = data[indx_Al]
        data_Fe = data[indx_Fe]
        """
        #plt.axes(xlim=(-0.01, 1.001), ylim=(8.5, 10.0))
        plt.axes(xlim=(-0.01, 1.001))
        plt.xlabel('Distance [m]')
        plt.ylabel('Elements amounts [molal]')
        plt.title(titlestr(t))
        for j in range(0, len(indices)):
            plt.plot(xcells, data[indices[j]], label=labels[j], **line(colors[j]))
        plt.legend(loc='lower right')
        plt.savefig('figures/elements/elements-{}.pdf'.format(i))
        plt.tight_layout()
        plt.show()
        plt.close()

def plot_figures_aqueous_species(params):

    plot_at_selected_steps = params["plot_at_selected_steps"]
    dt = params["dt"]
    folder = params["folder"]
    files = params["files"]
    xcells = params["xcells"]

    indices = [params["indx_Hcation"], params["indx_Clanoion"], params["indx_Cacation"], 
                params["indx_Mgcation"], params["indx_Nacation"], params["indx_HCO3anoion"], 
                params["indx_CO2aq"]]
    labels = ['H+', 'Cl-', 'Ca+2', 'Na+', 'HCO3-', 'CO2(aq)']
    colors = ['darkcyan', 'mediumaqumarine', 'aqumarine', 
              'bisque', 'burlywood', 'sienna', 'saddlebrown']
    
    xcells = params["xcells"]

    for i in plot_at_selected_steps:
        print("On pH figure at time step: {}".format(i))
        t = i * dt
        filearray = np.loadtxt(folder + '/' + files[i-1], skiprows=1)
        data = filearray.T
        #plt.axes(xlim=(-0.01, 1.001), ylim=(8.5, 10.0))
        plt.axes(xlim=(-0.01, 1.001))
        plt.xlabel('Distance [m]')
        plt.ylabel('Species amounts [molal]')
        plt.title(titlestr(t))
        for j in range(0, len(indices)):
            plt.plot(xcells, data[indices[j]], label=labels[j], **line(colors[j]))
        plt.legend(loc='lower right')
        plt.savefig('figures/aqueous_species/aqueous-species-{}.pdf'.format(i))
        plt.tight_layout()
        plt.show()
        plt.close()

def plot_figures_minerals(params):

    plot_at_selected_steps = params["plot_at_selected_steps"]
    dt = params["dt"]
    folder = params["folder"]
    files = params["files"]
    xcells = params["xcells"]

    indices = [params["indx_Hcation"], params["indx_Clanoion"], params["indx_Cacation"], 
                params["indx_Mgcation"], params["indx_Nacation"], params["indx_HCO3anoion"], 
                params["indx_CO2aq"]]
    labels = ['Calcite', 'Hydrotalcite', 'Portlandite', 'C4AH11', 'CSHQJenD', 'CSHQJenH', 
              'CSHQTobD', 'CSHQTobH', 'C3AFS', 'Brucite', 'Ettringite03_ss', 'Ettringite13', 'Ettringite9']
    colors = ['khaki', 'gold', 'orange', 'orangered', 'crimson', 'darkmagenta', 'indigo', 
              'black', 'darkblue', 'royalblue', 'cornflowerblue', 'lightsteelblue', 'slategrey']
    
    xcells = params["xcells"]

    for i in plot_at_selected_steps:
        print("On pH figure at time step: {}".format(i))
        t = i * dt
        filearray = np.loadtxt(folder + '/' + files[i-1], skiprows=1)
        data = filearray.T
        #plt.axes(xlim=(-0.01, 1.001), ylim=(8.5, 10.0))
        plt.axes(xlim=(-0.01, 1.001))
        plt.xlabel('Distance [m]')
        plt.ylabel('Species amounts [molal]')
        plt.title(titlestr(t))
        for j in range(0, len(indices)):
            plt.plot(xcells, data[indices[j]], label=labels[j], **line(colors[j]))
        plt.legend(loc='lower right')
        plt.savefig('figures/minerals/minerals-{}.pdf'.format(i))
        plt.tight_layout()
        plt.show()
        plt.close()

def plot_animation_ph(params):

    plot_at_selected_steps = params["plot_at_selected_steps"]
    dt = params["dt"]
    folder = params["folder"]
    files = params["files"]
    indx_ph = params["indx_ph"]
    xcells = params["xcells"]

    animation_starts_at_frame = params["animation_starts_at_frame"]  # the first frame index to be considered
    animation_ends_at_frame = params["animation_ends_at_frame"] # the last frame index to be considered
    animation_num_frames_to_jump = params["animation_num_frames_to_jump"] # the number of frames to jump between current and next
    animation_fps = params["animation_fps"] # the number of frames per second
    animation_interval_wait = params["animation_interval_wait"] # the time (in milliseconds) to wait between each frame

    # Auxiliary animation options
    animation_frame_range = range(animation_starts_at_frame, animation_ends_at_frame, animation_num_frames_to_jump)

    # Plot of mineral's volume the space coordinates
    fig = plt.figure()
    ax = plt.axes(xlim=(-0.01, 1.001), ylim=(8.5, 10.0))
    ax.set_xlabel('Distance [m]')
    ax.set_ylabel('pH')
    ax.set_title(titlestr(0.0))
    objects = [
        ax.plot([], [], label='pH', **line('teal'))[0],
    ]
    ax.legend(loc='lower right')

    def init():
        return tuple(objects)


    def animate(i):
        print("On pH animation index: {}".format(i))
        t = i * dt
        filearray = np.loadtxt(folder + '/' + files[i - 1], skiprows=1)
        data = filearray.T
        data_ph = data[indx_ph]
        objects[0].set_data(xcells, data_ph)
        ax.set_title(titlestr(t))
        return tuple(objects)


    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=animation_frame_range, interval=animation_interval_wait, blit=True)

    anim.save('videos/pH.mp4', fps=animation_fps, extra_args=['-vcodec', 'libx264'])

