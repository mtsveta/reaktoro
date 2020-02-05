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
        plt.axes(xlim=(-0.01, 2.001), ylim=(8.5, 10.0))
        plt.xlabel('Distance [m]')
        plt.ylabel('pH')
        plt.title(titlestr(t))
        plt.plot(xcells, data_ph, label='pH', **line('teal'))
        plt.legend(loc='lower right')
        plt.savefig('figures/ph/pH-{}.pdf'.format(i))
        plt.tight_layout()
        plt.close()

def plot_figures_Eh(params):

    plot_at_selected_steps = params["plot_at_selected_steps"]
    dt = params["dt"]
    folder = params["folder"]
    files = params["files"]
    indx_Eh = params["indx_Eh"]
    xcells = params["xcells"]

    for i in plot_at_selected_steps:
        print("On Eh figure at time step: {}".format(i))
        t = i * dt
        filearray = np.loadtxt(folder + '/' + files[i - 1], skiprows=1)
        data = filearray.T
        data_ph = data[indx_Eh]
        plt.axes(xlim=(-0.01, 2.001), ylim=(-0.5, 0.75))
        plt.xlabel('Distance [m]')
        plt.ylabel('Eh')
        plt.title(titlestr(t))
        plt.plot(xcells, data_ph, label='Eh', **line('darkseagreen'))
        plt.legend(loc='lower right')
        plt.savefig('figures/Eh/Eh-{}.pdf'.format(i))
        plt.tight_layout()
        plt.close()

def plot_figures_SI(params):

    plot_at_selected_steps = params["plot_at_selected_steps"]
    dt = params["dt"]
    folder = params["folder"]
    files = params["files"]
    indx_SI = params["indx_SI"]
    xcells = params["xcells"]

    for i in plot_at_selected_steps:
        print("On SI figure at time step: {}".format(i))
        t = i * dt
        filearray = np.loadtxt(folder + '/' + files[i - 1], skiprows=1)
        data = filearray.T
        data_SI = data[indx_SI]
        plt.axes(xlim=(-0.01, 2.001), ylim=(0, 2.5))
        plt.xlabel('Distance [m]')
        plt.ylabel('SI')
        plt.title(titlestr(t))
        plt.plot(xcells, data_SI, label='SI', **line('coral'))
        plt.legend(loc='lower right')
        plt.savefig('figures/SI/SI-{}.pdf'.format(i))
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

    # C, Cl, Fe, Mg, Na, Al
    small_range_element = [1, 1, 0, 1, 0, 1, 1, 0, 0, 1]

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
              'steelblue', 'salmon', 'lightseagreen', 'firebrick', 'olivedrab']
    for i in plot_at_selected_steps:
        print("On elements figure at time step: {}".format(i))
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
        plt.axes(xlim=(-0.01, 2.001), ylim=(1e-7, 1e1))
        #plt.axes(xlim=(-0.01, 2.001))
        plt.xlabel('Distance [m]')
        plt.ylabel('Elements molality [molal]')
        plt.yscale('log')
        plt.title(titlestr(t))
        for j in range(0, len(indices)):
            #if small_range_element[j]:
            plt.plot(xcells, data[indices[j]], label=labels[j], **line(colors[j]))
        plt.legend(loc='upper right')
        plt.savefig('figures/elements/elements-small-range-{}.pdf'.format(i))
        plt.tight_layout()
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
    labels = ['H+', 'Cl-', 'Ca+2', 'Mg+2', 'Na+', 'HCO3-', 'CO2(aq)']
    colors = ['darkcyan', 'mediumaquamarine', 'aquamarine',
              'bisque', 'burlywood', 'sienna', 'indianred']
    
    for i in plot_at_selected_steps:
        print("On aqueous species figure at time step: {}".format(i))
        t = i * dt
        filearray = np.loadtxt(folder + '/' + files[i-1], skiprows=1)
        data = filearray.T

        plt.axes(xlim=(-0.01, 2.001), ylim=(1e-11, 1e0))
        plt.xlabel('Distance [m]')
        plt.ylabel('Species molality [molal]')
        plt.yscale('log')
        plt.title(titlestr(t))
        for j in range(0, len(indices)):
            plt.plot(xcells, data[indices[j]], label=labels[j], **line(colors[j]))
        plt.legend(loc='lower right')
        plt.savefig('figures/aqueous_species/aqueous-species-{}.pdf'.format(i))
        plt.tight_layout()
        plt.close()

def plot_figures_solids(params):

    plot_at_selected_steps = params["plot_at_selected_steps"]
    dt = params["dt"]
    folder = params["folder"]
    files = params["files"]
    xcells = params["xcells"]

    indices = [params["indx_Cal"], params["indx_hydrotalcite"], params["indx_Portlandite"],  params["indx_C4AH11"],
               params["indx_CSHQ"], params["indx_C3AFS"], params["indx_Brc"], params["indx_ettringite"]]
    labels = ['Calcite', 'Hydrotalcite', 'Portlandite', 'C4AH11', 'CSHQ', 'C3AFS', 'Brucite', 'Ettringite']
    colors = ['lightsteelblue', 'orange', 'crimson', 'darkmagenta', 'indigo',
              'darkblue', 'cornflowerblue', 'khaki']
    
    for i in plot_at_selected_steps:
        print("On solids figure at time step: {}".format(i))
        t = i * dt
        filearray = np.loadtxt(folder + '/' + files[i-1], skiprows=1)
        data = filearray.T
        plt.axes(xlim=(-0.01, 2.001), ylim=(1e-26, 1e-4))
        #plt.axes(xlim=(-0.01, 2.001))
        plt.xlabel('Distance [m]')
        plt.ylabel('Mineral volum [%vol]')
        plt.yscale('log')
        plt.title(titlestr(t))
        for j in range(0, len(indices)):
            plt.plot(xcells, data[indices[j]], label=labels[j], **line(colors[j]))
        plt.legend(loc='center right')
        plt.savefig('figures/minerals/minerals-{}.pdf'.format(i))
        plt.tight_layout()
        plt.close()

def plot_animation_ph(params):

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
    ax = plt.axes(xlim=(-0.01, 2.001), ylim=(8.5, 10.0))
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
        print("On pH animation on step {}".format(i))
        t = i * dt
        filearray = np.loadtxt(folder + '/' + files[i - 1], skiprows=1)
        data = filearray.T
        data_ph = data[indx_ph]
        objects[0].set_data(xcells, data_ph)
        ax.set_title(titlestr(t))
        return tuple(objects)


    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=animation_frame_range, interval=animation_interval_wait, blit=True)

    anim.save('videos/pH.mp4', fps=animation_fps, extra_args=['-vcodec', 'libx264'])

def plot_animation_Eh(params):

    dt = params["dt"]
    folder = params["folder"]
    files = params["files"]
    indx_Eh = params["indx_Eh"]
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
    ax = plt.axes(xlim=(-0.01, 2.001), ylim=(-0.5, 0.75))
    ax.set_xlabel('Distance [m]')
    ax.set_ylabel('Eh')
    ax.set_title(titlestr(0.0))
    objects = [
        ax.plot([], [], label='Eh', **line('darkseagreen'))[0],
    ]
    ax.legend(loc='lower right')

    def init():
        return tuple(objects)


    def animate(i):
        print("On Eh animation on step {}".format(i))
        t = i * dt
        filearray = np.loadtxt(folder + '/' + files[i - 1], skiprows=1)
        data = filearray.T
        data_Eh = data[indx_Eh]
        objects[0].set_data(xcells, data_Eh)
        ax.set_title(titlestr(t))
        return tuple(objects)


    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=animation_frame_range, interval=animation_interval_wait, blit=True)

    anim.save('videos/Eh.mp4', fps=animation_fps, extra_args=['-vcodec', 'libx264'])

def plot_animation_SI(params):

    dt = params["dt"]
    folder = params["folder"]
    files = params["files"]
    indx_SI = params["indx_SI"]
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
    ax = plt.axes(xlim=(-0.01, 2.001), ylim=(0, 2.5))
    ax.set_xlabel('Distance [m]')
    ax.set_ylabel('SI')
    ax.set_title(titlestr(0.0))
    objects = [
        ax.plot([], [], label='SI', **line('coral'))[0],
    ]
    ax.legend(loc='lower right')

    def init():
        return tuple(objects)


    def animate(i):
        print("On SI animation on step {}".format(i))
        t = i * dt
        filearray = np.loadtxt(folder + '/' + files[i - 1], skiprows=1)
        data = filearray.T
        data_SI = data[indx_SI]
        objects[0].set_data(xcells, data_SI)
        ax.set_title(titlestr(t))
        return tuple(objects)


    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=animation_frame_range, interval=animation_interval_wait, blit=True)

    anim.save('videos/SI.mp4', fps=animation_fps, extra_args=['-vcodec', 'libx264'])

def plot_animation_elements(params):

    dt = params["dt"]
    folder = params["folder"]
    files = params["files"]
    xcells = params["xcells"]

    indices = [params["indx_C"], params["indx_Cl"], params["indx_Ca"], params["indx_Fe"], params["indx_K"],
               params["indx_Mg"], params["indx_Na"], params["indx_S"], params["indx_Si"], params["indx_Al"]]
    labels = ['C', 'Cl', 'Ca', 'Fe', 'K', 'Mg', 'Na', 'S', 'Si', 'Al']

    colors = ['grey', 'tan', 'gold', 'yellowgreen', 'palevioletred',
              'steelblue', 'salmon', 'lightseagreen', 'firebrick', 'olivedrab']

    animation_starts_at_frame = params["animation_starts_at_frame"]  # the first frame index to be considered
    animation_ends_at_frame = params["animation_ends_at_frame"] # the last frame index to be considered
    animation_num_frames_to_jump = params["animation_num_frames_to_jump"] # the number of frames to jump between current and next
    animation_fps = params["animation_fps"] # the number of frames per second
    animation_interval_wait = params["animation_interval_wait"] # the time (in milliseconds) to wait between each frame

    # Auxiliary animation options
    animation_frame_range = range(animation_starts_at_frame, animation_ends_at_frame, animation_num_frames_to_jump)

    # Plot of mineral's volume the space coordinates
    fig = plt.figure()
    ax = plt.axes(xlim=(-0.01, 2.001), ylim=(1e-7, 1e1))
    ax.set_xlabel('Distance [m]')
    ax.set_ylabel('Elements molality [molal]')
    plt.yscale('log')
    ax.set_title(titlestr(0.0))
    objects = [
        ax.plot([], [], label=labels[0], **line(colors[0]))[0],
        ax.plot([], [], label=labels[1], **line(colors[1]))[0],
        ax.plot([], [], label=labels[2], **line(colors[2]))[0],
        ax.plot([], [], label=labels[3], **line(colors[3]))[0],
        ax.plot([], [], label=labels[4], **line(colors[4]))[0],
        ax.plot([], [], label=labels[5], **line(colors[5]))[0],
        ax.plot([], [], label=labels[6], **line(colors[6]))[0],
        ax.plot([], [], label=labels[7], **line(colors[7]))[0],
        ax.plot([], [], label=labels[8], **line(colors[8]))[0],
        ax.plot([], [], label=labels[9], **line(colors[9]))[0],
    ]
    ax.legend(loc='lower right')

    def init(): return tuple(objects)

    def animate(i):
        print("On elements on step {}".format(i))
        t = i * dt
        filearray = np.loadtxt(folder + '/' + files[i - 1], skiprows=1)
        data = filearray.T
        for j in range(0, len(indices)):
            objects[j].set_data(xcells, data[indices[j]])
        ax.set_title(titlestr(t))
        return tuple(objects)

    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=animation_frame_range, interval=animation_interval_wait, blit=True)
    anim.save('videos/elements.mp4', fps=animation_fps, extra_args=['-vcodec', 'libx264'])

def plot_animation_aqueous_species(params):

    dt = params["dt"]
    folder = params["folder"]
    files = params["files"]
    xcells = params["xcells"]

    indices = [params["indx_Hcation"], params["indx_Clanoion"], params["indx_Cacation"],
               params["indx_Mgcation"], params["indx_Nacation"], params["indx_HCO3anoion"],
               params["indx_CO2aq"]]
    labels = ['H+', 'Cl-', 'Ca+2', 'Mg+2', 'Na+', 'HCO3-', 'CO2(aq)']
    colors = ['darkcyan', 'mediumaquamarine', 'aquamarine',
              'bisque', 'burlywood', 'sienna', 'indianred']

    animation_starts_at_frame = params["animation_starts_at_frame"]  # the first frame index to be considered
    animation_ends_at_frame = params["animation_ends_at_frame"] # the last frame index to be considered
    animation_num_frames_to_jump = params["animation_num_frames_to_jump"] # the number of frames to jump between current and next
    animation_fps = params["animation_fps"] # the number of frames per second
    animation_interval_wait = params["animation_interval_wait"] # the time (in milliseconds) to wait between each frame

    # Auxiliary animation options
    animation_frame_range = range(animation_starts_at_frame, animation_ends_at_frame, animation_num_frames_to_jump)

    # Plot of mineral's volume the space coordinates
    fig = plt.figure()
    ax =  plt.axes(xlim=(-0.01, 2.001), ylim=(1e-11, 1e0))
    ax.set_xlabel('Distance [m]')
    ax.set_ylabel('Species molality [molal]')
    plt.yscale('log')
    ax.set_title(titlestr(0.0))
    objects = [
        ax.plot([], [], label=labels[0], **line(colors[0]))[0],
        ax.plot([], [], label=labels[1], **line(colors[1]))[0],
        ax.plot([], [], label=labels[2], **line(colors[2]))[0],
        ax.plot([], [], label=labels[3], **line(colors[3]))[0],
        ax.plot([], [], label=labels[4], **line(colors[4]))[0],
        ax.plot([], [], label=labels[5], **line(colors[5]))[0],
        ax.plot([], [], label=labels[6], **line(colors[6]))[0],
    ]
    ax.legend(loc='lower right')

    def init(): return tuple(objects)

    def animate(i):
        print("On elements on step {}".format(i))
        t = i * dt
        filearray = np.loadtxt(folder + '/' + files[i - 1], skiprows=1)
        data = filearray.T
        for j in range(0, len(indices)):
            objects[j].set_data(xcells, data[indices[j]])
        ax.set_title(titlestr(t))
        return tuple(objects)

    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=animation_frame_range, interval=animation_interval_wait, blit=True)
    anim.save('videos/species.mp4', fps=animation_fps, extra_args=['-vcodec', 'libx264'])

def plot_animation_solids(params):

    dt = params["dt"]
    folder = params["folder"]
    files = params["files"]
    xcells = params["xcells"]

    indices = [params["indx_Cal"], params["indx_hydrotalcite"], params["indx_Portlandite"], params["indx_C4AH11"],
               params["indx_CSHQ"], params["indx_C3AFS"], params["indx_Brc"], params["indx_ettringite"]]
    labels = ['Calcite', 'Hydrotalcite', 'Portlandite', 'C4AH11', 'CSHQ', 'C3AFS', 'Brucite', 'Ettringite']
    colors = ['lightsteelblue', 'orange', 'crimson', 'darkmagenta', 'indigo',
              'darkblue', 'cornflowerblue', 'khaki']

    animation_starts_at_frame = params["animation_starts_at_frame"]  # the first frame index to be considered
    animation_ends_at_frame = params["animation_ends_at_frame"] # the last frame index to be considered
    animation_num_frames_to_jump = params["animation_num_frames_to_jump"] # the number of frames to jump between current and next
    animation_fps = params["animation_fps"] # the number of frames per second
    animation_interval_wait = params["animation_interval_wait"] # the time (in milliseconds) to wait between each frame

    # Auxiliary animation options
    animation_frame_range = range(animation_starts_at_frame, animation_ends_at_frame, animation_num_frames_to_jump)

    # Plot of mineral's volume the space coordinates
    fig = plt.figure()
    ax =  plt.axes(xlim=(-0.01, 2.001), ylim=(1e-26, 1e-4))
    ax.set_xlabel('Distance [m]')
    plt.ylabel('Mineral volum [%vol]')
    plt.yscale('log')
    ax.set_title(titlestr(0.0))
    objects = [
        ax.plot([], [], label=labels[0], **line(colors[0]))[0],
        ax.plot([], [], label=labels[1], **line(colors[1]))[0],
        ax.plot([], [], label=labels[2], **line(colors[2]))[0],
        ax.plot([], [], label=labels[3], **line(colors[3]))[0],
        ax.plot([], [], label=labels[4], **line(colors[4]))[0],
        ax.plot([], [], label=labels[5], **line(colors[5]))[0],
        ax.plot([], [], label=labels[6], **line(colors[6]))[0],
        ax.plot([], [], label=labels[7], **line(colors[7]))[0],
    ]
    ax.legend(loc='lower right')

    def init(): return tuple(objects)

    def animate(i):
        print("On solids on step {}".format(i))
        t = i * dt
        filearray = np.loadtxt(folder + '/' + files[i - 1], skiprows=1)
        data = filearray.T
        for j in range(0, len(indices)):
            objects[j].set_data(xcells, data[indices[j]])
        ax.set_title(titlestr(t))
        return tuple(objects)

    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=animation_frame_range, interval=animation_interval_wait, blit=True)
    anim.save('videos/solids.mp4', fps=animation_fps, extra_args=['-vcodec', 'libx264'])
