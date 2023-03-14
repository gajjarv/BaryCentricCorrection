import matplotlib.pyplot as plt
from blimpy import Waterfall as wf
import os
import matplotlib
matplotlib.use('MacOSX')
import subprocess
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

def bary_command(TEL, RA, DEC, Freq, fname):
    command = f"bary {TEL} {RA} {DEC} {Freq} < {fname} | awk '{{printf \"%.15f \",$3}}'"
    relvel_list = subprocess.check_output(command, shell=True).decode().strip().split()
    relvel_list = np.array(relvel_list, dtype=np.float64)
    obs_freq_bary = Freq * (1 - relvel_list)
    return obs_freq_bary


def plot_data_on_axis(data, time, axis, name):
    import matplotlib.pyplot as plt
    # convert the data from MHz to Hz
    data_hz = [num * 1e6 for num in data]
    datamean = np.mean(data_hz)
    datap = data_hz - datamean
    # set x-axis tick positions and labels
    nstep = 5
    ticks = np.linspace(-max(abs(datap)), max(abs(datap)), nstep)
    tick_labels = [f'{tick:.0f}' for tick in ticks]

    # plot the data with thicker border
    axis.plot(datap, time, linewidth=2.5, color='b', zorder=1)

    # set x-axis tick positions and labels with bigger font size
    axis.set_xticks(ticks)
    axis.set_xticklabels(tick_labels, fontsize=12,fontweight='bold')

    # set x-axis limits
    axis.set_xlim(-max(abs(datap)), max(abs(datap)))
    axis.set_ylim(0, max(time))

    axis.set_title(name, fontsize=16)

    # set x-axis label with bigger font size
    axis.set_xlabel('Frequency (Hz) - %f MHz' % (np.mean(data)), fontsize=12,fontweight='bold')

    # set tick label size
    axis.tick_params(axis='both', which='major', labelsize=12, width=1.5, length=6, pad=6)

    # draw a thicker border around the plot
    for axis_edge in ['left', 'bottom', 'right', 'top']:
        axis.spines[axis_edge].set_linewidth(2)

    num_points = 50
    noise = np.random.normal(size=(num_points, num_points))
    #extent = [-max(abs(datap)), max(abs(datap)), time[0], time[-1]]
    axis.imshow(noise,extent=[-max(abs(datap)), max(abs(datap)),0,max(time)], cmap='gray', alpha=0.2, aspect='auto', origin='upper')

    return axis

#ProximaSen
#RA = "14:29:42.9"
#DEC = "-61:59:53.8"
#Freq = "982.0"

RA="04:25:28.83"
DEC="46:21:57.25"
Freq=1420.0
#TEL="PK"
TEL="IL"

nVelInc = []
nVelIncMJD = 59054.7

nVelDec = []
nVelDecMJD = 58775.9

pVelInc = []
pVelIncMJD = 58885.7

pVelDec = []
pVelDecMJD = 58930.1

for i in range(16):
    nVelIncMJD += 0.00019418074; nVelInc.append(nVelIncMJD)
    nVelDecMJD += 0.00019418074; nVelDec.append(nVelDecMJD)
    pVelIncMJD += 0.00019418074; pVelInc.append(pVelIncMJD)
    pVelDecMJD += 0.00019418074; pVelDec.append(pVelDecMJD)

arrays = {"nVelInc": nVelInc, "nVelDec": nVelDec, "pVelInc": pVelInc, "pVelDec": pVelDec}

for name, array in arrays.items():
    filename = f"input_{name}"
    np.savetxt(filename, array, delimiter="\n",fmt='%.11f')

nVelInc_bary=bary_command(TEL,RA, DEC, Freq,"input_nVelInc")
nVelDec_bary=bary_command(TEL,RA, DEC, Freq,"input_nVelDec")
pVelInc_bary=bary_command(TEL,RA, DEC, Freq,"input_pVelInc")
pVelDec_bary=bary_command(TEL,RA, DEC, Freq,"input_pVelDec")

time = [i*16.777216 for i in range(0,16)]

fig, axs = plt.subplots(nrows=1, ncols=4,figsize=(16, 4),constrained_layout=True)

# plot each data array on a separate subplot
print(time)
print(pVelInc_bary)
plot_data_on_axis(pVelInc_bary,time, axs[0],"$+v_{rel}$\u2191")
axs[0].set_ylabel('Time (seconds)',fontweight='bold',fontsize=12)
plot_data_on_axis(pVelDec_bary,time, axs[1],"$+v_{rel}$\u2193")
plot_data_on_axis(nVelInc_bary,time, axs[2],"$-v_{rel}$\u2191")
plot_data_on_axis(nVelDec_bary,time, axs[3],"$-v_{rel}$\u2193")

plt.savefig("Barycentric_frequency_drift.pdf")

# display the plot
plt.show()

