import matplotlib.pyplot as plt
from blimpy import Waterfall as wf
import os
import matplotlib
matplotlib.use('MacOSX')
import subprocess
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


def dec_ddmmss_to_deg(dec_ddmmss):
    sign = 1 if dec_ddmmss[0] == "+" else -1
    dec_deg = sign * (
        int(dec_ddmmss[1:3])
        + int(dec_ddmmss[4:6]) / 60.0
        + float(dec_ddmmss[7:]) / 3600.0
    )
    return dec_deg

def bary_vel_plot(DEC, MJDarr, ax):
    RA = "04:25:28.83"
    Freq = "1420.0"
    TEL = "IL"

    DEC_deg = dec_ddmmss_to_deg(DEC)

    np.savetxt("input", MJDarr, delimiter="\n", fmt="%.11f")
    command = "bary %s %s %s %s < input | awk '{printf \"%%.7f \",$4}'" % (
        TEL,
        RA,
        DEC_deg,
        Freq,
    )
    output = subprocess.check_output(command, shell=True)
    output_str = output.decode("utf-8").strip()
    output_list = output_str.split(" ")
    data = np.array(output_list, dtype=np.float64)
    data = np.array(data)

    ax.plot(MJDarr, data, label=f"DEC={DEC}")
    ax.legend()

    return
RA = "04:25:28.83"
Freq = "1420"
# TEL="PK"
TEL = "IL"

fig, axs = plt.subplots()
MJDarr = np.arange(58600.0, 58965.0)

for DEC in ["-20:00:00", "-10:00:00", "00:00:00", "+10:00:00", "+20:00:00"]:
    bary_vel_plot(DEC, MJDarr, axs)

plt.show()

