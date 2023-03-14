import pylab as plt
from blimpy import Waterfall as wf
import os 
import matplotlib
matplotlib.use('MacOSX')

file1 = "nVelInc.fil"
nVelInc_bary = [981.9751610,	981.9751620,	981.9751630,	981.9751640,	981.9751650,	981.9751660,	981.9751670,	981.9751680,	981.9751690,	981.9751700,	981.9751710,	981.9751710,	981.9751720,	981.9751730,	981.9751740,	981.9751750]

time = [i*16.777216 for i in range(0,16)]
obs=wf(file1)
plt.figure(1)
obs.plot_waterfall(f_start=981.9751,f_stop=981.9752)
plt.plot(nVelInc_bary,time,color='red')

file2 = "nVelDec.fil"

nVelDec_bary = [981.9758250,	981.9758240,	981.9758240,	981.9758230,	981.9758230,	981.9758220,	981.9758220,	981.9758210,	981.9758210,	981.9758200,	981.9758190,	981.9758190,	981.9758180,	981.9758180,	981.9758170,	981.9758170]

plt.figure(2)
obs1=wf(file2)
obs1.plot_waterfall(f_start=981.97576,f_stop=981.9759)
plt.plot(nVelDec_bary,time,color='red')

file3 = "pVelInc.fil"

pValInc_bary = [982.0100700,	982.0100710,	982.0100730,	982.0100750,	982.0100770,	982.0100790,	982.0100800,	982.0100820,	982.0100840,	982.0100860,	982.0100880,	982.0100900,	982.0100910,	982.0100930,	982.0100950,	982.0100970]

plt.figure(3)
obs3=wf(file3)
obs3.plot_waterfall(f_start=982.0100,f_stop=982.0102)
plt.plot(pValInc_bary,time,color='red')

file4 = "pVelDec.fil"

pVelDec_bary = [982.0110370,	982.0110360,	982.0110350,	982.0110340,	982.0110320,	982.0110310,	982.0110300,	982.0110290,	982.0110270,	982.0110260,	982.0110250,	982.0110240,	982.0110230,	982.0110210,	982.0110200,	982.0110190]

plt.figure(4)
obs4=wf(file4)
obs4.plot_waterfall(f_start=982.0109,f_stop=982.0111)
plt.plot(pVelDec_bary,time,color='red')

plt.show()
