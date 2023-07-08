
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import mark_inset, inset_axes
import os
print(os.getcwd())

# This Python script allows us to show all the different properties and parameters for the Lennard-Jones proyect.
# To show the desired property or parameter, the rest of the subplots are framed in /* */

#------------------------------------------------------------
'''
with open('T_todos.txt') as f:
    lines = f.readlines()[0:]

MS = [0]
overlap = [0]


for line in lines:
    values = line.strip().split(' ')
    MS.append(float(values[0]))
    overlap.append(float(values[1]))
    
with open('overlapping.txt') as f:
    lines = f.readlines()[0:]

v_MS = [0]
overlap_1 = [0]
overlap_2 = [0]
overlap_3 = [0]
overlap_4 = [0]


for line in lines:
    values = line.strip().split(' ')
    v_MS.append(float(values[0]))
    overlap_1.append(float(values[1]))
    overlap_2.append(float(values[2]))
    overlap_3.append(float(values[3]))
    overlap_4.append(float(values[4]))
    
with open('varias_T.txt') as f:
    lines = f.readlines()[0:]

T = []
T1 = []
T2 = []
T3 = []
T4 = []


for line in lines:
    values = line.strip().split(' ')
    T.append(float(values[0]))
    T1.append(float(values[1]))
    T2.append(float(values[2]))
    T3.append(float(values[3]))   
    T4.append(float(values[4]))   
'''

with open('varios_datos.txt') as f:
    lines = f.readlines()[0:]

tv = []
s0 = []
s1 = []
s2 = []
s3 = []
s4 = []
s5 = []
s6 = []
s7 = []
s8 = []
s9 = []
s10 = []
    
for line in lines:
    values = line.strip().split(' ')
    tv.append(float(values[0]))
    s0.append(float(values[1]))
    s1.append(float(values[2]))
    s2.append(float(values[3]))
    s3.append(float(values[4]))
    s4.append(float(values[5]))
    s5.append(float(values[6]))
    s6.append(float(values[7]))
    s7.append(float(values[8]))
    s8.append(float(values[9]))
    s9.append(float(values[10]))
    s10.append(float(values[11]))
    
#------------------------------------------------------------

fig, ax = plt.subplots(1, 1, figsize=(20,14))

plt.style.use('default')

plt.rc('xtick', labelsize=25)
plt.rc('ytick', labelsize=25)

gs = gridspec.GridSpec(1, 1)
scatter_legend_size = 100

'''
#---Subplot 1------------
ax1 = plt.subplot(gs[0, 0])
ax1.tick_params(axis='x', direction='out', pad=15)
ax1.tick_params(axis='y', direction='out', pad=15)
plt.axhline(y=0, linestyle='--', linewidth=2, color='grey', dashes=(5, 5), alpha=0.6)
plt.axhline(y=1, linestyle='--', linewidth=2, color='grey', dashes=(5, 5), alpha=0.6)
plt.axhline(y=-1, linestyle='--', linewidth=2, color='grey', dashes=(5, 5), alpha=0.6)
plt.plot(MS, overlap, linewidth=5, label = 'm$^1$ para T = 10$^{-4}$ K')
plt.xlabel('Pasos de Montecarlo (PM)', fontsize=25, labelpad=20)
plt.ylabel('Solapamiento con $Salchicha$', fontsize=25, labelpad=15)
legend = plt.legend(loc=0, bbox_to_anchor=(0.95, 0.9), prop={'size': 25})
for legend_handle in legend.legendHandles:
    legend_handle._sizes = [scatter_legend_size]
#ax1.set_xlim(0, 4)
ax1.set_ylim(-1.2, 1.2)
plt.grid()
'''
'''
#---Subplot 1------------
ax1 = plt.subplot(gs[0, 0])
ax1.tick_params(axis='x', direction='out', pad=15)
ax1.tick_params(axis='y', direction='out', pad=15)
plt.axhline(y=0, linestyle='--', linewidth=2, color='grey', dashes=(5, 5), alpha=0.6)
plt.axhline(y=1, linestyle='--', linewidth=2, color='grey', dashes=(5, 5), alpha=0.6)
plt.axhline(y=-1, linestyle='--', linewidth=2, color='grey', dashes=(5, 5), alpha=0.6)
plt.plot(v_MS, overlap_1, linewidth=5, label = 'm$^1$')
plt.plot(v_MS, overlap_2, linewidth=5, label = 'm$^2$', color='#87138A')
plt.plot(v_MS, overlap_3, linewidth=5, label = 'm$^3$', color='#088F8F')
plt.plot(v_MS, overlap_4, linewidth=5, label = 'm$^4$', color='#89CFF0')
plt.xlabel('Pasos de Montecarlo (PM) para T = 10$^{-4}$ K', fontsize=25, labelpad=20)
plt.ylabel('Solapamiento de cada patrón', fontsize=25, labelpad=15)
legend = plt.legend(loc=0, prop={'size': 25})
for legend_handle in legend.legendHandles:
    legend_handle._sizes = [scatter_legend_size]
#ax1.set_xlim(0, 4)
ax1.set_ylim(-1.2, 1.2)
plt.grid()

axins = inset_axes(ax1, 2.5, 2.5, loc=1, bbox_to_anchor=(0.5, 0.95), bbox_transform=ax1.figure.transFigure)

axins.set_xlim(0, 6)
axins.set_ylim(0.5, 0.9)
plt.xticks(visible=False)
plt.yticks(visible=False)
mark_inset(ax1, axins, loc1=2, loc2=3, fc="none", ec="0.5")
axins.plot(v_MS, overlap_1, linewidth=5)
axins.plot(v_MS, overlap_2, linewidth=5, color='#87138A')
axins.plot(v_MS, overlap_3, linewidth=5, color='#088F8F')
axins.plot(v_MS, overlap_4, linewidth=5, color='#89CFF0')

axins2 = inset_axes(ax1, 2.5, 2.5, loc=1, bbox_to_anchor=(0.5, 0.55), bbox_transform=ax1.figure.transFigure)

axins2.set_xlim(0, 6)
axins2.set_ylim(-0.7, -0.2)
plt.xticks(visible=False)
plt.yticks(visible=False)
mark_inset(ax1, axins2, loc1=2, loc2=3, fc="none", ec="0.5")
axins2.plot(v_MS, overlap_1, linewidth=5)
axins2.plot(v_MS, overlap_2, linewidth=5, color='#87138A')
axins2.plot(v_MS, overlap_3, linewidth=5, color='#088F8F')
axins2.plot(v_MS, overlap_4, linewidth=5, color='#89CFF0')
'''

'''
#---Subplot 1------------
ax1 = plt.subplot(gs[0, 0])
ax1.tick_params(axis='x', direction='out', pad=15)
ax1.tick_params(axis='y', direction='out', pad=15)
plt.axhline(y=0, linestyle='--', linewidth=2, color='grey', dashes=(5, 5), alpha=0.6)
plt.axhline(y=1, linestyle='--', linewidth=2, color='grey', dashes=(5, 5), alpha=0.6)
plt.axhline(y=-1, linestyle='--', linewidth=2, color='grey', dashes=(5, 5), alpha=0.6)
plt.plot(MS, overlap, linewidth=10, label = 'm$^1$(k) para condicion inicial aleatoria', color='#40B5AD')
plt.xscale('log')
plt.xlabel('T (K)', fontsize=25, labelpad=20)
plt.ylabel('Solapamiento con $Salchicha$', fontsize=25, labelpad=15)
legend = plt.legend(loc=0, bbox_to_anchor=(0.95, 0.9), prop={'size': 25})
for legend_handle in legend.legendHandles:
    legend_handle._sizes = [scatter_legend_size]
#ax1.set_xlim(0, 4)
ax1.set_ylim(-1.2, 1.2)
plt.grid()
'''
'''
ax1 = plt.subplot(gs[0, 0])
ax1.tick_params(axis='x', direction='out', pad=15)
ax1.tick_params(axis='y', direction='out', pad=15)
plt.axhline(y=0, linestyle='--', linewidth=2, color='grey', dashes=(5, 5), alpha=0.6)
plt.axhline(y=1, linestyle='--', linewidth=2, color='grey', dashes=(5, 5), alpha=0.6)
plt.axhline(y=-1, linestyle='--', linewidth=2, color='grey', dashes=(5, 5), alpha=0.6)
plt.plot(T, T1, linewidth=5, label = 'm$^1$ condicion inicial aleatoria')
plt.plot(T, T2, linewidth=5, label = 'm$^1$ para def. 15%', color='#1F51FF')
plt.plot(T, T3, linewidth=5, label = 'm$^1$ para def. 30%', color='#088F8F')
plt.plot(T, T4, linewidth=5, label = 'm$^1$ para def. 40%', color='#87138A')
plt.plot(T, T5, linewidth=5, label = 'm$^1$ para def. 50%', color='#89CFF0')
plt.xscale('log')
plt.xlabel('T (K)', fontsize=25, labelpad=20)
plt.ylabel('Solapamiento con $Salchicha$', fontsize=25, labelpad=15)
legend = plt.legend(loc=0, prop={'size': 25})
for legend_handle in legend.legendHandles:
    legend_handle._sizes = [scatter_legend_size]
#ax1.set_xlim(0, 4)
ax1.set_ylim(-1.2, 1.2)
plt.grid()

axins = inset_axes(ax1, 3, 3, loc=1, bbox_to_anchor=(0.6, 0.8), bbox_transform=ax1.figure.transFigure)

axins.set_xlim(0.05, 0.11)
axins.set_ylim(-0.3, 0.3)
plt.xticks(visible=False)
plt.yticks(visible=False)
mark_inset(ax1, axins, loc1=1, loc2=4, fc="none", ec="0.5")
axins.plot(T, T1, linewidth=5)
axins.plot(T, T2, linewidth=5, color='#1F51FF')
axins.plot(T, T3, linewidth=5, color='#088F8F')
axins.plot(T, T4, linewidth=5, color='#87138A')
axins.plot(T, T5, linewidth=5, color='#89CFF0')
'''

'''
ax1 = plt.subplot(gs[0, 0])
ax1.tick_params(axis='x', direction='out', pad=15)
ax1.tick_params(axis='y', direction='out', pad=15)
plt.axhline(y=0, linestyle='--', linewidth=2, color='grey', dashes=(5, 5), alpha=0.6)
plt.axhline(y=1, linestyle='--', linewidth=2, color='grey', dashes=(5, 5), alpha=0.6)
plt.axhline(y=-1, linestyle='--', linewidth=2, color='grey', dashes=(5, 5), alpha=0.6)
plt.plot(T, T1, linewidth=5, label = 'm$^1$')
plt.plot(T, T2, linewidth=5, label = 'm$^2$', color='#1F51FF')
plt.plot(T, T3, linewidth=5, label = 'm$^3$', color='#088F8F')
plt.plot(T, T4, linewidth=5, label = 'm$^4$', color='#87138A')
#plt.plot(T, T5, linewidth=5, label = 'm$^1$ para def. 50%', color='#89CFF0')
plt.xscale('log')
plt.xlabel('T (K)', fontsize=25, labelpad=20)
plt.ylabel('Solapamiento de cada patrón', fontsize=25, labelpad=15)
legend = plt.legend(loc=1, prop={'size': 25})
for legend_handle in legend.legendHandles:
    legend_handle._sizes = [scatter_legend_size]
#ax1.set_xlim(0, 4)
ax1.set_ylim(-1.2, 1.2)
plt.grid()
'''
'''
axins = inset_axes(ax1, 2.8, 2.8, loc=1, bbox_to_anchor=(0.9, 0.55), bbox_transform=ax1.figure.transFigure)

axins.set_xlim(0.06, 0.18)
axins.set_ylim(-0.4, 0.5)
plt.xticks(visible=False)
plt.yticks(visible=False)
mark_inset(ax1, axins, loc1=2, loc2=3, fc="none", ec="0.5")
axins.plot(T, T1, linewidth=5)
axins.plot(T, T2, linewidth=5, color='#1F51FF')
axins.plot(T, T3, linewidth=5, color='#088F8F')
axins.plot(T, T4, linewidth=5, color='#87138A')
#axins.plot(T, T5, linewidth=5, color='#89CFF0')
'''
'''
valores = [1]*40
y = [0]*40
x = [0]*40
ideal = [0]*50
x_ideal = [0]*50


for i in range(1, 41):
    y[i-1] = 1

for i in range(1, 41):
    x[i-1] = i
'''
'''   
for i in range(14, 41):
    y[i-1] = 0
    
y[15] = y[18] = y[21] = y[22] = y[29] = 1
'''
'''
for i in range(1, 41):
    valores[i-1] = y[i-1]/x[i-1]
'''
'''   
for i in range(1, 51):
    ideal[i-1] = 1/i
    x_ideal[i-1] = i
'''
'''
ax1 = plt.subplot(gs[0, 0])
ax1.tick_params(axis='x', direction='out', pad=15)
ax1.tick_params(axis='y', direction='out', pad=15)
plt.plot(x, valores, linewidth=7, label = 'P$_R$ / N$_P$', color="green")
plt.xlabel('Área de la red (N$^2$)', fontsize=25, labelpad=20)
plt.ylabel('Nº patrones max. para recuerdo total (P$_c$)', fontsize=25, labelpad=15)
legend = plt.legend(loc=1, prop={'size': 25})
plt.xscale('log')
for legend_handle in legend.legendHandles:
    legend_handle._sizes = [scatter_legend_size]
#plt.xticks(np.arange(1, 100000, 5))
ax1.set_xlim(1, 100000)
ax1.set_ylim(-0.1, 1.1)
plt.grid()
'''
ax1 = plt.subplot(gs[0, 0])
ax1.tick_params(axis='x', direction='out', pad=15)
ax1.tick_params(axis='y', direction='out', pad=15)
plt.axhline(y=0, linestyle='--', linewidth=2, color='grey', dashes=(5, 5), alpha=0.6)
plt.axhline(y=1, linestyle='--', linewidth=2, color='grey', dashes=(5, 5), alpha=0.6)
plt.axhline(y=-1, linestyle='--', linewidth=2, color='grey', dashes=(5, 5), alpha=0.6)
plt.plot(tv, s10, linewidth=7, label = 'sesgo = -0.5', color="#f600ff")
plt.plot(tv, s9, linewidth=7, label = 'sesgo = -0.4', color="#ad32d7")
plt.plot(tv, s8, linewidth=7, label = 'sesgo = -0.3', color="#6e35a8")
plt.plot(tv, s7, linewidth=7, label = 'sesgo = -0.2', color="#3c2b74")
plt.plot(tv, s6, linewidth=7, label = 'sesgo = -0.1', color="#1c1b3e")
plt.plot(tv, s0, linewidth=7, label = 'sesgo = 0.0', color="black")
plt.plot(tv, s1, linewidth=7, label = 'sesgo = 0.1', color="#1f2635")
plt.plot(tv, s2, linewidth=7, label = 'sesgo = 0.2', color="#314865")
plt.plot(tv, s3, linewidth=7, label = 'sesgo = 0.3', color="#3a6e98")
plt.plot(tv, s4, linewidth=7, label = 'sesgo = 0.4', color="#3599cc")
plt.plot(tv, s5, linewidth=7, label = 'sesgo = 0.5', color="#00c6ff")
plt.xlabel('Pasos de Montecarlo (PM)', fontsize=25, labelpad=20)
plt.ylabel('Solapamiento con $Salchicha$', fontsize=25, labelpad=15)
legend = plt.legend(loc=7, prop={'size': 25})
for legend_handle in legend.legendHandles:
    legend_handle._sizes = [scatter_legend_size]
#plt.xticks(np.arange(1, 100000, 5))
ax1.set_xlim(0, 31)
ax1.set_ylim(-1.1, 1.1)
plt.grid()

axins2 = inset_axes(ax1, 3, 5, loc=1, bbox_to_anchor=(0.75, 0.55), bbox_transform=ax1.figure.transFigure)

axins2.set_xlim(0.95, 1.05)
axins2.set_ylim(-0.643, -0.604)
plt.xticks(visible=False)
plt.yticks(visible=False)
mark_inset(ax1, axins2, loc1=2, loc2=3, fc="none", ec="0.5")
axins2.axhline(y=0, linestyle='--', linewidth=2, color='grey', dashes=(5, 5), alpha=0.6)
axins2.axhline(y=1, linestyle='--', linewidth=2, color='grey', dashes=(5, 5), alpha=0.6)
axins2.axhline(y=-1, linestyle='--', linewidth=2, color='grey', dashes=(5, 5), alpha=0.6)
axins2.plot(tv, s10, linewidth=7, label = 'sesgo = -0.5', color="#f600ff")
axins2.plot(tv, s9, linewidth=7, label = 'sesgo = -0.4', color="#ad32d7")
axins2.plot(tv, s8, linewidth=7, label = 'sesgo = -0.3', color="#6e35a8")
axins2.plot(tv, s7, linewidth=7, label = 'sesgo = -0.2', color="#3c2b74")
axins2.plot(tv, s6, linewidth=7, label = 'sesgo = -0.1', color="#1c1b3e")
axins2.plot(tv, s0, linewidth=7, label = 'sesgo =  0.0', color="black")
axins2.plot(tv, s1, linewidth=7, label = 'sesgo =  0.1', color="#1f2635")
axins2.plot(tv, s2, linewidth=7, label = 'sesgo =  0.2', color="#314865")
axins2.plot(tv, s3, linewidth=7, label = 'sesgo =  0.3', color="#3a6e98")
axins2.plot(tv, s4, linewidth=7, label = 'sesgo =  0.4', color="#3599cc")
axins2.plot(tv, s5, linewidth=7, label = 'sesgo =  0.5', color="#00c6ff")

axins = inset_axes(ax1, 3, 6, loc=1, bbox_to_anchor=(0.5, 0.9), bbox_transform=ax1.figure.transFigure)

axins.set_xlim(0.99, 1.04)
axins.set_ylim(0.607, 0.645)
plt.xticks(visible=False)
plt.yticks(visible=False)
mark_inset(ax1, axins, loc1=2, loc2=3, fc="none", ec="0.5")
axins.axhline(y=0, linestyle='--', linewidth=2, color='grey', dashes=(5, 5), alpha=0.6)
axins.axhline(y=1, linestyle='--', linewidth=2, color='grey', dashes=(5, 5), alpha=0.6)
axins.axhline(y=-1, linestyle='--', linewidth=2, color='grey', dashes=(5, 5), alpha=0.6)
axins.plot(tv, s10, linewidth=7, label = 'sesgo = -0.5', color="#f600ff")
axins.plot(tv, s9, linewidth=7, label = 'sesgo = -0.4', color="#ad32d7")
axins.plot(tv, s8, linewidth=7, label = 'sesgo = -0.3', color="#6e35a8")
axins.plot(tv, s7, linewidth=7, label = 'sesgo = -0.2', color="#3c2b74")
axins.plot(tv, s6, linewidth=7, label = 'sesgo = -0.1', color="#1c1b3e")
axins.plot(tv, s0, linewidth=7, label = 'sesgo = 0.0', color="black")
axins.plot(tv, s1, linewidth=7, label = 'sesgo = 0.1', color="#1f2635")
axins.plot(tv, s2, linewidth=7, label = 'sesgo = 0.2', color="#314865")
axins.plot(tv, s3, linewidth=7, label = 'sesgo = 0.3', color="#3a6e98")
axins.plot(tv, s4, linewidth=7, label = 'sesgo = 0.4', color="#3599cc")
axins.plot(tv, s5, linewidth=7, label = 'sesgo = 0.5', color="#00c6ff")


