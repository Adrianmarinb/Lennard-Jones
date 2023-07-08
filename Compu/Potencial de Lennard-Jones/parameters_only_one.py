# -*- coding: utf-8 -*-
"""
Created on Sun Mar 12 18:42:06 2023

@author: adria
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import mark_inset, inset_axes
import os
print(os.getcwd())

# This Python script allows us to show all the different properties and parameters for the Lennard-Jones proyect.
# To show the desired property or parameter, the rest of the subplots are framed in /* */

#-------------------------------------------------------------
#We extract omega from .txt and remove points with no values


with open('parameters.txt') as f:
    lines = f.readlines()[0:]

t = []
potential = []
kinetic = []
total = []
velocity = []
momentum = []
fluctuations = []
distance_1p = []
distance_2p = []
distance_3p = []

for line in lines:
    values = line.strip().split(' ')
    t.append(float(values[0]))
    potential.append(float(values[1]))
    kinetic.append(float(values[2]))
    total.append(float(values[3]))
    velocity.append(float(values[4]))
    momentum.append(float(values[5]))
    fluctuations.append(float(values[6]))
    distance_1p.append(float(values[7]))
    distance_2p.append(float(values[8]))
    distance_3p.append(float(values[9]))
    
std_dev = np.std(total)
mean_e = np.mean(total)
print('Mean value of total energy: ', mean_e)
print('Standar deviation of total energy: ', std_dev)

#-------------------------------------------------------------
# Calculate the average value of velocity during a period of time

start_time = 5  # Specify the start time of the period
end_time = 30    # Specify the end time of the period

time_indices = np.where((np.array(t) >= start_time) & (np.array(t) <= end_time))
velocity_period = np.array(velocity)[time_indices]
average_velocity = np.mean(velocity_period)

# Print the average value
print("Average sum of square velocities during the period from", start_time, "to", end_time, "is:", round(average_velocity, 2), 'm2/s2')
Temperature = 0.5*average_velocity
print("The temperature of the system is", round(Temperature, 2), "K")

#-------------------------------------------------------------
# Comparison between temperatures

with open('comparison.tsv') as f:
    lines = f.readlines()[0:]

times = []
T1 = []
T2 = []
T3 = []
T4 = []

for line in lines:
    values = line.strip().split('\t')
    times.append(float(values[0]))
    T1.append(float(values[1]))
    T2.append(float(values[2]))
    T3.append(float(values[3]))
    T4.append(float(values[4]))

#-------------------------------------------------------------
# Plot evolution of temperature for each 0.5 second

start_time = 5  # Specify the start time of the period
end_time = 30    # Specify the end time of the period

# Assuming you have arrays t and velocity defined

interval = 1  # Specify the interval for calculating average velocity
time_indices_t = np.arange(start_time, end_time, interval)

average_velocities = []

for time in time_indices_t:
    time_range = np.where((np.array(t) >= time) & (np.array(t) < time + interval))
    velocity_period = np.array(velocity)[time_range]
    average_velocity = np.mean(velocity_period)
    average_velocities.append(average_velocity)

average_velocities = np.array(average_velocities)
temperatures = 0.5 * average_velocities
temp2 = temperatures[:-10]

stdt = np.std(temp2)
mean = np.mean(temp2)
print('Mean value of temperature: ', stdt)
print('Standar deviation of temperature: ', mean)

#-------------------------------------------------------------
# Study of the pressure of the system

def sum_array_values(arr):
    total = 0
    for value in arr:
        total += value
    return total

momentum_period = np.array(momentum)[time_indices]
Force = (1/(end_time - start_time)) * sum_array_values(momentum_period)

Area = 40

Pressure = Force / Area

Relation = Pressure / Temperature

print("The relationship Pressure/Temperature of the system is", round(Relation, 3), "K. Giving us a pressure of:", round(Pressure, 3), "Pa.")
print("With the Ideal gases equation, the pressure should be", round(Temperature * 16 / 100, 3), "units. This relationship would be", round(16 / 100, 3 ))

#-------------------------------------------------------------
# Study of velocities

v_distribution = []
v_module = []
vx = []
vy = []

v_new = np.linspace(0, 10, num=1000)
vx_new = np.linspace(-7, 7, num=1000)

def maxwell_distribution(v, T):
    exponent = -v**2 / (2 * T)
    return (v/T) * np.exp(exponent)

def gaussian_distribution(v, T):
    exponent = -v**2 / (2 * T)
    return np.sqrt(1/(2*np.pi*T)) * np.exp(exponent)

with open('v_statistics.txt') as f:
    lines = f.readlines()[0:]
  
for line in lines:
    values = line.strip().split(' ')
    v_distribution.append(float(values[0]))
    v_module.append(float(values[1]))
    vx.append(float(values[2]))
    vy.append(float(values[3]))
    
PDF = maxwell_distribution(np.array(v_new), Temperature)
PDF_vx = gaussian_distribution(np.array(vx_new), Temperature)
PDF_vy = gaussian_distribution(np.array(vx_new), Temperature)
    
#-------------------------------------------------------------
# Study of g_function

g_dist = []
g_values = []

with open('g_statistics.txt') as f:
    lines = f.readlines()[0:]
  
for line in lines:
    values = line.strip().split(' ')
    g_dist.append(float(values[0]))
    g_values.append(float(values[1]))

#-------------------------------------------------------------
# Plotting the data
    
#fig, ax = plt.subplots(1, 1, figsize=(20,8)) #Para energias y temperatura y maxwellV0
#fig, ax = plt.subplots(1, 1, figsize=(10,10))
fig, ax = plt.subplots(1, 1, figsize=(20,10))

plt.style.use('default')

plt.rc('xtick', labelsize=25)
plt.rc('ytick', labelsize=25)

gs = gridspec.GridSpec(1, 1)
scatter_legend_size = 60


#---Subplot 1------------
ax1 = plt.subplot(gs[0, 0])
ax1.tick_params(axis='x', direction='out', pad=15)
ax1.tick_params(axis='y', direction='out', pad=15)
plt.scatter(t, potential, s=5, label='Energía potencial media')
plt.scatter(t, kinetic, 5, label='Energía cinética media')
plt.scatter(t, total, 5, label='Energía total media')
#plt.plot(time_indices_t, 0.5*average_velocities, linewidth = 5, label='Evolución de T para intervalos de 1 s', color='orange')
plt.xlabel('t (s)', fontsize=25, labelpad=20)
plt.ylabel('E/ϵ', fontsize=25, labelpad=15)
#plt.title('Potencial de Lennard-Jones', fontsize=15, pad=10)  

plt.axvline(x=20, color='black', linestyle='--', lw=2, label= "Reescalado de velocidad (1.1)", alpha=0.3)
plt.axvline(x=25, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=30, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=35, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=40, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=45, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=50, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=55, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=60, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=65, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=70, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=75, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=80, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=85, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=90, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=95, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=100, color='black', linestyle='--', lw=2, alpha=0.3)

legend = ax1.legend(loc=9, prop={'size': 25})
for legend_handle in legend.legendHandles:
    legend_handle._sizes = [scatter_legend_size]
#ax1.set_xlim(-1, 11)
#ax1.set_ylim(-2.5,4 )
plt.xticks(np.arange(-0, 100.01, 5))
plt.grid()

axins = inset_axes(ax1, 2.8, 2.8, loc=1, bbox_to_anchor=(0.32, 0.9), bbox_transform=ax1.figure.transFigure)

axins.set_xlim(2, 6)
axins.set_ylim(-2.25, 1.35)
plt.xticks(visible=False)
plt.yticks(visible=False)
mark_inset(ax1, axins, loc1=2, loc2=4, fc="none", ec="0.5")
axins.scatter(t, potential, 5)
axins.scatter(t, kinetic, 5)
axins.scatter(t, total, 5)

'''
ax1 = plt.subplot(gs[0, 0])
#plt.hist(v_hist, bins= 5, range=(0, 2), label='|V|', color='green', linewidth=2, edgecolor = "white", width=0.2)
plt.bar(v_distribution, v_module, label='Distribución de velocidades registradas en t = 0.000 s', color='green', linewidth=2, edgecolor = "white", width=0.2)
plt.scatter(v_new, PDF, label='Distribución teórica de Maxwell para T', color='red')
plt.xlabel('v (r/σs)', fontsize=25, labelpad=20)
plt.ylabel('P(v)', fontsize=25, labelpad=20)
legend = ax1.legend(loc=1, prop={'size': 25})
for legend_handle in legend.legendHandles:
    legend_handle._sizes = [scatter_legend_size]
ax1.set_xlim(-0.2, 1.5)
#ax1.set_ylim(0, 0.16)
plt.xticks(np.arange(-0.2, 8, 0.5))
plt.grid()
'''
'''
ax1 = plt.subplot(gs[0, 0])
plt.grid()
#plt.hist(v_hist, bins= 5, range=(0, 2), label='|V|', color='green', linewidth=2, edgecolor = "white", width=0.2)
plt.bar(v_distribution, vx, label='Distribución registrada para v$_x$', linewidth=1, edgecolor = "white", width=0.2)
plt.scatter(vx_new, PDF_vx, label='Distribución teórica Gaussiana para T', color='red')
plt.xlabel('v$_x$ (r/σs)', fontsize=25, labelpad=20)
plt.ylabel('P(v$_x$)', fontsize=25, labelpad=20)
legend = ax1.legend(loc=1, prop={'size': 25})
for legend_handle in legend.legendHandles:
    legend_handle._sizes = [scatter_legend_size]
ax1.set_xlim(-4, 4)
#ax1.set_ylim(0, 26)
plt.xticks(np.arange(-7, 7, 1))
'''
'''
ax1 = plt.subplot(gs[0, 0])
plt.grid()
#plt.hist(v_hist, bins= 5, range=(0, 2), label='|V|', color='green', linewidth=2, edgecolor = "white", width=0.2)
plt.bar(v_distribution, vy, label='Distribución registrada para v$_y$', linewidth=1, edgecolor = "white", width=0.2, color='orange')
plt.scatter(vx_new, PDF_vy, label='Distribución teórica Gaussiana para T', color='red')
plt.xlabel('v$_y$ (r/σs)', fontsize=25, labelpad=20)
plt.ylabel('P(v$_y$)', fontsize=25, labelpad=20)
legend = ax1.legend(loc=1, prop={'size': 25})
for legend_handle in legend.legendHandles:
    legend_handle._sizes = [scatter_legend_size]
ax1.set_xlim(-4, 4)
#ax1.set_ylim(0, 26)
plt.xticks(np.arange(-7, 7, 1))
'''
'''
ax1 = plt.subplot(gs[0, 0])
ax1.tick_params(axis='x', direction='out', pad=15)
ax1.tick_params(axis='y', direction='out', pad=15)
plt.plot(times, T1, color='orange', linewidth=5, label='T para v = 1')
plt.plot(times, T2, times, T3, times, T4, linewidth=5)
plt.xlabel('t (s)', fontsize=25, labelpad=20)
plt.ylabel('E/ϵ', fontsize=25, labelpad=15)
legend = ax1.legend(['T para v = 1', 'T para v = 2', 'T para v = 3', 'T para v = 4'], loc=1, prop={'size': 25})
ax1.set_ylim(8, 18)
plt.grid()
'''

'''
ax1 = plt.subplot(gs[0, 0])
plt.scatter(t, fluctuations, 2, label='Fluctuación para la partícula 1', color='purple')
#plt.scatter(t, distance_2p, 2, label='Distance between 2 particles', color='purple')
plt.plot(time_indices_t, 0.5*average_velocities, linewidth = 3, label='Evolución de la temperatura T (K)', color='orange')
plt.xlabel('t (s)', fontsize=25, labelpad=20)
plt.ylabel('$<(r(t) - r(t=0))^2> (\sigma) $  ', fontsize=25, labelpad=20)

plt.axvline(x=20, color='black', linestyle='--', lw=2, label= 'Reescalado de velocidad (1.5)')
plt.axvline(x=30, color='black', linestyle='--', lw=2)
plt.axvline(x=35, color='black', linestyle='--', lw=2)
plt.axvline(x=45, color='black', linestyle='--', lw=2)

legend = ax1.legend(loc=0, prop={'size': 25})
for legend_handle in legend.legendHandles:
    legend_handle._sizes = [scatter_legend_size]
    
plt.xticks(np.arange(-0, 50.01, 5))
    
plt.grid()
'''
'''
interval_r = 0.5  # Specify the interval for calculating average velocity
time_indices_r = np.arange(start_time, end_time, interval_r)

r1 = []
r2 = []
r3 = []


for time in time_indices_r:
    time_range = np.where((np.array(t) >= time) & (np.array(t) < time + interval))
    period_d1 = np.array(distance_1p)[time_range]
    r1_t = np.mean(period_d1)
    r1.append(r1_t)
    
    period_d2 = np.array(distance_2p)[time_range]
    r2_t = np.mean(period_d2)
    r2.append(r2_t)
    
    period_d3 = np.array(distance_3p)[time_range]
    r3_t = np.mean(period_d3)
    r3.append(r3_t)



ax1 = plt.subplot(gs[0, 0])
#plt.plot(time_indices_r, r1, time_indices_r, r2, time_indices_r, r3, linewidth=5)
plt.plot(time_indices_r, r1, color='red', linewidth=5)
plt.plot(time_indices_r, r2, color='purple', linewidth=5)
plt.plot(time_indices_r, r3, color='pink', linewidth=5)
plt.plot(time_indices_t, 0.5*average_velocities, color='orange', linewidth=5)
plt.xlabel('t (s)', fontsize=25, labelpad=20)
plt.ylabel('<$(\Delta r_{ij}(t))^2> (\sigma)$', fontsize=25, labelpad=20)
plt.legend(loc=0, prop={'size': 25})

plt.axvline(x=20, color='black', linestyle='--', lw=2, label= "Reescalado de velocidad (1.1)", alpha=0.3)
plt.axvline(x=25, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=30, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=35, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=40, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=45, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=50, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=55, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=60, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=65, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=70, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=75, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=80, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=85, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=90, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=95, color='black', linestyle='--', lw=2, alpha=0.3)
plt.axvline(x=100, color='black', linestyle='--', lw=2, alpha=0.3)

legend = ax1.legend(['Separación cuadrática entre p1 y p2', 'Separación cuadrática entre p1 y p6', 'Separación cuadrática entre p1 y p8', 'Evolución de la temperatura T (K)'], loc=0, prop={'size': 25})
for legend_handle in legend.legendHandles:
    legend_handle._sizes = [scatter_legend_size]
    
plt.xticks(np.arange(-0, 100.01, 5))

plt.grid()
'''
'''
#---Subplot 6------------
ax5 = plt.subplot(gs[0, 0])
plt.bar(g_dist, g_values, label='Función de correlación de pares para p1 en estado gaseoso', color='brown', linewidth=1, edgecolor = "white", width=0.042)
plt.xlabel('r ($\sigma$)', fontsize=25, labelpad=20)
plt.ylabel('g(r)', fontsize=25, labelpad=20)
plt.legend(loc=0, prop={'size': 20})
plt.xticks(np.arange(1, 3.5, 0.25))
plt.grid()
'''









