# -*- coding: utf-8 -*-
"""
Created on Sun Mar 12 18:42:06 2023

@author: adria
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
print(os.getcwd())

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
distance_2p = []

for line in lines:
    values = line.strip().split(' ')
    t.append(float(values[0]))
    potential.append(float(values[1]))
    kinetic.append(float(values[2]))
    total.append(float(values[3]))
    velocity.append(float(values[4]))
    momentum.append(float(values[5]))
    fluctuations.append(float(values[6]))
    distance_2p.append(float(values[7]))     

#-------------------------------------------------------------
# Calculate the average value of velocity during a period of time
start_time = 20  # Specify the start time of the period
end_time = 50    # Specify the end time of the period

time_indices = np.where((np.array(t) >= start_time) & (np.array(t) <= end_time))
velocity_period = np.array(velocity)[time_indices]
average_velocity = np.mean(velocity_period)

# Print the average value
print("Average sum of square velocities during the period from", start_time, "to", end_time, "is:", round(average_velocity, 2), 'm2/s2')
Temperature = 0.5*average_velocity
print("The temperature of the system is", round(Temperature, 2), "K")

#-------------------------------------------------------------
# Plot evolution of temperature for each 0.5 second

start_time = 0  # Specify the start time of the period
end_time = 60    # Specify the end time of the period

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
print("With the Ideal gases equation, the pressure should be", round(Temperature / 100, 3), "units")

#-------------------------------------------------------------
# Study of velocities

v_distribution = []
v_module = []
vx = []
vy = []

v_new = np.linspace(0, 13, num=500)
vx_new = np.linspace(-15, 15, num=1000)

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

fig, ax = plt.subplots(2, 3, gridspec_kw={'height_ratios':[1, 2]}, figsize=(35,55))

plt.style.use('default')

plt.suptitle('Parameters', fontsize=30)
plt.rc('xtick', labelsize=25)
plt.rc('ytick', labelsize=25)

gs = gridspec.GridSpec(2, 3)


#---Subplot 1------------
ax1 = plt.subplot(gs[0, 0])
plt.scatter(t, potential, 5, label='Potential Energy')
plt.scatter(t, kinetic, 5, label='Kinetic Energy')
plt.scatter(t, total, 5, label='Total Energy')
plt.xlabel('Time (s)', fontsize=30, labelpad=20)
plt.ylabel('Energy of the System', fontsize=30, labelpad=20)
plt.title('Energy Evolution', fontsize=30, pad=20)  
plt.legend(loc=1, prop={'size': 20})
#ax1.set_ylim(-4, 4)
plt.grid()

#---Subplot 2------------

ax2 = plt.subplot(gs[1, 0])
#plt.hist(v_hist, bins= 5, range=(0, 2), label='|V|', color='green', linewidth=2, edgecolor = "white", width=0.2)
plt.bar(v_distribution, v_module, label='|V|', color='green', linewidth=2, edgecolor = "white", width=0.2)
plt.scatter(v_new, PDF, label='Theroretical Maxwell Distribution', color='red')
plt.xlabel('Velocities (m/s)', fontsize=30, labelpad=20)
plt.ylabel('Probabilites for each velocity', fontsize=30, labelpad=20)
plt.title('Maxwell distribution for |V|', fontsize=30, pad=20)  
plt.legend(loc=1, prop={'size': 15})
ax2.set_xlim(-1, 13)
#ax2.set_ylim(-0.01, 0.25)
plt.grid()

#---Subplot 3------------
ax3 = plt.subplot(gs[1, 1])
plt.bar(v_distribution, vx, label='Vx', linewidth=2, edgecolor = "white", width=0.2)
plt.scatter(vx_new, PDF_vx, label='Theroretical Gaussian Distribution', color='red')
plt.xlabel('Velocities (m/s)', fontsize=30, labelpad=20)
plt.ylabel('Probabilites for each velocity', fontsize=30, labelpad=20)
plt.title('Gaussian distribution for Vx', fontsize=30, pad=20)  
plt.legend(loc=1, prop={'size': 15})
ax3.set_xlim(-7.5, 13)
plt.grid()

#---Subplot 4------------
ax4 = plt.subplot(gs[1, 2])
plt.bar(v_distribution, vy, label='Vy', color = 'orange', linewidth=0.4, edgecolor = "white", width=0.2)
plt.scatter(vx_new, PDF_vy, label='Theroretical Gaussian Distribution', color='red')
plt.xlabel('Velocities (m/s)', fontsize=30, labelpad=20)
plt.ylabel('Probabilites for each velocity', fontsize=30, labelpad=20)
plt.title('Gaussian distribution for Vy', fontsize=30, pad=20)  
plt.legend(loc=1, prop={'size': 15})
ax4.set_xlim(-7.5, 13)
plt.grid()

#---Subplot 5------------
ax5 = plt.subplot(gs[0, 1])
plt.scatter(t, fluctuations, 2, label='Fluctuation', color='black')
plt.scatter(t, distance_2p, 2, label='Distance between 2 particles', color='purple')
plt.plot(time_indices_t, 0.5*average_velocities, linewidth = 3, label='Temperature evolution', color='orange')
plt.xlabel('Time (s)', fontsize=30, labelpad=20)
plt.ylabel('Distance', fontsize=30, labelpad=20)
plt.title('Distance from initial position', fontsize=30, pad=20)  
plt.legend(loc=0, prop={'size': 20})
plt.grid()

#---Subplot 6------------
ax5 = plt.subplot(gs[0, 2])
plt.bar(g_dist, g_values, label='g(r)', color='brown', linewidth=1, edgecolor = "white", width=0.05)
plt.xlabel('r (m)', fontsize=30, labelpad=20)
plt.ylabel('g(r)', fontsize=30, labelpad=20)
plt.title('g(r) function', fontsize=30, pad=20)  
plt.legend(loc=0, prop={'size': 20})
plt.grid()

plt.subplots_adjust(
top=0.91,
bottom=0.11,
left=0.125,
right=0.9,
hspace=0.305,
wspace=0.2
)

plt.show()












