import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sp


earth_radius = 6371                   # km
altitude = 500                        # km

G = 6.6743*(10**-11)               # Gravitational constant in m^3/ kgs^2
g = 9.81                           # gravity in m/s
r = (earth_radius + altitude)      # km
M = 5.9722 * (10**24)              # mass of earth in kg

earth_tilt = math.radians(23.5)       # rad
arg_peri = math.radians(73)           # argument of perigee (rad)
mean_anomaly = math.radians(271)      # Mean = True in a circular orbit (rad)

raan = math.radians(42)               # rad
inclination = math.radians(97)        # rad


orbit = 90                        	# Length of one orbit
orbit_num = 16*365              	# Number of orbit cycles
spacing = (orbit * orbit_num) + 1 

# total time of n amount of orbits
t_min = np.linspace(0, orbit * orbit_num, spacing)

yr2min = 365 * 24 * 60
lamb = 2 * math.pi * (t_min) / yr2min     # math.radians(2 * math.pi) * (time since march equinox (march 19, 2028)) / 1 year
# ***** this accounts for the earth changing position from the sun over time!!!!!

# argument of latitude = mean anomaly + argument of perigee (rad)
arg_lat = arg_peri + mean_anomaly

# arg_lat changes throughout the orbit as denoted by u
u = (arg_lat + (2 * math.pi * t_min / orbit)) % (2 * math.pi)


# Power
gomspace = 1.2  						   # W per cell
mspa41 = gomspace * 4                      # 4 cell 2u
mspa51 = gomspace * 5					   # 5 cell 2u
msp12 = gomspace * 12					   # 12 cell 6u deployable


# cos(beta) value
cos_beta = ((-np.sin(lamb) * math.cos(earth_tilt) * math.sin(raan)) - (np.cos(lamb) * math.cos(raan))) ** 2
cos_beta += (((np.sin(lamb) * math.cos(earth_tilt) * math.cos(raan)) - (np.cos(lamb) * math.sin(raan))) * math.cos(inclination) + (np.sin(lamb) * math.sin(earth_tilt) * math.sin(inclination)))** 2
cos_beta = np.sqrt(cos_beta)


# sin(beta) value
sin_beta = np.sqrt(np.maximum(0, 1 - (cos_beta ** 2)))


# eclipse angle
ec_ang = np.zeros_like(t_min)

for i in range(len(t_min)) :
    if (r * sin_beta[i] < earth_radius):
        eclipse_angle = (r**2) - (earth_radius**2)
        eclipse_angle /= (earth_radius**2) - ((r**2) * (sin_beta[i]**2))
        eclipse_angle = np.arctan(math.sqrt(eclipse_angle))
        ec_ang[i] = eclipse_angle

    else:
        ec_ang[i] = math.pi / 2


# theta angle
theta_num = (-np.sin(lamb) * math.sin(earth_tilt) * math.sin(inclination))
theta_num += (((np.cos(lamb) * math.sin(raan)) - (np.sin(lamb) * math.cos(earth_tilt) * math.cos(raan))) * math.cos(inclination))
theta_den = ((np.sin(lamb) * math.cos(earth_tilt) * math.sin(raan)) + (np.cos(lamb) * math.cos(raan)))
theta = np.arctan2(theta_num, theta_den)
theta += u + (math.pi / 2)
theta = theta % (2*np.pi)


# power
power = np.zeros_like(t_min)
for j in range(len(t_min)) :
    
    # Theta less than pi/2
    if ((theta[j] < (math.pi / 2)) and (theta[j] >= 0)):
        # y+ face
        power[j] += mspa51 * np.sin(theta[j]) * cos_beta[j] * 1 
        power[j] += msp12 * np.sin(theta[j]) * cos_beta[j] * 1

    # Theta between pi/2 and pi
    elif ((theta[j] >= (math.pi / 2)) and (theta[j] < math.pi)):
        # y+ face
        power[j] += mspa51 * np.sin(theta[j]) * cos_beta[j] * 1 
        power[j] += msp12 * np.sin(theta[j]) * cos_beta[j] * 1

    # Theta between pi and pi + eclipse angle
    elif ((theta[j] >= (math.pi)) and (theta[j] < (math.pi + ec_ang[j]))):
        power[j] += 0

    # Theta between pi + eclipse angle and 2pi - eclipse angle
    elif ((theta[j] >= (math.pi + ec_ang[j])) and (theta[j] < (2 * math.pi) - ec_ang[j])):
        power[j] = 0    # No power generated

    # Theta between 2pi - eclipse angle and 2pi
    elif ((theta[j] >= (2 * math.pi) - ec_ang[j]) and (theta[j] < (2 * math.pi))): 
        power[j] += 0
        
print(f"Max power: {np.max(power):.4f}")


''' watt hour calcs start'''
# Watt-hour = power * time where it can equate to the area under the curve using the riemann sum as we have an instantaneous power per minute
# cant use an integral since there isn't a function and the peaks vary sinusoidally throughout the year
# currently the trapezoid version is used as it gives a more acurate value compared to the low and high version
def trap_sum(string_val, interval_s, interval_e):
    r = 0
    for i in range(interval_s-1, interval_e):
        r += ((1/2)*(power[i] + power[i+1])) * (1/60)
        
    short_line = power[interval_s-1:interval_e]
    if (string_val == "Minimum"):
        
        sect = 3.7180
        proper_index = [193440]
    else:
        sect = np.max(short_line)
        peak, point = sp.find_peaks(short_line, (sect))
        proper_index = peak+interval_s-1
    print(string_val, sect, proper_index)
    print(string_val, "r =", r, "Watt-hour\n")
    


trap_sum("First orbit", 1, 90)
trap_sum("First day", 1, 1440)
trap_sum("Maximum", 61388, 62828)
trap_sum("Minimum", 192720, 194160)
trap_sum("Middle 1", 151483, 152923)
trap_sum("Middle 2", 263211, 264651)
    
'''watt hour calcs end'''



# P as a function of t
# plt.plot(t_min, power)
# plt.ylim(0, np.max(power))
# plt.title("Power Simulation")
# plt.ylabel("Power (W)")
# plt.xlabel("Time (min)")
# plt.yticks(np.arange(0, np.max(power), 1))
# plt.grid(True)

plt.figure()
plt.plot(t_min[:90], power[:90])				
plt.title("one orbit first one")
plt.ylabel("Power (W)")
plt.xlabel("Time (min)")
plt.grid(True)
  

plt.figure()
plt.plot(t_min[263211:264651], power[263211:264651])				
plt.title("Power over few orbits") # Middle 2m near start
plt.ylabel("Power (W)")
plt.xlabel("Time (min)")
plt.grid(True)


plt.figure()
plt.plot(t_min[61388:62828], power[61388:62828])
plt.title("max Power over few orbits")
plt.ylabel("Power (W)")
plt.xlabel("Time (min)")
plt.grid(True)

# plt.figure()
# plt.plot(t_min[150000:151440], power[150000:151440])
# plt.title("Power over few orbits")
# plt.ylabel("Power (W)")
# plt.xlabel("Time (min)")
# plt.grid(True)

plt.figure()
plt.plot(t_min[192720:194160], power[192720:194160])
plt.title("min Power over few orbits")
plt.ylabel("Power (W)")
plt.xlabel("Time (min)")
plt.grid(True)

plt.figure()
plt.plot(t_min[151483:152923], power[151483:152923])
plt.title("Power over few orbits") #middle # 1 , falling 
plt.ylabel("Power (W)")
plt.xlabel("Time (min)")
plt.grid(True)

plt.figure()
plt.plot(t_min[193440:193530], power[193440:193530])
plt.title("one orbit, hitting min Power over few orbits")
plt.ylabel("Power (W)")
plt.xlabel("Time (min)")
plt.grid(True)

plt.show()