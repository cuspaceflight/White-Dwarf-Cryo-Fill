import numpy as np
from CoolProp.CoolProp import PropsSI

# Tank dimensions
L = 100e-2                  # m
r_o = 15e-2                 # m
t = 14e-3                   # m


# Tank properties
k_tank = 236                # W/m/K
reflectivity_tank = 0.86    # http://rmico.com/bare-aluminum#:~:text=Metal%20Mirror%20%C2%BB%20Bare%20Aluminum,slight%20scattering%20throughout%20the%20spectrum.

# Ambient temperatures
T_air = 273.15 + 30         # K
p_air = 1.0125e5            # Pa

# Flow system and nitrous properties.
T_bulk_n2o = T_air          # K
m_desired_n2o = 5           # Desired mass of N2O in filled tank, kg
T_desired_n2o = 273.15 - 10 # Desired nitrous temperature
#mdot_fill_n2o = 0.1        # kg/s

T_boil_n2o = PropsSI('T','P', p_air,'Q', 0, 'NITROUSOXIDE')                                                             # K
h_boil_n2o = PropsSI('HMASS','P', p_air,'Q', 1, 'NITROUSOXIDE') - PropsSI('HMASS','P', p_air,'Q', 0, 'NITROUSOXIDE')    # J/kg
k_n2o = 103.0e-3                                                                                                        # W/m/K, http://edge.rit.edu/edge/P07106/public/Nox.pdf
cp_n2o_l = PropsSI('CPMASS','T', T_bulk_n2o,'Q', 0, 'NITROUSOXIDE')                                                     # J/kg/K

# Order of magnitude of internal conduction for boiling
Q_internal = k_n2o * (np.pi * r_o**2) * (T_bulk_n2o - T_boil_n2o)/L      # Assume axial conduction from bottom of tank to top
mdot_boiloff = Q_internal / h_boil_n2o

# Order of magnitude estimate for heat transfer through walls
r_i = r_o - t
R_walls = np.log(r_o / r_i) / (2 * np.pi * k_tank * L)
Q_walls = (T_air - T_bulk_n2o) / R_walls
T_bulk_eqlm = T_air - Q_internal * R_walls

dT_dt = (Q_internal - Q_walls) / (m_desired_n2o * cp_n2o_l)    # Q = m cp dT/dt, 
t_to_chill = (T_bulk_n2o - m_desired_n2o) / dT_dt

# Results
print(f"""
Boil-off rate = {mdot_boiloff} kg/s
dT/dt = {dT_dt} K/s
Time to chill = {t_to_chill/3600} hours
Equilibrium N2O temperature = {T_bulk_eqlm - 273.15} deg C
""")
