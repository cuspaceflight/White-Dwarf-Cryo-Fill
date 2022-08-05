import numpy as np
import scipy.optimize
from CoolProp.CoolProp import PropsSI

# Valve properties
D_fill = 8e-3
D_vent = 1.5e-3

# Tank dimensions
L = 100e-2                  # m
r_o = 15e-2                 # m
t = 14e-3                   # m

# Tank properties
k_tank = 236                # W/m/K
#reflectivity_tank = 0.86    # http://rmico.com/bare-aluminum#:~:text=Metal%20Mirror%20%C2%BB%20Bare%20Aluminum,slight%20scattering%20throughout%20the%20spectrum.

# Ambient temperatures
T_air = 273.15 + 30         # K
p_air = 1.0125e5            # Pa

# Supply tank - saturated
T_supply = T_air
p_v_supply = PropsSI('P','T', T_supply,'Q', 0, 'NITROUSOXIDE')   
p_supply = p_v_supply

def get_mdots(p_run):
    # Use Dyer to predict fill mass flow rate
    def mdot_dyer(p2):
        Cd = 0.7
        A = np.pi * (D_fill/2)**2
        k = np.sqrt((p_supply - p2) / (p_v_supply - p2))
        s1 = PropsSI('S','T', T_supply,'Q', 0, 'NITROUSOXIDE')   
        h1 = PropsSI('HMASS','T', T_supply,'Q', 0, 'NITROUSOXIDE')   
        rho_1 = PropsSI('DMASS','T', T_supply,'Q', 0, 'NITROUSOXIDE')
        rho_2 = PropsSI('DMASS','S', s1,'P', p2, 'NITROUSOXIDE')  
        h2 = PropsSI('HMASS','S', s1,'P', p2, 'NITROUSOXIDE')  
        mdot_hem = Cd * A * rho_2 * np.sqrt(2 * (h1 - h2))
        mdot_spi = Cd * A * np.sqrt(2 * rho_1 * (p_supply - p2))
        return (k / (1 + k) *  mdot_spi + 1 / (1 + k) * mdot_hem)

    mdot_crit = max(mdot_dyer(np.linspace(0, p_supply, 100)))
    mdot_fill = min(mdot_dyer(p_run), mdot_crit)

    # Use SPI to predict vent mass flow rate
    Cd = 0.7
    A = np.pi * (D_vent/2)**2 
    rho_1 = PropsSI('DMASS','P', p_run,'Q', 1, 'NITROUSOXIDE')
    mdot_vent = Cd * A * np.sqrt(2 * rho_1 * (p_run - p_air))

    return mdot_vent, mdot_fill

def mdot_error(p_run):
    mdot_vent, mdot_fill = get_mdots(p_run)
    return mdot_vent - mdot_fill

p_run = scipy.optimize.root_scalar(mdot_error, bracket = [p_air + 100, p_supply - 100]).root
mdot_vent, mdot_fill = get_mdots(p_run)

# Calculate heat loss through venting
dh_boiloff = PropsSI('HMASS','P', p_run,'Q', 1, 'NITROUSOXIDE') - PropsSI('HMASS','P', p_run,'Q', 0, 'NITROUSOXIDE')   
Q_vent = mdot_vent * dh_boiloff

# Order of magnitude estimate for heat transfer through walls
T_run = PropsSI('T','P', p_run,'Q', 0, 'NITROUSOXIDE')   
#r_i = r_o - t
#R_walls = np.log(r_o / r_i) / (2 * np.pi * k_tank * L)          # Thermal resistance for radial conduction
#Q_walls = (T_air - T_run) / R_walls
#T_bulk_eqlm = T_air - Q_vent * R_walls

# Time to chill
m_desired_n2o = 5           # Desired mass of N2O in filled tank, kg
T_desired_n2o = 273.15 - 10 # Desired nitrous temperature

cp_n2o_l = PropsSI('CPMASS','P', p_run,'Q', 0, 'NITROUSOXIDE')   
dT_dt = (Q_vent) / (m_desired_n2o * cp_n2o_l)       # Q = m cp dT/dt, 
t_to_chill = (T_supply - T_desired_n2o) / dT_dt     # Time to chill nitrous only (ignoring tank thermal mass)

# Results
print(f"""
Boil-off rate = {mdot_vent*1000:.03} g/s
Heat loss by vent = {Q_vent/1e3:.03} kW
Run pressure = {p_run/1e5:.04} bar
Supply pressure = {p_supply/1e5:.04} bar
dT/dt = {dT_dt:.03} K/s
Time to chill = {t_to_chill/60:.03} min
""")
