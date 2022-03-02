# -*- coding: utf-8 -*-
"""
Een volledige re-write van viskotter simulatie.
"""
import numpy as np
import matplotlib.pyplot as plt

# Simulation control variables
tmax = 3600
dt = 1
time = np.linspace(0, tmax, 1 + round((tmax)/dt))

in_p = 3.2830  # initial rpm
v_s0 = 6.5

def X(t):
    global tmax
    t_control = np.array([0, 0.1*tmax, 0.2*tmax, 0.5*tmax,
                         0.6*tmax, 0.7*tmax, tmax])
    X_parms = np.array([0.85, 0.85, 0.3, 0.4, 0.4, 1, 1])

    for i in range(len(t_control)):
        if t_control[i] < t:
            return X_parms[i]
    return X_parms[-1]

def Y(t):
    global tmax
    t_control = np.array([0, 0.1*tmax, 0.2*tmax, 0.5*tmax,
                         0.6*tmax, 0.7*tmax, tmax])
    Y_parms = np.array([1, 1, 1, 1, 1, 1, 1])

    for i in range(len(t_control)):
        if t_control[i] < t:
            return Y_parms[i]
    return Y_parms[-1]

# fuel properties
LHV = 42700             # Lower Heating Value [kJ/kg]

# water properties
rho_sw = 1025           # density of seawater [kg/m3]
viscositeit = 1.2872E-06

# ship data
m_ship = 358000         # ship mass [kg]
c1 = 1500               # resistance coefficient c1 in R = c1*vs^2
v_s0 = 6.5430           # ship design speed [m/s]
t = 0.1600              # thrust deduction factor[-]
w = 0.2000              # wake factor [-]
l = 31.5                # length of ship [m]
eta_h = (1-t)/(1-w)

# propellor data
D_p = 3                 # diameter of propellor [m]
K_T_a = -0.3821         # factor a in K_T = a*J + b [-]
K_T_b = 0.2885          # factor b in K_T = a*J + b [-]
K_T_factor = np.array((-0.3821, 0.2885))
K_Q_a = -0.03346        # factor a in K_Q = a*J + b [-]
K_Q_b = 0.0308          # factor b in K_Q = a*J + b [-]
eta_R = 1.0100          # relative rotative efficiency [-]

# engine data
m_f_nom = 1.314762      # nominal fuel injection [g]
eta_e = 0.3800          # nominal engine efficiency [-]
i = 6                   # number of cylinders [-]
k_es = 2                # k-factor for engines based on nr.of strokes per cycle
P_b = np.zeros(tmax)    # engine power [kW]
P_b[0] = 960            # Nominal engine power [kW]
M_b = np.zeros(tmax)    # engine torque [Nm]
M_b[0] = P_b[0]*1000/2/np.pi/(900/60)  # ([P_b*1000/2/np.pi/n_eng_nom])

# gearbox data
eta_TRM = 0.9500        # transmission efficiency [-]
i_gb = 4.2100           # gearbox ratio [-]
I_tot = 200             # total mass of inertia of propulsion system [kg*m^2]

def Resistance(speed, Y):
    global C_Ws, k, l_s, viscositeit, rho_sw
    C_Ws = np.array([[0. , 0.43588989, 0.87177979, 1.74355958, 2.61533937,
       3.48711915, 4.35889894, 5.23067873, 5.66656863, 6.10245852,
       6.53834842, 6.97423831], [0, -0.00096818, -0.00128629, -0.00125762,  0.00016455,
       0.00061338,  0.00130253,  0.00298342,  0.00390739,  0.00489266,
       0.00713291,  0.01095216]])
    k = 0.23261693620133994
    l_s = 33.5

    C_Ws_cur = np.interp(speed, C_Ws[0], C_Ws[1])
    Re_s_cur = (speed*l_s)/viscositeit
    C_Fs_cur = 0.075 / (((np.log10(Re_s_cur)-2))**2)
    C_Ts_cur = (1+k)*C_Fs_cur + C_Ws_cur
    weerstand = Y*C_Ts_cur*0.5*rho_sw*(snelheid_schip**2)*snelheid_schip
    return weerstand

def engine():

simulation_length = tmax/dt +
v_s = np.zeros(simulation_length)