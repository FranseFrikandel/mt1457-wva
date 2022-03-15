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

def X_func(t):
    global tmax
    t_control = np.array([0, 0.1*tmax, 0.2*tmax, 0.5*tmax,
                         0.6*tmax, 0.7*tmax, tmax])
    X_parms = np.array([0.85, 0.85, 0.3, 0.4, 0.4, 1, 1])

    for i in range(len(t_control)):
        if t_control[i] < t:
            return X_parms[i]
    return X_parms[-1]

def Y_func(t):
    global tmax
    t_control = np.array([0, 0.1*tmax, 0.2*tmax, 0.5*tmax,
                         0.6*tmax, 0.7*tmax, tmax])
    Y_parms = np.array([1, 1, 1, 1, 1, 1, 1])

    for i in range(len(t_control)):
        if t_control[i] < t:
            return Y_parms[i]
    return Y_parms[-1]

# water properties
rho_sw = 1025           # density of seawater [kg/m3]
viscositeit = 1.2872E-06

# ship data
m_ship = 358000         # ship mass [kg]
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

def resistance(speed, Y):
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
    weerstand = Y*C_Ts_cur*0.5*rho_sw*(speed**2)
    return weerstand

def engine(n_e, X):
    n_cyl = 6                   # number of cylinders [-]
    k_es = 2 
    m_f_nom = 1.314762      # nominal fuel injection [g]
    eta_td = 0.52
    LHV = 42700
    n_nom = 900/60

    Q_f = X * m_f_nom * LHV
    Q_loss_cooling = 1908.8 + 7635.2 * X
    W_loss_mech = 711.1 + 1659.3 * (n_e/n_nom)
    W_i = (Q_f - Q_loss_cooling) * eta_td
    W_e = W_i - W_loss_mech
    P_b = W_e * n_e * n_cyl / k_es
    M_b = P_b / (2 * np.pi * n_e)

    return M_b, P_b

def gearbox(n_e, M_e):
    eta_TRM = 0.9500        # transmission efficiency [-]
    i_gb = 4.2100           # gearbox ratio [-]
    I_tot = 200             # total mass of inertia of propulsion system [kg*m^2]
    return n_e/i_gb, M_e*i_gb*eta_TRM

def propeller(n_p, v_a):
    pass

simulation_length = tmax/dt + 1
v_s = np.zeros(simulation_length)
n_e = np.zeros(simulation_length)
a_s = np.zeros(simulation_length)
X = np.zeros(simulation_length)
Y = np.zeros(simulation_length)

# Initial values
v_s[0] = 0
n_e[0] = 200/60

for i, t in enumerate(time):
    # if i == 0:
    #     continue

    X[i] = X_func(t)
    Y[i] = Y_func(t)

    M_b[i], P_b[i] = engine(n_e[i], X[i])

    
    v_s[i+1] = v_s[i] + a[i] * dt

    if not (0.2 < X(t) < 1):
        print("Fuel rack buiten toegestaan bereik.")
    
    if n_e[i] > 900/60:
        print("Maximaal toerental overschreden.")