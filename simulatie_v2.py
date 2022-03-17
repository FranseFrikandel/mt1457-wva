# -*- coding: utf-8 -*-
"""
Een volledige re-write van viskotter simulatie.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

# Simulation control variables
tmax = 36000
dt = 1
time = np.linspace(0, tmax, 1 + round((tmax)/dt))

def X_func(t):
    return 0.85
    # global tmax
    # t_control = np.array([0, 0.1*tmax, 0.2*tmax, 0.5*tmax,
    #                      0.6*tmax, 0.7*tmax, tmax])
    # X_parms = np.array([0.85, 0.85, 0.3, 0.4, 0.4, 1, 1])

    # for i in range(len(t_control)):
    #     if t_control[i] < t:
    #         return X_parms[i]
    # return X_parms[-1]

def Y_func(t):
    return 1
    # global tmax
    # t_control = np.array([0, 0.1*tmax, 0.2*tmax, 0.5*tmax,
    #                      0.6*tmax, 0.7*tmax, tmax])
    # Y_parms = np.array([1, 1, 1, 1, 1, 1, 1])

    # for i in range(len(t_control)):
    #     if t_control[i] < t:
    #         return Y_parms[i]
    # return Y_parms[-1]

def clutch_func(t, X, n_e):
    """
    Functie die bepaalt of de motor wel of niet moet worden ontkoppelt. Kan op tijdsbasis of basis van andere
    parameters.
    """
    return True
    if X < 0.21 and n_e < 300/60:
        return False
    return True

# water properties
rho_sw = 1025         # density of seawater [kg/m3]
viscositeit = 1.2872E-06

# gearbox data
i_gb = 4.2100           # gearbox ratio [-]

# ship data
m_ship = 358000       # ship mass [kg]
v_s0 = 6.5430         # ship design speed [m/s]
t = 0.1600            # thrust deduction factor[-]
w = 0.2000            # wake factor [-]
l = 31.5              # length of ship [m]
I_prop = 140           # total mass of inertia of propulsion system [kg*m^2]
I_eng = 60 / i_gb
eta_h = (1-t)/(1-w)   # Hull efficiency [-]
D_p = 2.53

def resistance(speed, Y):
    global viscositeit, rho_sw
    C_Ws = np.array([[0. , 0.43588989, 0.87177979, 1.74355958, 2.61533937,
       3.48711915, 4.35889894, 5.23067873, 5.66656863, 6.10245852,
       6.53834842, 6.97423831], 
       
       [0, -0.00096818, -0.00128629, -0.00125762,  0.00016455,
       0.00061338,  0.00130253,  0.00298342,  0.00390739,  0.00489266,
       0.00713291,  0.01095216]])
    k = 0.23261693620133994
    l_s = 33.5
    S_s = 272.6

    C_Ws_cur = np.interp(speed, C_Ws[0], C_Ws[1])
    Re_s_cur = (speed*l_s)/viscositeit
    C_Fs_cur = 0.075 / (((np.log10(Re_s_cur)-2))**2)
    C_Ts_cur = (1+k)*C_Fs_cur + C_Ws_cur
    weerstand = Y*C_Ts_cur*0.5*rho_sw*(speed**2)*S_s
    return weerstand

def engine(n_e, X):
    n_cyl = 6                   # number of cylinders [-]
    k_es = 2 
    m_f_nom = 1.314762      # nominal fuel injection [g]
    eta_td = 0.52
    LHV = 42700
    n_nom = 900/60

    flux_fuel = (X * m_f_nom * n_e) * n_cyl / k_es

    Q_f = X * m_f_nom * LHV
    Q_loss_cooling = 1908.8 + 7635.2 * X
    W_loss_mech = 711.1 + 1659.3 * (n_e/n_nom)
    W_i = (Q_f - Q_loss_cooling) * eta_td
    W_e = W_i - W_loss_mech
    P_b = W_e * n_e * n_cyl / k_es
    M_b = P_b / (2 * np.pi * n_e)

    return M_b, P_b, flux_fuel

def gearbox(n_e, M_e):
    eta_TRM = 0.9500        # transmission efficiency [-]
    i_gb = 4.2100           # gearbox ratio [-]
    n_p = n_e / i_gb
    M_trm = M_e*i_gb*eta_TRM
    return n_p, M_trm

def propeller(n_p, J):
    global rho_sw
    D_p = 2.53              # diameter of propellor [m]
    eta_R = 1.0100          # relative rotative efficiency [-]

    K_T = - (0.164 * J*J) - (0.257 * J) + 0.276
    K_Q = - (0.0187 * J*J) - (0.0193 * J) + 0.0295
    F_prop = K_T * n_p ** 2 * rho_sw * D_p ** 4
    M_prop = (K_Q * n_p ** 2 * rho_sw * D_p ** 5) / eta_R

    return F_prop, M_prop


sim_length = int(tmax/dt + 1)
v_s = np.zeros(sim_length)
v_a = np.zeros(sim_length)
n_e = np.zeros(sim_length)
a_s = np.zeros(sim_length)
M_b = np.zeros(sim_length)
M_trm = np.zeros(sim_length)
n_p = np.zeros(sim_length)
P_b = np.zeros(sim_length)
f_flux = np.zeros(sim_length)
F_prop = np.zeros(sim_length)
M_prop = np.zeros(sim_length)
R = np.zeros(sim_length)
alpha_p = np.zeros(sim_length)
alpha_e = np.zeros(sim_length)
clutch = np.zeros(sim_length)
J = np.zeros(sim_length)
X = np.zeros(sim_length)
Y = np.zeros(sim_length)

# Initial values
v_s[0] = 2.
n_e[0] = 200/60
n_p[0] = n_e[0] / i_gb

for i, t in enumerate(time):
    if i >= sim_length-1:
        continue
    X[i] = X_func(t)
    Y[i] = Y_func(t)
    clutch[i] = clutch_func(t, X[i], n_e[i])

    M_b[i], P_b[i], f_flux[i] = engine(n_e[i], X[i])

    if clutch[i]:
        n_p[i], M_trm[i] = gearbox(n_e[i], M_b[i])

    v_a[i] = v_s[i] * (1 - w)
    J[i] = (v_a[i] / (n_p[i] * D_p))

    F_prop[i], M_prop[i] = propeller(n_p[i], J[i])

    if clutch[i]:
        alpha_p[i] = (M_trm[i] - M_prop[i])/(2*np.pi*(I_prop+I_eng*i_gb))
        alpha_e[i] = alpha_p[i] * i_gb
        n_e[i+1] = n_e[i] + alpha_e[i] * dt
    # else:
    #     alpha_p[i] = M_prop/(2*np.pi*I_prop)
    #     alpha_e[i] = M_b/(2*np.pi*I_eng)
    #     n_e[i+1] = n_e[i] + alpha_e[i] * dt
    #     n_p[i+1] = n_p[i] + alpha_p[i] * dt

    R[i] = resistance(v_s[i], Y[i])

    a_s[i] = (F_prop[i] - R[i]) / m_ship
    v_s[i+1] = v_s[i] + a_s[i] * dt

    # if not (0.2 < X[i] < 1):
    #     print("Fuel rack buiten toegestaan bereik.")
    
    # if n_e[i] > 900/60:
    #     print("Maximaal toerental overschreden.")

dist_travelled = integrate.cumulative_trapezoid(v_s, dx=dt, initial=0)
fuel_used = integrate.cumulative_trapezoid(f_flux, dx=dt, initial=0)

fig = plt.figure(figsize=(10, 20))
ax1 = fig.add_subplot(6, 1, 1)  # fig.add_subplot(#rows, #cols, #plot)
#shipspeed
ax1.plot(time[5:sim_length-5], v_s[5:sim_length-5])
ax1.set(title='Ship Propulsion Output',
        ylabel='Ship speed [m/s]',
        xlabel='Time [s]')
ax1.grid()
# distance traveled
ax2 = fig.add_subplot(6, 1, 2)  # fig.add_subplot(#rows, #cols, #plot)
ax2.plot(time[5:sim_length-5], dist_travelled[5:sim_length-5])
ax2.set(title='Ship Distance Travelled',
        ylabel='Distance travelled [m]',
        xlabel='Time [s]')
ax2.grid()
# brandstofverbruik
ax3 = fig.add_subplot(6, 1, 3)  # fig.add_subplot(#rows, #cols, #plot)
ax3.plot(time[5:sim_length-5], fuel_used[5:sim_length-5])
ax3.set(title='Fuel Consumption over Time',
        ylabel='fuel consumption [g]',
        xlabel='Time [s]')
ax3.grid()
# fuelrack
ax4 = fig.add_subplot(6, 1, 4)
ax4.plot(time[5:sim_length-5], X[5:sim_length-5]*100)
ax4.set(title='Fuel Rack over Time',
        ylabel='Fuel rack [%]',
        xlabel='Time [s]')
ax4.grid()
fig.tight_layout()
fig.savefig('grafieken/resultaat_vaarsim_mt1457.png')

fig6, ax6 = plt.subplots()
ax6.plot(time[5:sim_length-5], n_e[5:sim_length-5], label="Engine RPS")
ax6.plot(time[5:sim_length-5], n_p[5:sim_length-5], label="Propeller RPS")
ax6.set(title='RPS over Time',
        xlabel='Time [s]',
        ylabel='RPS [Hz]')
ax6.legend()
ax6.grid()
fig6.tight_layout()
fig6.savefig("grafieken/n_p-n_e.png")

fig7, ax7 = plt.subplots()
ax7.plot(v_s[5:sim_length-5], R[5:sim_length-5]/1000)
ax7.set(title='Resistance over ship velocity',
        xlabel='Ship velocity [m/s]',
        ylabel='Resistance [kN]')
ax7.grid()
fig7.tight_layout()
fig7.savefig("grafieken/Snelheid-weerstand.png")

fig8, ax8 = plt.subplots()
ax8.plot(v_a[5:sim_length-5], F_prop[5:sim_length-5]/1000)
ax8.set(title='Thrust over advance velocity',
        xlabel='Advance Velocity [m/s]',
        ylabel='Thrust [kN]')
ax8.grid()
fig8.tight_layout()
fig8.savefig("grafieken/v_advance-thrust.png")

fig9, ax9 = plt.subplots()
ax9.plot(n_p[5:sim_length-5], M_prop[5:sim_length-5])
ax9.plot(n_p[5:sim_length-5], M_trm[5:sim_length-5])
ax9.plot(n_p[5:sim_length-5], alpha_p[5:sim_length-5])
ax9.set(title='Propellor torque over propellor RPM',
        xlabel='propellor RPS [Hz]',
        ylabel='Torque [kNm]')
ax9.grid()
fig9.tight_layout()
fig9.savefig("grafieken/n_p-M_prop.png")

fig10, ax10 = plt.subplots()
ax10.plot(n_e[5:sim_length-5], M_b[5:sim_length-5]/1000)
ax10.set(title='Engine torque over engine RPM',
        xlabel='RPS [Hz]',
        ylabel='Torque [kNm]')
ax10.grid()
fig10.tight_layout()
fig10.savefig("grafieken/n_e-M_b.png")

# fig11, ax11 = plt.subplots()
# ax11.plot(time[5:sim_length-5], P_E[5:sim_length-5]/1000, label="Towing power")
# ax11.plot(time[5:sim_length-5], P_d[5:sim_length-5]/1000, label="Propeller power")
# ax11.plot(time[5:sim_length-5], P_b[5:sim_length-5]/1000, label="Brake power")
# ax11.plot(time[5:sim_length-5], Q_f_l[5:sim_length-5]/1000, label="Thermal energy per ignition")
# ax11.set(title='Power over Time',
#         xlabel='Time [s]',
#         ylabel='Power [kW]')
# ax11.legend()
# ax11.grid()
# fig11.tight_layout()
# fig11.savefig("grafieken/P-t.png")

# fig12, ax12 = plt.subplots()
# ax12.plot(time[5:sim_length-5], np.zeros(sim_length-10) + eta_e*100, label="Engine efficiency")
# ax12.plot(time[5:sim_length-5], np.zeros(sim_length-10) + eta_h*100, label="Hull efficiency")
# ax12.plot(time[5:sim_length-5], np.zeros(sim_length-10) + eta_TRM*100, label="Transmission efficiency")
# ax12.plot(time[5:sim_length-5], eta_o[5:sim_length-5]*100, label="Open water propeller efficiency")
# ax12.set(title='Efficiency over Time',
#         xlabel='Time [s]',
#         ylabel='Efficiency [%]')
# ax12.legend()
# ax12.grid()
# fig12.tight_layout()
# fig12.savefig("grafieken/efficiency-time.png")

# fig13, ax13 = plt.subplots()
# ax13.plot(P_p[5:sim_length-5], eta_o[5:sim_length-5]*100, label="Open water propeller efficiency")
# ax13.plot(P_b[5:sim_length-5], np.zeros(sim_length-10) + eta_e*100, label="Engine efficiency")
# ax13.set(title='Efficiency over Power',
#         xlabel='Power [W]',
#         ylabel='Efficiency [%]')
# ax13.grid()
# ax13.legend()
# fig13.tight_layout()
# fig13.savefig("grafieken/efficiency-power.png")

fig14, ax14 = plt.subplots()
ax14.plot(time[5:sim_length-5], X[5:sim_length-5]*100, color="tab:blue")
ax14_2 = ax14.twinx()
ax14_2.plot(time[5:sim_length-5], n_e[5:sim_length-5], color="tab:orange")
ax14.set(title='Fuel rack and RPM over time')
ax14.set_xlabel('Time [s]')
ax14.set_ylabel('Fuel rack [%]', color="tab:blue")
ax14_2.set_ylabel("Engine RPM", color="tab:orange")
ax14.grid()
fig14.tight_layout()
fig14.savefig("grafieken/P-F-t.png")
