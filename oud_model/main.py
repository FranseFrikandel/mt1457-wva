# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 11:02:26 2017

    Viskotter simulation file
    Version 1.0H
    J. Rodrigues Monteiro, based on matlab code from P. de Vos
    Delft University of Technology
    3ME / MTT / SDPO / ME

History:
    20171108I: initial python version    JRM
    20200108J: simps integration         EU
              no graphing endpoits...   EU
    20200228H: simps correctie           EU
"""
import numpy as np
import matplotlib.pyplot as plt

from scipy import integrate

# ----------- parameters for simulation --------------------------------------
tmax = 36000            # simulation time [s]
dt = 1                  # timestep [s]
# Time
mytime = np.linspace(0, tmax-1, tmax) # TODO: Stappen zijn niet correct opgedeeld in dt

# fuel properties
LHV = 42700             # Lower Heating Value [kJ/kg]
print('fuel properties loaded')

# water properties
rho_sw = 1025           # density of seawater [kg/m3]
viscositeit = 1.2872E-06
print('water properties loaded')

# ship data
m_ship = 358000         # ship mass [kg]
c1 = 1500               # resistance coefficient c1 in R = c1*vs^2
v_s0 = 6.5430           # ship design speed [m/s]
t = 0.1600              # thrust deduction factor[-]
w = 0.2000              # wake factor [-]
l = 31.5                # length of ship [m]
eta_h = (1-t)/(1-w)
print('ship data loaded')

# propellor data
D_p = 3                 # diameter of propellor [m]
K_T_a = -0.3821         # factor a in K_T = a*J + b [-]
K_T_b = 0.2885          # factor b in K_T = a*J + b [-]
K_T_factor = np.array((-0.3821, 0.2885))
K_Q_a = -0.03346        # factor a in K_Q = a*J + b [-]
K_Q_b = 0.0308          # factor b in K_Q = a*J + b [-]
eta_R = 1.0100          # relative rotative efficiency [-]
print('propellor data loaded')

# engine data
m_f_nom = 1.314762      # nominal fuel injection [g]
eta_e = 0.3800          # nominal engine efficiency [-]
i = 6                   # number of cylinders [-]
k_es = 2                # k-factor for engines based on nr.of strokes per cycle
P_b = np.zeros(tmax)    # engine power [kW]
P_b[0] = 960            # Nominal engine power [kW]
M_b = np.zeros(tmax)    # engine torque [Nm]
M_b[0] = P_b[0]*1000/2/np.pi/(900/60)  # ([P_b*1000/2/np.pi/n_eng_nom])
print('engine data loaded')

# gearbox data
eta_TRM = 0.9500        # transmission efficiency [-]
i_gb = 4.2100           # gearbox ratio [-]
I_tot = 200             # total mass of inertia of propulsion system [kg*m^2]
print('gearbox data loaded')

# initial values
in_p = 3.2830           # initial rpm
iv_t_control = np.array([0, 0.1*tmax, 0.2*tmax, 0.5*tmax,
                         0.6*tmax, 0.7*tmax, tmax])
X_parms = np.array([0.85, 0.85, 0.3, 0.4, 0.4, 1, 1])   # % maximum fuelrack
Y_parms = np.array([1, 1, 1, 1, 1, 1, 1])               # disturbance factor

# simulation control parameters
xvals = np.linspace(0, tmax-1, tmax)
ov_X_set = np.interp(xvals, iv_t_control, X_parms)
ov_Y_set = np.interp(xvals, iv_t_control, Y_parms)

# # initial values
# in_p = 3.2830                                                           # initial rpm
# iv_t_control = np.array([0, 100, 200, 300, 400, 500, 600, 700, 800])
# X_parms = (np.sin(iv_t_control)+1)*0.5                                         # % maximum fuelrack
# Y_parms = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1])                         # disturbance factor

# # simulation control parameters
# xvals = np.linspace(0, tmax-1, tmax)
# # ov_X_set = np.interp(xvals, iv_t_control, X_parms)
# # ov_Y_set = np.interp(xvals, iv_t_control, Y_parms)
# ov_X_set = (np.sin(0.05*mytime)+1)*0.5
# ov_Y_set = np.interp(xvals, iv_t_control, Y_parms)

# Tabel van schip weerstandscoefficienten vs snelheid
# C_Ws = np.array([[0. , 0.43588989, 0.87177979, 1.74355958, 2.61533937,
#        3.48711915, 4.35889894, 5.23067873, 5.66656863, 6.10245852,
#        6.53834842, 6.97423831], [0, -0.00096818, -0.00128629, -0.00125762,  0.00016455,
#        0.00061338,  0.00130253,  0.00298342,  0.00390739,  0.00489266,
#        0.00713291,  0.01095216]])
# k = 0.23261693620133994

# --------- Start van de funtie definities

# def R_schip(snelheid_schip):
#     global Y, C_Ws, k, l_s, viscositeit, rho_sw
#     C_Ws_cur = np.interp(snelheid_schip, C_Ws[0], C_Ws[1])
#     Re_s_cur = (snelheid_schip*l_s)/viscositeit
#     C_Fs_cur = 0.075 / (((np.log10(Re_s_cur)-2))**2)
#     C_Ts_cur = (1+k)*C_Fs_cur + C_Ws_cur
#     weerstand = Y*C_Ts_cur*0.5*rho_sw*(snelheid_schip**2)
#     return weerstand

def R_schip(snelheid_schip):
      global Y, c1
      weerstand =  Y * c1 * snelheid_schip**2
      return weerstand

# -------- Make arrays -------------------------------------------------------

# Velocity of the ship [m/s]
v_s = np.zeros(tmax)
v_s[0] = v_s0
# Distance traveled [m]
s = np.zeros(tmax)
# Advance velocity [m/s]
v_a = np.zeros(tmax)
v_a[0] = (1-w) * v_s0
# Rpm propellor [Hz]
n_p = np.zeros(tmax)
n_p[0] = in_p
# Rpm diesel engine [Hz]
n_e = np.zeros(tmax)
n_e[0] = 900/60         # Nominal engine speed in rotations per second [Hz]
# Resistance [N]
R = np.zeros(tmax)
Y = ov_Y_set[0]
R[0] = R_schip(v_s0)
# Acceleration ship [m/s^2]
sum_a = np.zeros(tmax)
# Acceleration propellor[1/s^2]
sum_dnpdt = np.zeros(tmax)
m_flux_f = np.zeros(tmax)
out_fc = np.zeros(tmax)

M_Trm = np.zeros(tmax)            # M_B * i_gb * eta_TRM
KT = np.zeros(tmax)               # Thrust coefficient [-]
KQ = np.zeros(tmax)               # Torque coefficient [-]
Rsp = np.zeros(tmax)              # Resistance propelled situation [N]
F_prop = np.zeros(tmax)           # Thrust power [N]
M_prop = np.zeros(tmax)           # Torque [Nm]
P_O = np.zeros(tmax)              # Open water propellor power
P_p = np.zeros(tmax)              # Propellor power [kW]
P_b = np.zeros(tmax)              # Engine brake power [kW]
P_T = np.zeros(tmax)              # Thrust power [kW]
P_E = np.zeros(tmax)              # Engine power [kW]
P_d = np.zeros(tmax)              # Delivered propellor power
J = np.zeros(tmax)                # Advance ratio [-]
Q_f_l = np.zeros(tmax)            # Heat energy per power stroke

#---------- Run simulation -----------------------------------------------

for k in range(tmax-1):
    # advance ratio
    J[k+1] = ((v_a[k] / n_p[k]) / D_p)
    # Thrust and torque
    F_prop[k] = ((((J[k+1] * K_T_a) + K_T_b) *
                  n_p[k] ** 2) * rho_sw * D_p ** 4)
    M_prop[k] = (((((J[k+1] * K_Q_a) + K_Q_b) *
                  n_p[k] ** 2) * rho_sw * D_p ** 5) / eta_R)
    KT[k+1] = J[k+1] * K_T_a + K_T_b
    KQ[k+1] = J[k+1] * K_Q_a + K_Q_b
    P_O[k+1] = ((((J[k+1] * K_Q_a) + K_Q_b) *
                n_p[k] ** 2) * rho_sw * D_p ** 5) * n_p[k] * 2 * np.pi
    P_p[k+1] = M_prop[k] * n_p[k] * 2 * np.pi
    # Calculate acceleration from resulting force --> ship speed & tr.distance
    sum_a[k+1] = ((F_prop[k] - (R[k] / (1-t)))/m_ship)
    #v_s_new = (np.trapz(sum_a[k:k+2], dx=0.01)) + v_s[k]
    v_s[k+1]  = integrate.simps(sum_a[:k+2], dx=0.01) + v_s0 # TODO: Waarom dx=0.01?
    #v_s[k+1] = v_s_new
    Rsp[k+1] = R[k] / (1-t)
    # Traveled distance
    s[k+1] = s[k] + v_s[k+1] * dt
    # Advance velocity
    v_a[k+1] = v_s[k+1] * (1 - w)
    P_T[k+1] = F_prop[k] * v_a[k+1]
    # Resistance
    Y = ov_Y_set[k]
    R[k+1] = R_schip( v_s[k+1])
    P_E[k+1] = v_s[k+1] * R[k+1]
    # Calculate acceleration from resulting force --> propellor np
    sum_dnpdt[k+1] = ((M_b[k] * i_gb * eta_TRM) - M_prop[k])/(2*np.pi*I_tot)
    n_p[k+1] = integrate.simps(sum_dnpdt[:k+2], dx=0.01) + n_p[0]
    # Engine speed
    n_e[k+1] = n_p[k+1] * i_gb
    # Fuel rack
    X = ov_X_set[k]
    m_flux_f[k+1] = (X * m_f_nom * n_e[k+1]) * i / k_es
    # Fuel consumption
    out_fc[k+1] =  integrate.simps(m_flux_f[:k+2], dx=0.01)+out_fc[0]
    Q_f = X * m_f_nom * LHV
    Q_f_l[k] = Q_f
    W_e = Q_f * eta_e
    # Brake power
    P_b[k+1] = (W_e * n_e[k+1] * i) / k_es
    # Engine torque
    M_b[k+1] = P_b[k+1] / (2 * np.pi * n_e[k+1])
    M_Trm[k+1] = M_b[k+1] * i_gb * eta_TRM

# EU just to be sure
v_s[0]=v_s0
v_s[1]=v_s0

# Extra calculations (Not necessary to do inside the loop, and way faster when doing directly too)
eta_o = P_T / P_O
P_d = P_O/eta_R

# -------------- Plot Figure -------------------------------------------------

# create figure with four subplots
fig = plt.figure(figsize=(10, 20))
ax1 = fig.add_subplot(6, 1, 1)  # fig.add_subplot(#rows, #cols, #plot)
#shipspeed
ax1.plot(mytime[5:tmax-5], v_s[5:tmax-5])
ax1.set(title='Ship Propulsion Output',
        ylabel='Ship speed [m/s]',
        xlabel='Time [s]')
ax1.grid()
# distance traveled
ax2 = fig.add_subplot(6, 1, 2)  # fig.add_subplot(#rows, #cols, #plot)
ax2.plot(mytime[5:tmax-5], s[5:tmax-5])
ax2.set(title='Ship Distance Travelled',
        ylabel='Distance travelled [m]',
        xlabel='Time [s]')
ax2.grid()
# brandstofverbruik
ax3 = fig.add_subplot(6, 1, 3)  # fig.add_subplot(#rows, #cols, #plot)
ax3.plot(mytime[5:tmax-5], out_fc[5:tmax-5])
ax3.set(title='Fuel Consumption over Time',
        ylabel='fuel consumption [g]',
        xlabel='Time [s]')
ax3.grid()
# fuelrack
ax4 = fig.add_subplot(6, 1, 4)
ax4.plot(mytime[5:tmax-5], ov_X_set[5:tmax-5])
ax4.set(title='Fuel Rack over Time',
        ylabel='Fuel rack [%]',
        xlabel='Time [s]')
ax4.grid()
fig.tight_layout()
fig.savefig('grafieken/resultaat_vaarsim_mt1457.png')

fig6, ax6 = plt.subplots()
ax6.plot(mytime[5:tmax-5], n_e[5:tmax-5], label="Engine RPS")
ax6.plot(mytime[5:tmax-5], n_p[5:tmax-5], label="Propeller RPS")
ax6.set(title='RPS over Time',
        xlabel='Time [s]',
        ylabel='RPS [Hz]')
ax6.legend()
ax6.grid()
fig6.tight_layout()
fig6.savefig("grafieken/n_p-n_e.png")

fig7, ax7 = plt.subplots()
ax7.plot(v_s[5:tmax-5], R[5:tmax-5]/1000)
ax7.set(title='Resistance over ship velocity',
        xlabel='Ship velocity [m/s]',
        ylabel='Resistance [kN]')
ax7.grid()
fig7.tight_layout()
fig7.savefig("grafieken/Snelheid-weerstand.png")

fig8, ax8 = plt.subplots()
ax8.plot(v_a[5:tmax-5], F_prop[5:tmax-5]/1000)
ax8.set(title='Thrust over advance velocity',
        xlabel='Advance Velocity [m/s]',
        ylabel='Thrust [kN]')
ax8.grid()
fig8.tight_layout()
fig8.savefig("grafieken/v_advance-thrust.png")

fig9, ax9 = plt.subplots()
ax9.plot(n_p[5:tmax-5], M_prop[5:tmax-5]/1000)
ax9.set(title='Propellor torque over propellor RPM',
        xlabel='propellor RPS [Hz]',
        ylabel='Torque [kNm]')
ax9.grid()
fig9.tight_layout()
fig9.savefig("grafieken/n_p-M_prop.png")

fig10, ax10 = plt.subplots()
ax10.plot(n_e[5:tmax-5], M_b[5:tmax-5]/1000)
ax10.set(title='Engine torque over engine RPM',
        xlabel='RPS [Hz]',
        ylabel='Torque [kNm]')
ax10.grid()
fig10.tight_layout()
fig10.savefig("grafieken/n_e-M_b.png")

fig11, ax11 = plt.subplots()
ax11.plot(mytime[5:tmax-5], P_E[5:tmax-5]/1000, label="Towing power")
ax11.plot(mytime[5:tmax-5], P_d[5:tmax-5]/1000, label="Propeller power")
ax11.plot(mytime[5:tmax-5], P_b[5:tmax-5]/1000, label="Brake power")
ax11.plot(mytime[5:tmax-5], Q_f_l[5:tmax-5]/1000, label="Thermal energy per ignition")
ax11.set(title='Power over Time',
        xlabel='Time [s]',
        ylabel='Power [kW]')
ax11.legend()
ax11.grid()
fig11.tight_layout()
fig11.savefig("grafieken/P-t.png")

fig12, ax12 = plt.subplots()
ax12.plot(mytime[5:tmax-5], np.zeros(tmax-10) + eta_e*100, label="Engine efficiency")
ax12.plot(mytime[5:tmax-5], np.zeros(tmax-10) + eta_h*100, label="Hull efficiency")
ax12.plot(mytime[5:tmax-5], np.zeros(tmax-10) + eta_TRM*100, label="Transmission efficiency")
ax12.plot(mytime[5:tmax-5], eta_o[5:tmax-5]*100, label="Open water propeller efficiency")
ax12.set(title='Efficiency over Time',
        xlabel='Time [s]',
        ylabel='Efficiency [%]')
ax12.legend()
ax12.grid()
fig12.tight_layout()
fig12.savefig("grafieken/efficiency-time.png")

fig13, ax13 = plt.subplots()
ax13.plot(P_p[5:tmax-5], eta_o[5:tmax-5]*100, label="Open water propeller efficiency")
ax13.plot(P_b[5:tmax-5], np.zeros(tmax-10) + eta_e*100, label="Engine efficiency")
ax13.set(title='Efficiency over Power',
        xlabel='Power [W]',
        ylabel='Efficiency [%]')
ax13.grid()
ax13.legend()
fig13.tight_layout()
fig13.savefig("grafieken/efficiency-power.png")

fig14, ax14 = plt.subplots()
ax14.plot(mytime[5:tmax-5], ov_X_set[5:tmax-5]*100, color="tab:blue")
ax14_2 = ax14.twinx()
ax14_2.plot(mytime[5:tmax-5], n_e[5:tmax-5], color="tab:orange")
ax14.set(title='Fuel rack and RPM over time')
ax14.set_xlabel('Time [s]')
ax14.set_ylabel('Fuel rack [%]', color="tab:blue")
ax14_2.set_ylabel("Engine RPM", color="tab:orange")
ax14.grid()
fig14.tight_layout()
fig14.savefig("grafieken/P-F-t.png")
