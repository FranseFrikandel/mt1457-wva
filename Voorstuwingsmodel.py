import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import integrate

# ----------- parameters for simulation --------------------------------------
tmax = 36000           # simulation time [s]
dt = 1                  # timestep [s]

# fuel properties
LHV = 42700             # Lower Heating Value [kJ/kg]
print('fuel properties loaded')

# water properties
rho_sw = 1025           # density of seawater [kg/m3]
print('water properties loaded')

# ship data
m_ship = 358000         # ship mass [kg]
c1 = 1500               # resistance coefficient c1 in R = c1*vs^2
v_s0 = 6.5430           # ship design speed [m/s]
t = 0.1600              # thrust deduction factor[-]
w = 0.2000              # wake factor [-]
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
M_b[0] = P_b[0]*1000/2/math.pi/(900/60)  # ([P_b*1000/2/math.pi/n_eng_nom])
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

# --------- Start van de funtie definities

def R_schip(snelheid_schip):
      global Y, c1
      weerstand =  Y * c1 * snelheid_schip**2
      return weerstand


# -------- Make arrays -------------------------------------------------------

# Time
mytime = np.linspace(0, tmax-1, tmax)
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
J = np.zeros(tmax)                # Advance ratio [-]


# ------------- Run simulation -----------------------------------------------
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
                n_p[k] ** 2) * rho_sw * D_p ** 5) * n_p[k] * 2 * math.pi
    P_p[k+1] = M_prop[k] * n_p[k] * 2 * math.pi
    # Calculate acceleration from resulting force --> ship speed & tr.distance
    sum_a[k+1] = ((F_prop[k] - (R[k] / (1-t)))/m_ship)
    #v_s_new = (np.trapz(sum_a[k:k+2], dx=0.01)) + v_s[k]
    v_s[k+1]  = integrate.simps(sum_a[:k+2], dx=0.01)+v_s0
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
    sum_dnpdt[k+1] = ((M_b[k] * i_gb * eta_TRM) - M_prop[k])/(2*math.pi*I_tot)
    n_p[k+1] = integrate.simps(sum_dnpdt[:k+2], dx=0.01) + n_p[0]    
    # Engine speed
    n_e[k+1] = n_p[k+1] * i_gb
    # Fuel rack
    X = ov_X_set[k]
    m_flux_f[k+1] = (X * m_f_nom * n_e[k+1]) * i / k_es
    # Fuel consumption
    out_fc[k+1] =  integrate.simps(m_flux_f[:k+2], dx=0.01)+out_fc[0]
    Q_f = X * m_f_nom * LHV
    W_e = Q_f * eta_e
    # Brake power
    P_b[k+1] = (W_e * n_e[k+1] * i) / k_es
    # Engine torque
    M_b[k+1] = P_b[k+1] / (2 * math.pi * n_e[k+1])
    M_Trm[k+1] = M_b[k+1] * i_gb * eta_TRM

# EU just to be sure
v_s[0]=v_s0
v_s[1]=v_s0

# -- Deel 2 -- 

v = np.array([0.000,0.200,0.400,0.600,0.800,1.000,1.200,1.400,1.500])           #m/s
n = np.array([15, 15, 15, 15, 15, 15, 15, 15, 15])                              #omw/s
T = np.array([19.459,17.472,15.413,12.987,10.465,7.402,4.265,1.053,-0.605])     #N
Q = np.array([0.277,0.255,0.234,0.207,0.176,0.141,0.104,0.063,0.042])           #Nm
rho = 1000                                                                      #kg/m^3
D = 0.133                                                                       #m

J_nieuw = v/(n*D)
K_T = T/(rho * (n**2) * (D**4))
K_Q = Q/(rho * (n**2) * (D**5))
eta = (1/(2*math.pi))*((K_T*J_nieuw)/K_Q)

fit_eta = np.polyfit(J_nieuw, eta,2)
fit_K_T = np.polyfit(J_nieuw, K_T,2)
fit_K_Q = np.polyfit(J_nieuw, K_Q,2)
#################### oorspronkelijke lineaire plot
K_T_a = -0.3821         # factor a in K_T = a*J + b [-]
K_T_b = 0.2885          # factor b in K_T = a*J + b [-]
K_Q_a = -0.03346        # factor a in K_Q = a*J + b [-]
K_Q_b = 0.0308          # factor b in K_Q = a*J + b [-]
#eta_R_lin = np.full((36000),eta_R )
K_T_lin = (K_T_a*J) + K_T_b
K_Q_lin = (K_Q_a*J) + K_Q_b
####################

J_plot = np.linspace(0, 0.75, 20)

K_T_Exp = fit_K_T[0] * J_plot**2 + fit_K_T[1] * J_plot + fit_K_T[2]
K_Q_Exp = fit_K_Q[0] * J_plot**2 + fit_K_Q[1] * J_plot + fit_K_Q[2]

# print("fit eta: {}x^2 + {}x + {}".format(fit_eta))
# print("fit K_t: {}x^2 + {}x + {}".format(fit_K_T))
# print("fit K_q: {}x^2 + {}x + {}".format(fit_K_Q))

plt.figure(figsize=())
plt.scatter(J_nieuw, K_T, color = 'r')
plt.scatter(J_nieuw, 10*K_Q, color = 'r')
plt.scatter(J_nieuw, eta, color = 'r')
plt.ylim(0,(0.6))
plt.xlim(0,0.8)
############################## moeten we de eta_R ook plotten >>> is dit de zelfde eta als die we uitrekenen?
plt.plot(J,K_T_lin, 'g', label ='K_T_lineair', linestyle = 'dashed')
plt.plot(J,10*K_Q_lin, 'g', label = 'K_Q_lineair', linestyle = ':')
#plot.plot(J,eta_R_lin,'g', label = '\u03B7_linear')


plt.plot(J_nieuw,K_T, 'r', label = 'K_T_meet',linestyle='dashed')
plt.plot(J_nieuw,10*K_Q,'r', label ='K_Q_meet',linestyle=':')
plt.plot(J_nieuw,eta,'r', label = '\u03B7_meet')
plt.plot(J_plot, K_T_Exp)
plt.plot(J_plot, K_Q_Exp*10)
plt.xlabel('Snelheids graad J')
plt.ylabel('K_T, 10*K_Q, \u03B7')
plt.title('Openwater schroef diagram') 
plt.grid()
plt.legend()
plt.savefig("Openwater schroef diagram.png")
