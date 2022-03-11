import numpy as np
import matplotlib.pyplot as plt

def nieuw_motormodel(n_e, X):
    n_cyl = 6                   # number of cylinders [-]
    k_es = 2
    m_f_nom = 1.314762      # nominal fuel injection [g]
    eta_td = 0.52
    LHV = 42700
    n_nom = 900/60

    Q_f = X * m_f_nom * LHV
    Q_cooling = 1908.8 + 7635.2 * X
    Q_in = Q_f - Q_cooling
    W_loss_mech = 711.1 + 1659.3 * (n_e/n_nom)
    W_i = Q_in * eta_td
    W_e = W_i - W_loss_mech
    P_b = W_e * n_e * n_cyl / k_es
    M_b = P_b / (2 * np.pi * n_e)

    eta_e = W_e / Q_f
    eta_q = (Q_f - Q_cooling) / Q_f
    eta_mech = W_e / W_i

    return M_b, P_b, eta_e, eta_q, eta_mech

def oud_motormodel(n_e, X):
    n_cyl = 6                   # number of cylinders [-]
    k_es = 2 
    m_f_nom = 1.314762      # nominal fuel injection [g]
    eta_e = 0.38
    LHV = 42700

    Q_f = X * m_f_nom * LHV
    W_e = Q_f * eta_e
    P_b = (W_e * n_e * n_cyl) / k_es
    # Engine torque
    M_b = P_b / (2 * np.pi * n_e)

    return M_b, P_b

ySize = 100
xSize = 100

# Sweep over alle fuel racks
M_b_n = np.zeros((ySize, xSize))
P_b_n = np.zeros((ySize, xSize))
eta_e_n = np.zeros((ySize, xSize))
eta_q_n = np.zeros((ySize, xSize))
eta_mech_n = np.zeros((ySize, xSize))

X_array = np.linspace(0.2, 1, xSize)
n_e_array = np.linspace(200/60, 900/60, ySize)

for i, X in enumerate(X_array):
    for j, n_e in enumerate(n_e_array):
        M_b_n[j, i], P_b_n[j, i], eta_e_n[j, i], eta_q_n[j, i], eta_mech_n[j, i] = nieuw_motormodel(n_e, X)

X2D, Y2D = np.meshgrid(X_array*100, n_e_array*60)

fig_power = plt.figure()
ax_power = fig_power.add_subplot(111)
ax_power.set_xlabel("Fuel rack [%]")
ax_power.set_ylabel("Engine RPM")
im_power = ax_power.pcolormesh(X2D, Y2D, P_b_n/1000, cmap="hot", vmin=0, vmax=1000, shading="gouraud")
fig_power.colorbar(im_power, ticks=[100, 200, 300, 400, 500, 600, 700, 800, 900], orientation='vertical')
fig_power.savefig("motormodel/vermogen.pdf")

fig_torque = plt.figure()
ax_torque = fig_torque.add_subplot(111)
ax_torque.set_xlabel("Fuel rack [%]")
ax_torque.set_ylabel("Engine RPM")
im_torque = ax_torque.pcolormesh(X2D, Y2D, M_b_n, cmap="hot", vmin=0, vmax=15000, shading="gouraud")
fig_torque.colorbar(im_torque, ticks=[1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000], 
                    orientation='vertical')
fig_torque.savefig("motormodel/koppel.pdf")

# Rendementen op 550 RPM.

cur_rpm_idx = int(ySize * 350/700)

fig_efficencies = plt.figure()
ax_efficencies = fig_efficencies.add_subplot(111)
ax_efficencies.plot(X_array*100, eta_e_n[cur_rpm_idx, :]*100, label=r"$\eta_e$")
ax_efficencies.plot(X_array*100, eta_q_n[cur_rpm_idx, :]*100, label=r"$\eta_q$")
ax_efficencies.plot(X_array*100, eta_mech_n[cur_rpm_idx, :]*100, label=r"$\eta_{mech}$")
ax_efficencies.plot(X_array*100, np.zeros(xSize) + 52, label=r"$\eta_{td}$")
ax_efficencies.plot(X_array*100, np.zeros(xSize) + 100, label=r"$\eta_{comb}$")
ax_efficencies.plot(X_array*100, np.zeros(xSize) + 38, label=r"$\eta_{e}$ origineel")
ax_efficencies.set_xlabel("Fuel rack [%]")
ax_efficencies.set_ylabel("Efficiency [%]")
ax_efficencies.set_title("Engine efficiency vs fuel rack")
ax_efficencies.grid()
box = ax_efficencies.get_position()
ax_efficencies.set_position([box.x0, box.y0 + box.height * 0.15,
                 box.width, box.height * 0.9])
ax_efficencies.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
          fancybox=True, shadow=True, ncol=6)
fig_efficencies.savefig("motormodel/rendement.pdf")
