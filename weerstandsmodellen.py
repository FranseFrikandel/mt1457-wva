import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json

V_m = np.array([0, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.3, 1.4, 1.5, 1.6])   #m/s
V_s = V_m * np.sqrt(19)
R_Tm_offset = 0.4979
R_Tm = np.array([0.4979, 0.4672, 0.4026, 0.1850, -0.3232, -0.9847, -1.9831, -3.8821, -5.1792, -6.7597, -9.6771, -14.7063])  #N
d_voor = np.array([-12.1711, -12.8886, -13.3017, -14.2612, -16.6058, -19.0421, -22.8294, -27.7378, -33.9478, -34.7307, -37.5462, -38.6969])
d_achter = np.array([-4.2874, -5.0345, -5.3309, -5.5694, -6.7273, -6.2312, -8.0008, -9.6798, -14.0001, -11.4722, -14.4286, -19.5151])
d_gem = (d_voor + d_achter)/2
R_Tm = - R_Tm + R_Tm_offset
g = 9.81            #m/s**2

l_s=31.5
l_m = 31.50/19        #m
viscositeit_m = 1.0811E-06 
viscositeit_s = 1.2872E-06
rho_m = 998.7780      # kg/m^3
rho_s = 1026.6376 
S_s = 272.6
S_m = S_s / (19*19)

# omschaling
F_r = V_m/(np.sqrt(g*l_m))

Re_s= V_s * l_s / viscositeit_s
Re_m = V_m * l_m / viscositeit_m
C_Fm = 0.075 / (((np.log10(Re_m)-2))**2)
C_Tm = R_Tm/(0.5*rho_m*(V_m**2)*S_m)

#-=--=-=-=-=-=-=-=-=-==-=-=-=-=-=-==-=-==-=-=-=-=-==-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=--=-=-=-=-=-=-=-=-==-=-=-=-=-=-==-

fit_from = 4
fit_to = -3

x = (F_r**4)/(C_Fm)
y = C_Tm/C_Fm
fit = np.polyfit(x[fit_from:fit_to], y[fit_from:fit_to], 1)
fity = fit[0]*x + fit[1]

k = fit[1] - 1

#-=--=-=-=-=-=-=-=-=-==-=-=-=-=-=-==-=-==-=-=-=-=-==-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=--=-=-=-=-=-=-=-=-==-=-=-=-=-=-==-

C_Wm = C_Tm-(1+k)*C_Fm
C_Ws = C_Wm
Re_s = (V_s*l_s)/viscositeit_s
C_Fs = 0.075 / (((np.log10(Re_s)-2))**2)
C_Ts = (1+k)*C_Fs + C_Ws
R_Ts = C_Ts*0.5*rho_s*(V_s**2)*S_s


tabel_model = pd.DataFrame({"$V_m$": V_m, "$R_{tm}$": R_Tm, "$F_r$": F_r, "$Re$": Re_m, "$C_{Fm}$": C_Fm, 
    "$C_{Wm}$": C_Tm - C_Fm, "$C_{Tm}$": C_Tm, "$\\frac{F_r^4}{C_fm}$": x,"$\\frac{C_tm}{C_fm}$": y})

tabel_schip = pd.DataFrame({"$V_s$": V_s, "$R_{ts}$": R_Ts, "$F_r$": F_r, "$Re$": Re_s, "$C_{fs}$": C_Fs, 
    "$C_{ts}$": C_Ts})#, "$\\frac{F_r^4}{C_fm}$": x,"$\\frac{C_tm}{C_fm}$": y})

tabel_model.to_latex("tabel_model.txt", escape = False, index=False)
tabel_schip.to_latex("tabel_schip.txt", escape = False, index=False)

print(f"Fit is: {fit[0]}x + {fit[1]}")
print(f"De vormfactor is: {k}")

fig1, ax1 = plt.subplots()
# ax1.set_xlim([0, 2])
# ax1.set_ylim([0.95, 1.8])
ax1.set_xlabel(r"$F_r^4 / C_{fm}$")
ax1.set_ylabel(r"$C_{tm} / C_{fm}$")
ax1.scatter(np.concatenate((x[0:fit_from], x[fit_to:-1])), 
           np.concatenate((y[0:fit_from], y[fit_to:-1])), color="tab:orange")
ax1.scatter(x[fit_from:fit_to], y[fit_from:fit_to], color="tab:blue")
ax1.plot(x, fity)
ax1.grid()
ax1.set_title("Prohaska plot")
fig1.savefig("weerstandsproef/prohaska.pdf")

fig2, ax2 = plt.subplots()
ax2.plot(V_m, d_voor, label="Inzinking voor")
ax2.plot(V_m, d_achter, label="Inzinking achter")
ax2.plot(V_m, d_gem, label="Inzinking gemiddeld")
ax2.set_xlabel("Modelsnelheid [m/s]")
ax2.set_ylabel("Inzinking [mm]")
ax2.legend()
ax2.grid()
ax2.set_title("Inzinking van het model")
fig2.savefig("weerstandsproef/inzinking.pdf")

# fig3, ax3 = plt.subplots()
# ax3.plot(V_s, R_Ts)
# ax3.legend()
# ax3.grid()
# fig3.savefig("weerstandsproef/weerstand.pdf")

fig4, ax4 = plt.subplots()
ax4.plot(V_m, R_Tm)
# ax4.legend()
ax4.grid()
ax4.set_ylabel("Weerstand [N]")
ax4.set_xlabel("Snelheid [m/s]")
ax4.set_title("Totale weerstand van het model")
fig4.savefig("weerstandsproef/weerstand_model.pdf")

Y = 1
c1 = 1500

# with fileHandle as open("weerstandsdata.dat"):
#     fileHandle.write(C_Ws)
#     fileHandle.write(k)

def R_schip_oud(snelheid_schip):
    global Y, c1
    weerstand =  Y * c1 * snelheid_schip**2
    return weerstand

def R_schip_nieuw(snelheid_schip):
    global Y, C_Ws, V_s, k, l_s, viscositeit_s, rho_s
    C_Ws_cur = np.interp(snelheid_schip, V_s, C_Ws)
    Re_s_cur = (snelheid_schip*l_s)/viscositeit_s
    C_Fs_cur = 0.075 / (((np.log10(Re_s_cur)-2))**2)
    C_Ts_cur = (1+k)*C_Fs_cur + C_Ws_cur
    weerstand = Y*C_Ts_cur*0.5*rho_s*(snelheid_schip**2)*S_s
    return weerstand

R_oud = np.zeros(70)
R_nieuw = np.zeros(70)
snelheden = np.linspace(0.1, 7, 70)

for i, snelheid in enumerate(snelheden):
    R_oud[i] = R_schip_oud(snelheid)
    R_nieuw[i] = R_schip_nieuw(snelheid)

fig5, ax5 = plt.subplots()
ax5.plot(snelheden, R_oud, label="Oud weerstandsmodel")
ax5.plot(snelheden, R_nieuw, label="Nieuw weerstandsmodel")
ax5.legend()
ax5.grid()
ax5.set_xlabel("Schipsnelheid [m/s]")
ax5.set_ylabel("Weerstand [N]")
ax5.set_title("Het oude en nieuwe weerstandsmodel")
fig5.savefig("weerstandsproef/weerstand_modellen.pdf")

fig6, ax6 = plt.subplots()
ax6.plot(V_m[1:-1], C_Wm[1:-1], label=r"$C_Wm$")
ax6.plot(V_m[1:-1], C_Tm[1:-1], label=r"$C_Tm$")
ax6.plot(V_m[1:-1], C_Fm[1:-1], label=r"$C_Fm$")
ax6.set_xlabel("Modelsnelheid [m/s]")
ax6.set_ylabel("Weerstandscoefficient")
ax6.legend()
ax6.grid()
fig6.savefig("weerstandsproef/weerstandscoefficient.pdf")

fig6, ax6 = plt.subplots()
ax6.plot(V_m[1:-1], C_Wm[1:-1], label=r"$C_Wm$")
ax6.plot(V_m[1:-1], C_Tm[1:-1], label=r"$C_Tm$")
ax6.plot(V_m[1:-1], C_Fm[1:-1], label=r"$C_Fm$")
ax6.set_xlabel("Modelsnelheid [m/s]")
ax6.set_ylabel("Weerstandscoefficient")
ax6.legend()
ax6.grid()
fig6.savefig("weerstandsproef/weerstandscoefficient.pdf")

fig7, ax7 = plt.subplots()
ax7.bar(1, C_Ws[-1], label=r"$C_{Ws}$")
ax7.bar(1, C_Fs[-1], label=r"$C_{Fs}$", bottom=C_Wm[-1])
ax7.bar(2, C_Ts[-1], label=r"$C_{Ts}$")
ax7.set_ylabel("Weerstandscoefficient")
ax7.axes.xaxis.set_ticklabels([])
ax7.legend()
# ax7.grid()

fig7.savefig("weerstandsproef/weerstandscoefficient-bar.pdf")
