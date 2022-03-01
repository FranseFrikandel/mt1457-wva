import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

V_m = np.array([0, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.3, 1.4, 1.5, 1.6])   #m/s
R_tm_offset = 0.4979
R_tm = np.array([0.4979, 0.4672, 0.4026, 0.1850, -0.3232, -0.9847, -1.9831, -3.8821, -5.1792, -6.7597, -9.6771, -14.7063])  #N
R_tm = - R_tm + R_tm_offset
g = 9.81            #m/s**2
l = 31.50/19        #m
viscositeit= 1.0811E-06
rho = 998.7780      # kg/m^3
S_m = 272.6 / (19*19)
Re = V_m * l / viscositeit
C_fm = 0.075 / (((np.log10(Re)-2))**2)
F_r = V_m/(np.sqrt(g*l))
C_tm = R_tm/(0.5*rho*(V_m**2)*S_m)

x = (F_r**4)/(C_fm)
y = C_tm/C_fm

tabel = pd.DataFrame({"$V_m$": V_m, "$R_{tm}$": R_tm, "$F_r$": F_r, "$Re$": Re, "$C_{fm}$": C_fm, 
                      "$C_{tm}$": C_tm, "$\\frac{F_r^4}{C_fm}$": x,"$\\frac{C_tm}{C_fm}$": y})

tabel.to_latex("tabel.txt", escape = False, index=False)

fit_from = 4
fit_to = -3

fit = np.polyfit(x[fit_from:fit_to], y[fit_from:fit_to], 1)

fity = fit[0]*x + fit[1]

fig, ax = plt.subplots()
# ax.set_xlim([0, 2])
# ax.set_ylim([0.95, 1.8])
ax.set_xlabel(r"$F_r^4 / C_{fm}$")
ax.set_ylabel(r"$C_{tm} / C_{fm}$")
ax.scatter(np.concatenate((x[0:fit_from], x[fit_to:-1])), 
           np.concatenate((y[0:fit_from], y[fit_to:-1])), color="tab:orange")

ax.scatter(x[fit_from:fit_to], y[fit_from:fit_to], color="tab:blue")
ax.plot(x, fity)

print("Fit is: {}x + {}".format(fit[0], fit[1]))

print("De vormfactor is: {}".format(fit[1]-1))

plt.show()
