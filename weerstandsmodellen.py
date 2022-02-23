import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

V_m = np.array([0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1])   #m/s
R_tm = np.array([0.000,0.098,0.200,0.340,0.570,0.830,1.130,1.510,1.980,2.590,3.400])  #N
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
ax.set_xlim([0, 1.4])
ax.set_ylim([0.95, 1.8])
ax.set_xlabel(r"$F_r^4 / C_{fm}$")
ax.set_ylabel(r"$C_{tm} / C_{fm}$")
ax.scatter(x, y)
ax.plot(x, fity)

print("Fit is: {}x + {}".format(fit[0], fit[1]))

print("De vormfactor is: {}".format(fit[1]-1))

plt.show()
