"""
Dit experiment bestand runt de simulatie met de parameters zoals origineel aangeleverd.
"""
import numpy as np
import simulatie_v2

tmax = 3000
dt = 0.1

v_s_init = 6
n_e_init = 800/60

def X_func(t):
    global tmax
    t_control = np.array([0, 0.1*tmax, 0.2*tmax, 0.4*tmax,
                         0.8*tmax, tmax])
    X_parms = np.array([0.5, 0.5, 0.25, 0.01, 0.01, 0.01])

    return max(np.interp(t, t_control, X_parms), 0.2)

def Y_func(t):
    global tmax
    t_control = np.array([0, 0.1*tmax, 0.2*tmax, 0.5*tmax,
                         0.6*tmax, 0.7*tmax, tmax])
    Y_parms = np.array([1, 1, 1, 1, 1, 1, 1])

    return np.interp(t, t_control, Y_parms)

def clutch_func(t, X, n_e):
    """
    Functie die bepaalt of de motor wel of niet moet worden ontkoppelt. Kan op tijdsbasis of basis van andere
    parameters.
    """
    if X < 0.205:
        return False
    return True

simulatie_v2.simulation(X_func, Y_func, clutch_func, tmax, dt, v_s_init, n_e_init, "grafieken/exp_koppeling_1/")