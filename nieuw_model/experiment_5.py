"""
Dit experiment bestand runt de simulatie met de parameters zoals origineel aangeleverd.
"""
import numpy as np
import simulatie_v2

tmax = 800
dt = 0.1

v_s_init = 6.5430
n_e_init = 800/60

def X_func(t):
    global tmax
    t_control = np.array([0, 0, 0.2*tmax, 0.4*tmax, 0.6*tmax,
                         0.8*tmax, tmax])
    X_parms = np.array([1, 1, 0.8, 0.6, 0.4, 0.2, 0.0])

    return max(np.interp(t, t_control, X_parms), 0.2)

def Y_func(t):
    global tmax
    t_control = np.array([0, 0, 0.2*tmax, 0.4*tmax, 0.6*tmax,
                         0.8*tmax, tmax])
    Y_parms = np.array([2, 4, 8, 16, 32, 64, 128])

    return np.interp(t, t_control, Y_parms)

def clutch_func(t, X, n_e):
    """
    Functie die bepaalt of de motor wel of niet moet worden ontkoppelt. Kan op tijdsbasis of basis van andere
    parameters.
    """
    return True

simulatie_v2.simulation(X_func, Y_func, clutch_func, tmax, dt, v_s_init, n_e_init, "grafieken/exp5/")