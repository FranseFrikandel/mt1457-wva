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
    return 0.5 * (np.sin(t/20) + 1)

def Y_func(t):
    global tmax
    return 1

def clutch_func(t, X, n_e):
    """
    Functie die bepaalt of de motor wel of niet moet worden ontkoppelt. Kan op tijdsbasis of basis van andere
    parameters.
    """
    return True

simulatie_v2.simulation(X_func, Y_func, clutch_func, tmax, dt, v_s_init, n_e_init, "grafieken/exp4/")