#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Blah, blah
"""

from math import exp
import numpy as np
import sys
import parameters as p

__author__  = "Martin De Kauwe"
__version__ = "1.0 (18.04.2019)"
__email__   = "mdekauwe@gmail.com"


def photosynthesis(temp):

    # Calculate CO2 compensation point (partial pressure)
    # Eqn 8, Haxeltine & Prentice 1996a
    gamma_star = p.PO2 / 2.0 / lookup_tau(0.57, 2600., temp)

    return(gamma_star)


def lookup_tau(q10, base25, temp):

    # Constants required for Q10 lookup tables used by photosynthesis
    Q10_MINTEMP = -70.	# minimum temperature ever (deg C)
    Q10_MAXTEMP = 70.	# maximum temperature ever (deg C)
    Q10_PRECISION = 0.01	# rounding precision for temperature
    Q10_NDATA = int( (Q10_MAXTEMP - Q10_MINTEMP) / Q10_PRECISION + 1.5 )

    data = np.zeros(Q10_NDATA)
    for i in range(0, Q10_NDATA):
        data[i] = base25 * q10**(Q10_MINTEMP + i*Q10_PRECISION - 25.0) / 10.0

    if temp < Q10_MINTEMP:
        temp = Q10_MINTEMP
    elif temp > Q10_MAXTEMP:
        temp = Q10_MAXTEMP

    i = int( (temp - Q10_MINTEMP) / Q10_PRECISION + 0.5 )

    return ( data[i] )

if __name__ == "__main__":


    temp = 15.0
    x = photosynthesis(temp)

    xx = []
    temps = np.arange(1, 40)
    for t in temps:
        xx.append(photosynthesis(t))

    import matplotlib.pyplot as plt
    plt.plot(temps, xx)
    plt.show()
