#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Blah, blah
"""

from math import exp
import numpy as np
import sys

import parameters as p
import constants as c

__author__  = "Martin De Kauwe"
__version__ = "1.0 (18.04.2019)"
__email__   = "mdekauwe@gmail.com"


def photosynthesis(temp, co2, lambdax):

    # Q10 temperature response of CO2/O2 specificity ratio
    tau = lookup_Q10(0.57, 2600., temp)

    # Q10 temperature response of Michaelis constant for O2
    ko = lookup_Q10(1.2, 3.0e4, temp)

    # Q10 temperature response of Michaelis constant for CO2
    kc = lookup_Q10(2.1, 30.0, temp)

    # Calculate CO2 compensation point (partial pressure)
    # Eqn 8, Haxeltine & Prentice 1996a
    gamma_star = p.p02 / 2.0 / tau

    # Intercellular partial pressure of CO2 given stomatal opening (Pa)
    # Eqn 7, Haxeltine & Prentice 1996a
    pi_co2 = lambdax * co2 * c.PATMOS * c.CO2_CONV

    # Calculation of C1_C3, Eqn 4, Haxeltine & Prentice 1996a
    #
    # High-temperature inhibition modelled by suppression of LUE by decreased
    # relative affinity of rubisco for CO2 with increasing temperature
    # (Table 3.7, Larcher 1983)
    #
    # Notes:
    # - there is an error in Eqn 4, Haxeltine & Prentice 1996a (missing
    #   2.0* in denominator) which is fixed here (see Eqn A2, Collatz
    #   et al 1991)
    # - the explicit low temperature inhibition function has been removed
    #   and replaced by a temperature-dependent upper limit on V_m, see
    #   below
    # - the reduction in maximum photosynthesis due to leaf age (phi_c)
    #   has been removed
    # - alpha_a, accounting for reduction in PAR utilisation efficiency
    #   from the leaf to ecosystem level, appears in the calculation of
    #   apar (above) instead of here
    # - C_mass, the atomic weight of carbon, appears in the calculation
    #   of V_m instead of here
    c1 = (pi_co2 - gamma_star) / (pi_co2 + 2.0 * gamma_star) * p.alpha_c3

    # Calculation of C2_C3, Eqn 6, Haxeltine & Prentice 1996a
    c2 = (pi_co2 - gamma_star) / (pi_co2 + kc * (1.0 + p.p02 / ko))

    print(gamma_star, pi_co2, c1, c2)



def lookup_Q10(q10, base25, temp):

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


    temp = 15.0  # deg C
    co2 = 400.0  # umol mol-1
    lambdax = 1.0

    photosynthesis(temp, co2, lambdax)

    """
    xx = []
    temps = np.arange(1, 40)
    for t in temps:
        xx.append(photosynthesis(t))

    import matplotlib.pyplot as plt
    plt.plot(temps, xx)
    plt.show()
    """
