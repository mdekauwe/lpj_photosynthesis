#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Blah, blah
"""

from math import exp
import numpy as np
import sys
import matplotlib.pyplot as plt

import parameters as p
import constants as c

__author__  = "Martin De Kauwe"
__version__ = "1.0 (18.04.2019)"
__email__   = "mdekauwe@gmail.com"


def photosynthesis(temp, apar, co2, day_length, lambdax, vm=None):
    """
    Total daily gross photosynthesis

    Calculation of total daily gross photosynthesis and leaf-level net daytime
    photosynthesis given degree of stomatal closure (as parameter lambda).
    Includes implicit scaling from leaf to plant projective area basis.
    Adapted from Farquhar & von Caemmerer (1982) photosynthesis model, as
    simplified by Collatz et al (1991), Collatz et al (1992),
    Haxeltine & Prentice (1996a,b) and Sitch et al. (2000)
    """
    # Q10 temperature response of CO2/O2 specificity ratio
    tau = lookup_Q10(0.57, 2600., temp)

    # Q10 temperature response of Michaelis constant for O2
    ko = lookup_Q10(1.2, 3.0e4, temp)

    # Q10 temperature response of Michaelis constant for CO2
    kc = lookup_Q10(2.1, 30.0, temp)

    tscal = calc_temp_inhibition(temp)

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

    if vm is None:
        vm = vmax(temp, apar, day_length, c1, c2, tscal)

    # Calculation of daily leaf respiration
    # Eqn 10, Haxeltine & Prentice 1996a
    rd_g = vm * p.BC3

    # PAR-limited photosynthesis rate (gC/m2/h)
    # Eqn 3, Haxeltine & Prentice 1996a
    je = c1 * tscal * apar * c.CMASS * c.CQ #/ day_length

    # Rubisco-activity limited photosynthesis rate (gC/m2/h)
    # Eqn 5, Haxeltine & Prentice 1996a
    jc = c2 * vm #/ 24.0

    # Calculation of daily gross photosynthesis
    # Eqn 2, Haxeltine & Prentice 1996a
    # Notes: - there is an error in Eqn 2, Haxeltine & Prentice 1996a (missing
    # 			theta in 4*theta*je*jc term) which is fixed here
    # g c m-2 h-1
    agd_g = (je + jc - \
                np.sqrt((je + jc) * (je + jc) - 4.0 * p.theta * je * jc)) / \
                (2.0 * p.theta) #* day_length

    return agd_g

def vmax(temp, apar, day_length, c1, c2, tscal):
    # Calculation of non-water-stressed rubisco capacity assuming leaf nitrogen
    # not limiting (Eqn 11, Haxeltine & Prentice 1996a)
    #
    s =  24.0 / day_length * p.BC3

    # Calculation of sigma is based on Eqn 12 Haxeltine & Prentice 1996a
    sigma = np.sqrt(max(0., 1.0 - (c2 - s) / (c2 - p.theta * s)))

    # maximum daily rate of net photosynthesis, g C m-2 d-1
    arg1 = 1.0 / p.BC3
    arg2 = c.CMASS * c.CQ * c1 / c2 * tscal * apar
    arg3 = 2.0 * p.theta * s * (1. - sigma) - s + c2 * sigma
    vm = arg1 * arg2 * arg3

    # Conversion factor in calculation of leaf nitrogen: includes conversion of:
    #   - Vm from gC/m2/day to umolC/m2/sec
    #   - nitrogen from mg/m2 to kg/m2
    conv = 1.0 / (c.SEC_TO_HR * day_length * c.CMASS)

    vm *= conv

    return vm

def calc_temp_inhibition(temp):
    # Calculate temperature-inhibition coefficient
    #
    # This function (tscal) is mathematically identical to function tstress
	# in LPJF. In contrast to earlier versions of modular LPJ and LPJ-GUESS,
	# it includes both high- and low-temperature inhibition.
    k1 = (p.pstemp_min + p.pstemp_low) / 2.

    tscal = (1.0 - 0.01 * \
             np.exp(4.6 / (p.pstemp_max - p.pstemp_high) * \
             (temp - p.pstemp_high))) / \
			 (1.0 + np.exp((k1 - temp) / (k1 - p.pstemp_min) * 4.6))

    return tscal

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

    from get_days_met_forcing import get_met_data

    doy = 180.
    #
    ## Met data ...
    #
    (par, tair, vpd) = get_met_data(p.lat, p.lon, doy)

    par = np.mean(par.reshape(-1, 2), axis=1)
    par *= 1800.0/par.max() * c.UMOL_TO_J * c.SEC_TO_HR
    tair = np.mean(tair.reshape(-1, 2), axis=1)

    #day_length = 12.0
    #temp = 15.0  # deg C
    #par = 1000.0 * (c.SEC_TO_HR * day_length) # total daily PAR today (J/m2/day)
    co2 = 400.0  # umol mol-1


    # ratio of intercellular to ambient partial pressure of CO2
    lambda_max = 0.8
    lambdax = lambda_max


    # umol m-2 s-1 -> g m-2 h-1
    vm = 40. * c.CMASS * c.SEC_TO_HR

    fpar = 0.6

    # In sub-daily mode daylength should be 24h, to obtain values in daily units
    day_length = 24.
    a = np.zeros(len(par))
    for i in range(len(par)):

        # Scale fractional PAR absorption at plant projective area level (FPAR)
        # to fractional absorption at leaf level (APAR)
        # Eqn 4, Haxeltine & Prentice 1996a
        apar = par[i] * fpar;

        a[i] = photosynthesis(tair[i], apar, co2, day_length, lambdax, vm)

    print(np.sum(a))
    plt.plot(a)
    plt.show()
