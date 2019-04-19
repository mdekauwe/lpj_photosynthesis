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


def photosynthesis(Tleaf, apar, co2, lambdax, vm=None):
    """
    Total daily gross photosynthesis

    Calculation of total daily gross photosynthesis and leaf-level net daytime
    photosynthesis given degree of stomatal closure (as parameter lambda).
    Includes implicit scaling from leaf to plant projective area basis.
    Adapted from Farquhar & von Caemmerer (1982) photosynthesis model, as
    simplified by Collatz et al (1991), Collatz et al (1992),
    Haxeltine & Prentice (1996a,b) and Sitch et al. (2000)
    """

    # Inhibition for low and high temperatures. High-temperature inhibition
    # modelled by suppression of LUE by decreased relative affinity of rubisco
    # for CO2 with increasing temperature
    # (Table 3.7, Larcher 1983)
    tscal = calc_temp_inhibition(Tleaf)
    #tscal = 1.0

    # Kinetic parameters (Kc, Ko & tau) are modelled using a Q10 reln

    # Q10 temperature response of CO2/O2 specificity ratio
    # units: -
    tau = lookup_Q10(p.q10tau, p.tau25, Tleaf)

    # Q10 temperature response of Michaelis constant for O2
    # units: kPa
    ko = lookup_Q10(p.q10ko, p.ko25, Tleaf)

    # Q10 temperature response of Michaelis constant for CO2
    # units: Pa
    kc = lookup_Q10(p.q10kc, p.kc25, Tleaf)

    # Calculate CO2 compensation point (partial pressure)
    # Eqn 8, Haxeltine & Prentice 1996a
    # units: Pa
    gamma_star = p.O2 / 2.0 / tau

    #ko = p.ko25 * p.q10ko**((Tleaf - 25.0) / 10.0)
    #kc = p.kc25 * p.q10kc**((Tleaf - 25.0) / 10.0)
    #tau = p.tau25 * p.q10tau**((Tleaf - 25.0) / 10.0)

    # Intercellular partial pressure of CO2 given stomatal opening
    # Eqn 7, Haxeltine & Prentice 1996a
    # units: Pa
    pi_co2 = lambdax * co2

    # Calculation of C1_C3, Eqn 4, Haxeltine & Prentice 1996a

    # Notes:
    # - there is an error in Eqn 4, Haxeltine & Prentice 1996a (missing
    #   2.0* in denominator) which is fixed here (see Eqn A2, Collatz
    #   et al 1991)
    # units: Pa
    c1 = (pi_co2 - gamma_star) / (pi_co2 + 2.0 * gamma_star) * p.alpha_c3

    # Calculation of C2_C3, Eqn 6, Haxeltine & Prentice 1996a
    # units: Pa
    c2 = (pi_co2 - gamma_star) / (pi_co2 + kc * (1.0 + p.O2 / ko))

    if vm is None:
        vm = vmax(Tleaf, apar, c1, c2, tscal)

    # Daily leaf respiration
    # Eqn 10, Haxeltine & Prentice 1996a
    # units: mol m-2 s-1
    Rd = vm * p.BC3

    # PAR-limited photosynthesis rate
    # Eqn 3, Haxeltine & Prentice 1996a
    # units: mol m-2 s-1
    je = c1 * tscal * apar

    # Rubisco-activity limited photosynthesis rate
    # Eqn 5, Haxeltine & Prentice 1996a
    # units: mol m-2 s-1
    jc = c2 * vm

    # Gross photosynthesis, A
    # Eqn 2, Haxeltine & Prentice 1996a
    # Notes: - there is an error in Eqn 2, Haxeltine & Prentice 1996a (missing
    # 			theta in 4*theta*je*jc term) which is fixed here
    # units: mol m-2 s-1
    A = (je + jc - \
            np.sqrt((je + jc) * (je + jc) - 4.0 * p.theta * je * jc)) / \
            (2.0 * p.theta)

    # Net photosynthesis, An
    # units: mol m-2 s-1
    An = A - Rd

    return An * c.MOL_TO_UMOL

def vmax(temp, apar, c1, c2, tscal):
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

def arrh(k25, Ea, Tk):
    """ Temperature dependence of kinetic parameters is described by an
    Arrhenius function.

    Parameters:
    ----------
    k25 : float
        rate parameter value at 25 degC or 298 K
    Ea : float
        activation energy for the parameter [J mol-1]
    Tk : float
        leaf temperature [deg K]

    Returns:
    -------
    kt : float
        temperature dependence on parameter

    References:
    -----------
    * Medlyn et al. 2002, PCE, 25, 1167-1179.
    """
    return k25 * np.exp((Ea * (Tk - 298.15)) / (298.15 * c.RGAS * Tk))


if __name__ == "__main__":

    from get_days_met_forcing import get_met_data

    #
    ## Met data ...
    #
    doy = 180.
    (par, tair, vpd) = get_met_data(p.lat, p.lon, doy)

    #par = np.mean(par.reshape(-1, 2), axis=1)
    #par *= 1800.0/par.max() #* c.SEC_TO_HR
    #tair = np.mean(tair.reshape(-1, 2), axis=1)
    #tair = np.ones(len(par)) * 25.

    # umol mol-1 to Pa
    co2 = 400.0 * p.patm * c.CO2_CONV

    # Ratio of intercellular to ambient partial pressure of CO2
    lambdax = p.lambda_max

    # Convert Vcmax mol m-2 s-1
    vm = 60. * c.UMOL_TO_MOL

    fpar = 1.0 - np.exp(-p.k * p.LAI)

    A = np.zeros(len(par))
    An = np.zeros(len(par))
    je = np.zeros(len(par))
    jc = np.zeros(len(par))

    for i in range(len(par)):

        # Scale fractional PAR absorption at plant projective area level (FPAR)
        # to fractional absorption at leaf level (APAR)
        # Eqn 4, Haxeltine & Prentice 1996a
        # units: mol m-2 s-1
        apar = (par[i] * c.UMOL_TO_MOL) * fpar * p.alpha_a

        An[i] = photosynthesis(tair[i], apar, co2, lambdax, vm)

    plt.plot(An, label="An")
    plt.legend(numpoints=1, loc="best")
    plt.ylabel("Photosynthesis ($\mathrm{\mu}$mol m$^{-2}$ s$^{-1}$)")
    plt.xlabel("Hour of day")
    plt.show()
