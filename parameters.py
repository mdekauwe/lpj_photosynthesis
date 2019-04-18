# O2 partial pressure (Pa)
O2 = 2.09e3#2.09e4

# atmospheric pressure (Pa)
PATMOS = 1e5

# intrinsic quantum efficiency of CO2 uptake, C3 plants
alpha_c3 = 0.08

# leaf respiration as fraction of maximum rubisco, C3 plants
BC3 = 0.015

# colimitation (shape) parameter
theta = 0.7

# optimal (maximum) lambda in C3 plants
lambda_max = 0.8

ko25 = 3.0e4   # value of ko at 25 deg C
kc25 = 30.0    # value of kc at 25 deg C
tau25 = 2600.0 # value of tau at 25 deg C
q10ko = 1.2    # q10 for temperature-sensitive parameter ko
q10kc = 2.1    # q10 for temperature-sensitive parameter kc
q10tau = 0.57  # q10 for temperature-sensitive parameter tau

# Parameters common to all temperate trees
pstemp_min = -2.  # approximate low temperature limit for photosynthesis (deg C)
pstemp_low = 15.  # approximate lower range of temperature optimum for
                  # photosynthesis (deg C)
pstemp_high = 25. # approximate upper range of temperature optimum for
                  # photosynthesis (deg C)
pstemp_max = 38.  # maximum temperature limit for photosynthesis (deg C)

lat = 44.
lon = -120.0
