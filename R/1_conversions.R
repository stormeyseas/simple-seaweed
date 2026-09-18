# Unit conversions ------------------------------------------------------------
#
# Plain numeric conversion factors; the app deliberately avoids the `units`
# package so it can run in the browser via Shinylive/webR.

# Nitrogen: the app (and users) work in umol N L-1; macrogrow expects mg N m-3.
# Both NO3 and NH4 carry one N atom, so the conversion is the same for both.
#
#   1 umol L-1 = 1 mmol m-3 = 14.067 mg N m-3
N_molar_mass <- 14.067 # mg N per mmol

umol_to_mgN <- function(x) x * N_molar_mass
mgN_to_umol <- function(x) x / N_molar_mass

# Biomass: 1 g L-1 = 1e6 mg m-3
gL_to_mgm3 <- function(x) x * 1e6
mgm3_to_gL <- function(x) x / 1e6

# Mass: mg -> g
mg_to_g <- function(x) x / 1000
