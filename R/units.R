# Nitrogen unit conversions ---------------------------------------------------
#
# The app (and users) work in umol N L-1; macrogrow expects mg N m-3.
# Both NO3 and NH4 carry one N atom, so the conversion is the same for both.
#
#   1 umol L-1 = 1 mmol m-3 = 14.067 mg N m-3

N_molar_mass <- set_units(14.067, "mg mmol-1")

umol_to_mgN <- function(x) {
  set_units(x, "umol L-1") |>
    set_units("mmol m-3") |>
    (\(v) v * N_molar_mass)() |>
    set_units("mg m-3") |>
    drop_units()
}

mgN_to_umol <- function(x) {
  set_units(x, "mg m-3") |>
    (\(v) v / N_molar_mass)() |>
    set_units("umol L-1") |>
    drop_units()
}
