# Pre-calculation helpers -----------------------------------------------------

# Thin wrapper: grow_macroalgae() needs every environmental driver as a
# vector of length(t), so recycle the constant conditions here.
#
# `nitrate` and `ammonium` are taken in umol N L-1 (user units) and converted
# to mg N m-3 for macrogrow. The returned conc_* columns are converted back
# to umol N L-1; uptake columns (up_Ni, up_Am) stay in mg N.
run_growth <- function(sp, temperature, salinity, light, velocity, nitrate, ammonium, kW = kW_default, site_params = site_params_default) {
  n <- sp$days
  grow_macroalgae(
    t = seq_len(n),
    temperature = rep(temperature, n),
    salinity = rep(salinity, n),
    light = rep(light, n),
    kW = rep(kW, n),
    velocity = rep(velocity, n),
    nitrate = rep(umol_to_mgN(nitrate), n),
    ammonium = rep(umol_to_mgN(ammonium), n),
    site_params = site_params,
    spec_params = sp$spec_params,
    initials = sp$initials
  ) |>
    as.data.frame() |>
    transform(
      conc_nitrate = mgN_to_umol(conc_nitrate),
      conc_ammonium = mgN_to_umol(conc_ammonium)
    )
}

# Surface light (umol photons m-2 s-1) that maximises I_lim at the starting biomass. I_lim peaks below 1 and then declines (photoinhibition), so we maximise rather than solve for I_lim == 1.
optimal_light <- function(sp, kW = kW_default, site_params = site_params_default, interval = c(1, 5000)) {
  opt <- optimize(
    function(I) {
      I_lim(
        Nf = sp$initials[["Nf"]],
        kW = kW,
        I = I,
        spec_params = sp$spec_params,
        site_params = site_params
      )
    },
    interval = interval, maximum = TRUE
  )
  c(light = opt$maximum, I_lim = opt$objective)
}

# Velocity (m/s) that maximises final biomass under the species' optimal
# temperature and salinity, with reference nutrient levels. Higher velocity
# reduces flow attenuation but increases biomass loss, hence the trade-off.
optimal_velocity <- function(sp, light, interval = c(0.01, 1)) {
  final_biomass <- function(v) {
    out <- run_growth(
      sp,
      temperature = sp$spec_params[["T_opt"]],
      salinity = sp$spec_params[["S_opt"]],
      light = light,
      velocity = v,
      nitrate = ref_nitrate,
      ammonium = ref_ammonium
    )
    tail(out$B_ww.mg, 1)
  }
  opt <- optimize(final_biomass, interval = interval, maximum = TRUE)
  c(velocity = opt$maximum)
}

# Attach pre-calculated settings to every species. Velocity is fixed at
# velocity_default rather than optimised (see optimal_velocity() above).
precalc_species <- function(species_list, velocity = velocity_default) {
  lapply(species_list, function(sp) {
    lt <- optimal_light(sp)
    sp$light <- unname(lt[["light"]])
    sp$I_lim_at_start <- unname(lt[["I_lim"]])
    sp$velocity <- velocity
    sp
  })
}

# Summary table of the fixed settings, for display in the app
settings_table <- function(species) {
  do.call(rbind, lapply(species, function(sp) {
    data.frame(
      Species = sp$label,
      `Days` = sp$days,
      `Start Nf` = sp$initials[["Nf"]],
      `Light` = round(sp$light),
      `I_lim (start)` = round(sp$I_lim_at_start, 3),
      `Velocity (m/s)` = round(sp$velocity, 3),
      `kW` = kW_default,
      check.names = FALSE
    )
  })) |>
    `rownames<-`(NULL)
}
