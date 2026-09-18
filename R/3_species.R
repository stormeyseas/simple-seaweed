# Species definitions ---------------------------------------------------------
#
# Each entry holds:
#   label       display name
#   spec_params named numeric vector in the format expected by macrogrow
#   days        length of the growing period (days)
#   initials    starting state passed to grow_macroalgae()

site_params_default <- c(hz = 30, d_top = 1, hc = 3, farmA = 2500)

# Starting biomass for all species: 5 g/L, expressed in mg m-3 for macrogrow
starting_biomass <- gL_to_mgm3(5)

# Get species parameters
get_species_params <- function(file) {
  df <- read.csv(file)
  vals <- df$value
  names(vals) <- df$parameter
  vals
}

armata <- get_species_params("data/asparagopsis_armata_gametophyte.csv")
taxiformis <- get_species_params("data/asparagopsis_taxiformis_gametophyte.csv")
ecklonia <- get_species_params("data/ecklonia.csv")
ulva <- get_species_params("data/ulva.csv")

# To simplify a bit - I am setting self-shading to very low and consisten across species

species_list <- list(
  a_armata = list(
    label = "A. armata",
    spec_params = armata,
    days = 42,
    initials = c(
      Q_rel = 0.5,
      biomass_to_Nf(spec_params = armata, biomass = starting_biomass, Q_rel = 0.5, dry = F)
    )
  ),
  a_taxiformis = list(
    label = "A. taxiformis",
    spec_params = taxiformis,
    days = 42,
    initials = c(
      Q_rel = 0.5,
      biomass_to_Nf(spec_params = taxiformis, biomass = starting_biomass, Q_rel = 0.5, dry = F)
    )
  ),
  ulva = list(
    label = "Ulva",
    spec_params = ulva,
    days = 30,
    initials = c(
      Q_rel = 0.5,
      biomass_to_Nf(spec_params = ulva, biomass = starting_biomass, Q_rel = 0.5, dry = F)
    )
  ),
  ecklonia = list(
    label = "Ecklonia",
    spec_params = ecklonia,
    days = 90,
    initials = c(
      Q_rel = 0.5,
      biomass_to_Nf(spec_params = ecklonia, biomass = starting_biomass, Q_rel = 0.5, dry = F)
    )
  )
)

# Species display
species_pal <- c(
  "ecklonia" = "darkgoldenrod",
  "ulva" = "chartreuse4",
  "a_armata" = "steelblue",
  "a_taxiformis" = "brown2",
  "Ecklonia" = "darkgoldenrod",
  "Ulva" = "chartreuse4",
  "A. armata" = "steelblue",
  "A. taxiformis" = "brown2"
)

# Fixed (non-user) environmental settings
kW_default <- 0.6
# 10 cm/s in m/s
velocity_default <- 0.1

# Reference nutrient levels (umol N L-1, converted inside run_growth()) used
# only when pre-calculating the optimum velocity, so the result does not
# depend on user inputs.
ref_nitrate <- 10
ref_ammonium <- 5
