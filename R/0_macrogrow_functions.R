# ---- Source: C:/Users/treimer/Documents/R-temp-files/macrogrow/R/algae_height.R ----
#' Macroalgae height
#' 
#' @description
#' Calculates macroalgae height as 
#' \deqn{h_m = \left(\frac{N_f}{h_a}\right)^{h_b} + h_c}
#' where `h_m` is limited by a maximum of `h_max`. 
#'
#' @inheritParams Q_int
#' @param spec_params A vector of named numbers. Must include the parameters:
#' * `h_max`, maximum species height. Can be `NA`.
#' * `h_a`, `h_b` and `h_c`, parameters governing height change with `N_f`. If not supplied algae height will always be `h_max`.
#'
#' @details
#' Defaults are \eqn{h_a=1000}, \eqn{h_b=1} and \eqn{h_c=0}. 
#' Algae height therefore defaults to \eqn{N_f \times 10^{-3}} if no parameters are supplied.
#' 
#' @return a scalar of macroalgae height (m)
#' @export
#'
#' @examples 
#' my_species <- c(h_a = 750, h_b = 0.5, h_c = 0.01, h_max = 1)
#' Nf <- seq(100, 1000, 10)
#' 
#' # Calculate height for a range of fixed nitrogen values
#' height <- sapply(X = Nf, FUN = height, spec_params = my_species)
#' \dontrun{
#'   plot(Nf, height, type = "l")
#' }

height <- function(Nf, spec_params) {
  # Check which parameters are supplied
  h_a <- ifelse(is.na(spec_params['h_a']), 1000, spec_params['h_a'])
  h_b <- ifelse(is.na(spec_params['h_b']), 1000, spec_params['h_b'])
  h_c <- ifelse(is.na(spec_params['h_c']), 1000, spec_params['h_c'])
  hm <- (Nf/h_a)^h_b + h_c
  
  if (!is.na(spec_params['h_max'])) {
    hm <- min(hm, spec_params['h_max'])
  }
  return(unname(hm))
}

# ---- Source: C:/Users/treimer/Documents/R-temp-files/macrogrow/R/biomass_Nf_conversions.R ----
#' Convert biomass to Nf and Ns

#' @description
#' Converts wet or dry biomass to `N_f` and `N_s` via:
#' \deqn{B = \frac{N_f + N_s}{Q_{min}}}
#' and where the ratio between `N_f` and `N_s` is calculated via:
#' \deqn{\frac{Ns}{Nf} = \frac{Q_{int}}{Q_{min}} - 1} 
#'
#' @inheritParams Nf_to_biomass
#' @inheritParams Q_rel 
#' @param biomass starting biomass, mg m-3
#' @param Q_rel the non-dimensionalised relative internal nutrient quotient (\eqn{Q_{rel}}). Only one of \eqn{Q_{int}} or \eqn{Q_{rel}} need to be provided. If neither \eqn{Q_{int}} or \eqn{Q_{rel}} are provided the default if \eqn{Q_{rel}=0.5} will be used.
#' @param spec_params a vector of named numbers. Must include:
#'  * `DWWW` (if dry = F), the conversion from dry weight to wet weight
#'  * `Q_min`, the minimum internal nutrient quotient (mg gDW-1)
#'  * `Q_max`, the maximum internal nutrient quotient
#' @param dry whether dry (default) or wet biomass is provided, mg m-3

#' @return Nf, mg m-3
#' @export
#'
#' @examples 
#' my_species <- c(DWWW = 7.5, Q_min = 20, Q_max = 45)
#' starting_biomass <- 250 # mg m-3
#' 
#' # Using default Q_rel = 0.5
#' biomass_to_Nf(biomass = starting_biomass, spec_params = my_species, dry = T)
#' 
#' # Using a specific Q_int
#' biomass_to_Nf(biomass = starting_biomass, Q_int = 30, spec_params = my_species, dry = T)
#' 
#' @seealso [Nf_to_biomass()], [Q_rel()], [Q_int()]
#' 
biomass_to_Nf <- function(biomass, Q_int = NULL, Q_rel = 0.5, spec_params, dry = T) {
  # If only Q_rel is given, convert to Q_int
  if (is.null(Q_int)) {Q_int<- Q_int(Q_rel = Q_rel, spec_params = spec_params)}

  # Biomass must be dry
  if (dry == F) {biomass <- biomass/unname(spec_params['DWWW'])}

  Nf <- biomass * unname(spec_params['Q_min']) * 10^-3 # Nf comes directly from biomass
  Nf_Ns <- Q_int / unname(spec_params['Q_min']) - 1 # Ratio of Nf to Ns
  Ns <- Nf_Ns * Nf

  return(c(Nf = Nf, Ns = Ns))
}


#' Convert Nf to biomass
#'
#' @inheritParams Q_int
#' @inheritParams Q_rel 
#' @param dry logical, return dry or wet biomass. If dry = F, `spec_params['DWWW']` must be provided
#'
#' @details
#' Not all parameters need to be provided to this function. 
#' 
#' @return dry (or wet) biomass, mg m-3
#' @export
#'
#' @examples 
#' my_species <- c(DWWW = 7.5, Q_min = 20, Q_max = 45)
#' starting_Nf <- 150 # mg m-3
#' starting_Ns <- 50 # mg m-3
#' 
#' # Using default Q_rel = 0.5
#' Nf_to_biomass(Nf = starting_Nf, Ns = starting_Ns, spec_params = my_species, dry = T)
#' 
#' # Using a specific Q_int
#' Nf_to_biomass(biomass = starting_biomass, Q_int = 30, spec_params = my_species, dry = T)
#' 
#' @seealso [biomass_to_Nf()], [Q_rel()], [Q_int()]
#' 
Nf_to_biomass <- function(Nf, Ns, Q_int = NULL, Q_rel = 0.5, spec_params, dry = T) {
  # If only Q_rel is given, convert to Q_int
  if (is.null(Q_int)) {Q_int <- Q_int(Nf = Nf, Ns = Ns, Q_rel = Q_rel, spec_params = spec_params)}
  biomass <- ((Nf + Ns) / Q_int) * 10^3
  # If biomass is dry, convert
  if (dry == F) {biomass <- biomass * unname(spec_params['DWWW'])}
  return(biomass)
}

# ---- Source: C:/Users/treimer/Documents/R-temp-files/macrogrow/R/check_grow.R ----
#' @title Pre-check all parameters for grow_macroalgae()
#' 
#' @description
#' Check that all parameters, inputs and settings are correct for the `grow_macroalgae()` function. Gives a report on what needs to be modified for the main function to run smoothly. This is to avoid the main function slowing down to give endless warnings and messages. 
#'
#' @inheritParams grow_macroalgae
#'
#' @importFrom glue glue
#' @import rlang cli
#' 
#' @return printout of potential errors for main function
#' @export 
#' 
#' @details
#' Example csv with all the spec_params & site_params required?
#'
#' @seealso [grow_macroalgae()]
#' 
#' @examples "see here" link?
check_grow <- function(
    t = 1:30,
    temperature,
    salinity,
    light,
    kW,
    velocity,
    nitrate,
    ammonium,
    ni_uptake,
    am_uptake,
    site_params,
    spec_params,
    initials,
    sparse_output = T,
    other_constants = c(s = 0.0045, gam = 1.13, a2 = 0.2^2, Cb = 0.0025)
  ) {
  rlang::inform("Starting all checks...")
  
  # Start date
  if (any(!is.integer(t), !is.numeric(t))) {
    inform(c("x" = "Variable 't' must be a vector of integers, or numbers coercible to integers. (Hint: Convert dates using lubridate::yday)"))
  } else {
    inform(c("v" = "Timeseries looks good."))
  }

  # Check that input variables are present and look good
  helpcheck(t, temperature, "temperature")
  helpcheck(t, salinity, "salinity")
  helpcheck(t, light, "light")
  helpcheck(t, kW, "kW")
  helpcheck(t, light, "light")
  helpcheck(t, velocity, "velocity")
  helpcheck(t, nitrate, "nitrate")
  helpcheck(t, ammonium, "ammonium")

  # Check that site params is good
  essential_site_params <- c("farmA", "hz", "hc", "d_top")
  if (!all(essential_site_params %in% names(site_params))) {
    rlang::inform(paste0("x" = "Parameter '", essential_site_params[which(!essential_site_params %in% names(site_params))], "' is missing from site_params."))
  } else if (any(is.na(site_params[essential_site_params]))) {
    rlang::inform(paste0("x" = "Parameter '", essential_site_params[which(is.na(site_params[essential_site_params]))], "' in site_params cannot be NA."))
  } else {
    inform(c("v" = "Site params looks good."))
  }
  
  # Check validity of initial variables
  if (!'Nf' %in% names(initials)) {
    rlang::inform("x" = "Parameter 'Nf' is missing from initials. Calculate it using biomass_to_Nf() first if required.")
  } else if (is.na(initials['Nf'])) {
    rlang::inform("x" = "Parameter 'Nf' in initials cannot be NA.")
  } else {
    inform(c("v" = "Variable Nf in initials looks good."))
  }
  
  if (!'Q_int' %in% names(initials) & !'Q_rel' %in% names(initials)) {
    rlang::inform("x" = "Parameters 'Q_int' and 'Q_rel' are missing from initials. You must provide one of them to get macroalgae initial state. You can calculate Q_int using Q_int() from Nf and Ns.")
  } else if (is.na(initials['Q_int']) & is.na(initials['Q_rel'])) {
    rlang::inform("x" = "Parameters 'Q_int' and 'Q_rel' in initials cannot both be NA.")
  } else {
    inform(c("v" = "Variable Q_int/Q_rel in initials looks good."))
  }
  
  # Check that all essential species parameters are present
  essential_spec_params <- c('Q_min', 'Q_max', 'K_c', 'mu', 'a_cs', 'I_o', 'T_opt', 'T_min', 'T_max', 'DWWW')
  if (!all(essential_spec_params %in% names(spec_params))) {
    rlang::inform(paste0("x" = "Parameter '", essential_spec_params[which(!essential_spec_params %in% names(spec_params))], "' is missing from spec_params."))
  } else if (any(is.na(spec_params[essential_spec_params]))) {
    rlang::inform(paste0("x" = "Parameter '", essential_spec_params[which(is.na(spec_params[essential_spec_params]))], "' in spec_params cannot be NA."))
  } else {
    ess_spec_params <- T
  }
  
  # Optional parts, which depend on other inputs -------------------------------------------------------------------------------------------------------  
  # Validate nitrogen uptake parameters
  validate_uptake_params(ni_uptake, spec_params, "nitrate", "ni")
  validate_uptake_params(am_uptake, spec_params, "ammonium", "am")
  
  if (is.na(site_params['turbulence'])) {rlang::inform(
    ">" = "Parameter 'turbulence' is missing from site_params. Turbulence loss will not work regardless of spec_params supplied.")
  }
  
  
  # "D_m"
  # "D_ve"
  # "D_lo"
  # "D_mi"
  # "D_hi"
  if (is.na(spec_params['S_opt'])) {rlang::inform("x" = "Parameter 'S_opt' is missing from spec_params")}
  if (is.na(spec_params['S_min'])) {rlang::inform("x" = "Parameter 'S_min' is missing from spec_params")}
  if (is.na(spec_params['S_max'])) {rlang::inform("x" = "Parameter 'S_max' is missing from spec_params")}
  if (is.na(spec_params['h_max'])) {rlang::inform("x" = "Parameter 'h_max' is missing from spec_params")}
  if (any(is.na(c(spec_params['h_a'], spec_params['h_b'], spec_params['h_c'])))) {
    rlang::inform(">" = "Some height parameters are missing from spec_params. Unless supplied, parameters will default to h_a = 1000, h_b = 1, h_c = 0.")
  }
  
  # Get checks and warnings from u_c()
  # Get checks and warnings from T_lim()
  # Get checks and warnings from S_lim()
  # Get checks and warnings from Q_lim()
  # Get checks and warnings from I_lim()
  # Get checks and warnings from loss()
  # Get checks and warnings from get_uptake()
  
  if(!is.logical(sparse_output)) {
    rlang::inform("x" = "sparse_output is not logical - must be 'T' (default) or 'F'.")
  } else if (sparse_output == F) {
    rlang::inform(">" = "sparse_output is 'F', full outputs will be given.")
  } else {
    rlang::inform(">" = "sparse_output is 'T', truncated outputs will be given.")
  }
}

#' Check input vector compatibility with time series
#'
#' @description
#' Internal helper function that validates input vectors against the time series vector.
#' Performs comprehensive checks for NULL values, numeric type, length compatibility,
#' and missing values. Used by [`check_grow()`](R/check_grow.R:20) to validate all input variables.
#'
#' @param t Numeric vector representing the time series
#' @param vec Numeric vector to be validated against the time series
#' @param name Character string specifying the variable name for error messages
#'
#' @return Invisible boolean indicating validation success (TRUE) or failure (FALSE)
#' @keywords internal
#'
#' @details
#' The function performs the following validation steps:
#' \itemize{
#'   \item Checks for NULL values in inputs
#'   \item Validates that vec is numeric
#'   \item Ensures vec length matches t length
#'   \item Handles missing values using [`handle_missing_values()`](R/check_grow.R:183)
#' }
#'
#' @seealso [`handle_missing_values()`](R/check_grow.R:183), [`check_grow()`](R/check_grow.R:20)
helpcheck <- function(t, vec, name) {
  # Input validation for edge cases
  if (is.null(vec)) {
    rlang::inform(c("x" = glue::glue("Variable '{name}' is NULL.")))
    return(invisible(FALSE))
  }
  
  if (is.null(t) || length(t) == 0) {
    rlang::inform(c("x" = "Timeseries vector 't' is NULL or empty."))
    return(invisible(FALSE))
  }
  
  if (is.null(name) || nchar(name) == 0) {
    rlang::inform(c("x" = "Variable name is NULL or empty."))
    return(invisible(FALSE))
  }
  
  # Check if vector is numeric
  if (!is.numeric(vec)) {
    rlang::inform(c("x" = glue::glue("Variable '{name}' must be a numeric vector.")))
    return(invisible(FALSE))
  }
  
  # Check length compatibility
  vec_length <- length(vec)
  t_length <- length(t)
  
  if (vec_length != t_length) {
    rlang::inform(c("x" = glue::glue(
      "Variable '{name}' length ({vec_length}) does not match timeseries vector 't' length ({t_length})."
    )))
    return(invisible(FALSE))
  }
  
  # Check for missing values and handle based on variable importance
  na_count <- sum(is.na(vec))
  if (na_count > 0) {
    handle_missing_values(name, na_count, vec_length)
    return(invisible(FALSE))
  }
  
  # All checks passed
  rlang::inform(c("v" = glue::glue("Variable '{name}' looks good.")))
  return(invisible(TRUE))
}

#' Handle missing value messages based on variable type
#'
#' @description
#' Internal helper function that generates appropriate warning messages for missing values
#' based on the variable type and importance. Categorizes variables as essential or optional
#' and provides context-specific guidance for handling missing data.
#'
#' @param name Character string specifying the variable name
#' @param na_count Integer number of missing values in the variable
#' @param total_length Integer total length of the variable vector
#'
#' @return No return value. Function is called for its side effects of printing informative messages.
#' @keywords internal
#'
#' @details
#' Variables are categorized as:
#' \itemize{
#'   \item \strong{Essential variables}: temperature, ammonium, nitrate - cannot have missing values
#'   \item \strong{Optional variables}: salinity, light, kW, velocity - missing values disable specific limitations
#'   \item \strong{Other variables}: general warning about potential model performance impact
#' }
#'
#' The function calculates the percentage of missing values and provides variable-specific
#' guidance on the consequences of missing data.
#'
#' @seealso [`helpcheck()`](R/check_grow.R:136)
handle_missing_values <- function(name, na_count, total_length) {
  # Define variable categories and their handling
  essential_vars <- c("temperature", "ammonium", "nitrate")
  optional_vars <- list(
    "salinity" = "Salinity limitation will not be factored into growth.",
    "light" = "Light limitation will not be factored into growth.",
    "kW" = "Light limitation will not be factored into growth.",
    "velocity" = "Biomass loss due to current speed will not be factored into growth."
  )
  
  na_percentage <- round(na_count / total_length * 100, 1)
  base_message <- glue::glue("Variable '{name}' has {na_count} missing values.")
  
  if (name %in% essential_vars) {
    rlang::inform(c("x" = glue::glue(
      "{base_message} This input is essential and cannot have missing values."
    )))
  } else if (name %in% names(optional_vars)) {
    rlang::inform(c(">" = glue::glue(
      "{base_message} {optional_vars[[name]]} To use this limitation, ensure the input vector has no missing values."
    )))
  } else {
    rlang::inform(c(">" = glue::glue(
      "{base_message} This may affect model performance."
    )))
  }
}



#' Validate uptake parameters for nitrate or ammonium
#'
#' @param uptake_type Character string specifying uptake type ("MM" for Michaelis-Menton, "linear" for linear, or NA)
#' @param spec_params Named list containing species parameters
#' @param nutrient_name Character string for the nutrient name (e.g., "nitrate", "ammonium")
#' @param param_suffix Character string for parameter suffix (e.g., "ni", "am")
#'
#' @return Invisible boolean indicating success
#' @keywords internal
validate_uptake_params <- function(uptake_type, spec_params, nutrient_name, param_suffix) {
  # Define required parameters for each uptake type
  mm_params <- paste0(c("V_", "K_"), param_suffix)
  linear_params <- paste0(c("M_", "C_"), param_suffix)
  
  # Helper function to check if all parameters are present
  has_all_params <- function(params) {
    all(all(params %in% names(spec_params)), !is.na(spec_params[params]))
  }
  
  # Handle different uptake type scenarios
  if (is.na(uptake_type)) {
    # Auto-detect uptake type based on available parameters
    has_mm <- has_all_params(mm_params)
    has_linear <- has_all_params(linear_params)
    
    if (has_mm && has_linear) {
      rlang::inform(c(">" = glue::glue(
        "spec_params has provided parameters for both Michaelis-Menton and linear uptake for {nutrient_name}. ",
        "Uptake will default to Michaelis-Menton kinetics."
      )))
      return(invisible(TRUE))
    } else if (has_mm) {
      rlang::inform(c(">" = glue::glue(
        "spec_params has provided Michaelis-Menton parameters for {nutrient_name}. ",
        "Uptake will use Michaelis-Menton kinetics."
      )))
      return(invisible(TRUE))
    } else if (has_linear) {
      rlang::inform(c(">" = glue::glue(
        "spec_params has provided linear parameters for {nutrient_name}. ",
        "Uptake will use linear kinetics."
      )))
      return(invisible(TRUE))
    } else {
      rlang::inform(c("x" = glue::glue(
        "No uptake parameters provided for {nutrient_name}. ",
        "Either ({paste(mm_params, collapse = ', ')}) for Michaelis-Menton or ",
        "({paste(linear_params, collapse = ', ')}) for linear uptake must be provided in spec_params."
      )))
      return(invisible(FALSE))
    }
  } else if (uptake_type == "MM") {
    # Validate Michaelis-Menton parameters
    if (!has_all_params(mm_params)) {
      missing_params <- mm_params[!mm_params %in% names(spec_params)]
      rlang::inform(c("x" = glue::glue(
        "Variable '{param_suffix}_uptake' (uptake of {nutrient_name}) is set to Michaelis-Menton kinetics ",
        "but required parameters are missing: {paste(missing_params, collapse = ', ')}. ",
        "All parameters ({paste(mm_params, collapse = ', ')}) must be provided in spec_params."
      )))
      return(invisible(FALSE))
    }
  } else if (uptake_type == "linear") {
    # Validate linear parameters
    if (!has_all_params(linear_params)) {
      missing_params <- linear_params[!linear_params %in% names(spec_params)]
      rlang::inform(c("x" = glue::glue(
        "Variable '{param_suffix}_uptake' (uptake of {nutrient_name}) is set to linear kinetics ",
        "but required parameters are missing: {paste(missing_params, collapse = ', ')}. ",
        "All parameters ({paste(linear_params, collapse = ', ')}) must be provided in spec_params."
      )))
      return(invisible(FALSE))
    }
  } else {
    # Invalid uptake type
    rlang::inform(c("x" = glue::glue(
      "Invalid uptake type '{uptake_type}' for {nutrient_name}. ",
      "Must be 'MM' for Michaelis-Menton, 'linear' for linear kinetics, or NA for auto-detection."
    )))
    return(invisible(FALSE))
  }
  
  # All validations passed
  rlang::inform(c("v" = glue::glue("Variable '{param_suffix}_uptake' parameters look good.")))
  return(invisible(TRUE))
}

# ---- Source: C:/Users/treimer/Documents/R-temp-files/macrogrow/R/data.R ----
#' Asparagopsis armata model parameters
#'
#' A named vector containing physiological and morphological parameters for the macroalgae species *Asparagopsis armata* used in the growth model.
#'
#' @format A named numeric vector with 34 elements:
#' \describe{
#'   \item{`V_am`}{Maximum ammonium uptake rate (Michaelis-Menton uptake, mg gDW-1 d-1)}
#'   \item{`K_am`}{Half-saturation constant for ammonium uptake (Michaelis-Menton uptake, mg m-3)}
#'   \item{`M_am`}{Ammonium uptake rate slope (linear uptake, m3 gDW-1 d-1)}
#'   \item{`C_am`}{Ammonium uptake rate constant (linear uptake, m3 gDW-1 d-1)}
#'   \item{`V_ni`}{Maximum nitrate uptake rate (Michaelis-Menton uptake, mg gDW-1 d-1)}
#'   \item{`K_ni`}{Half-saturation constant for nitrate uptake (Michaelis-Menton uptake, mg m-3)}
#'   \item{`M_ni`}{Nitrate uptake rate slope (linear uptake, m3 gDW-1 d-1)}
#'   \item{`C_ni`}{Nitrate uptake rate constant (linear uptake, m3 gDW-1 d-1)}
#'   \item{`V_ot`}{Maximum other nitrogen uptake rate (Michaelis-Menton uptake, mg gDW-1 d-1)}
#'   \item{`K_ot`}{Half-saturation constant for other nitrogen uptake (Michaelis-Menton uptake, mg m-3)}
#'   \item{`M_ot`}{Other nitrogen uptake rate slope (linear uptake, m3 gDW-1 d-1)}
#'   \item{`C_ot`}{Other nitrogen uptake rate constant (linear uptake, m3 gDW-1 d-1)}
#'   \item{`Q_min`}{Minimum nitrogen quotient (mg gDW-1)}
#'   \item{`Q_max`}{Maximum nitrogen quotient (mg gDW-1)}
#'   \item{`K_c`}{Half-saturation constant for growth (mg gDW-1)}
#'   \item{`mu`}{Maximum growth rate (d-1)}
#'   \item{`D_m`}{Base loss rate (d-1)}
#'   \item{`D_ve`}{Loss rate with velocity (m-1 s d-1)}
#'   \item{`D_lo`}{Loss rate at low turbidity (d-1)}
#'   \item{`D_mi`}{Loss rate at medium turbidity (d-1)}
#'   \item{`D_hi`}{Loss rate at high turbidity (d-1)}
#'   \item{`a_cs`}{Carbon-specific self-shading constant (m mg-1)}
#'   \item{`I_o`}{Light saturation parameter (μmol photons m-2 s-1)}
#'   \item{`T_opt`}{Optimal temperature (°C)}
#'   \item{`T_min`}{Minimum temperature (°C)}
#'   \item{`T_max`}{Maximum temperature (°C)}
#'   \item{`S_opt`}{Optimal salinity (g L-1)}
#'   \item{`S_min`}{Minimum salinity (g L-1)}
#'   \item{`S_max`}{Maximum salinity (g L-1)}
#'   \item{`h_a`}{Height-controlling parameter a (d-1)}
#'   \item{`h_b`}{Height-controlling parameter b (d-1)}
#'   \item{`h_c`}{Height-controlling parameter c (d-1)}
#'   \item{`h_max`}{Maximum height (m)}
#'   \item{`DWWW`}{Dry weight to wet weight ratio (gWW gDW-1)}
#' }
#'
#' @source Parameters derived from literature values and used in *Asparagopsis armata* growth modeling here (insert citation).
#'
#' @seealso [grow_macroalgae()] for the main growth model function.
#'
#' @examples
#' # View all parameters
#' a_armata
#' 
#' # Access specific parameters
#' a_armata["mu"]  # Maximum growth rate
#' a_armata["h_max"]  # Maximum height
"a_armata"

#' Default site parameters
#'
#' A named vector containing default hydrodynamic and site characteristics used in the macroalgae growth model.
#'
#' @format A named numeric vector with 4 elements:
#' \describe{
#'   \item{`hz`}{Water depth (m)}
#'   \item{`d_top`}{Distance from surface to top of the macroalgae canopy (m)}
#'   \item{`hc`}{Canopy vertical width in the water column (m)}
#'   \item{`farmA`}{Total farm area (m2)}
#' }
#'
#' @examples
#' # View all site parameters
#' site_params
#' 
#' # Access specific parameters
#' site_params["hz"]  # Water depth
#' site_params["kW"]  # Water attenuation coefficient
"site_params"

#' Example environmental time series data
#'
#' A data frame containing synthetic environmental data over a full year (365 days) with seasonal patterns and random variation. This dataset demonstrates the expected format for environmental forcing data used in the macroalgae growth model.
#'
#' @format A data frame with 365 rows and 8 columns:
#' \describe{
#'   \item{t}{Time step (day of year)}
#'   \item{temperature}{Water temperature (°C)}
#'   \item{light}{Surface light availability (μmol photons m-2 s-1)}
#'   \item{kW}{Light attenuation coefficient in water (m-1)}
#'   \item{velocity}{Mean water velocity (m s-1)}
#'   \item{salinity}{Salinity (g L-1)}
#'   \item{nitrate}{Nitrate concentration (mg m-3)}
#'   \item{ammonium}{Ammonium concentration (mg m-3)}
#' }
#'
#' @details
#' The environmental data includes:
#' - Seasonal temperature variation (15°C ± 3.5°C with random noise)
#' - Light availability with seasonal patterns
#' - Seasonal variation in light attenuation coefficient
#' - Seasonal variation in water velocity 
#' - Salinity variation around 35 g/L
#' - Nitrogen concentrations (nitrate and ammonium) with seasonal patterns
#'
#' All variables include random variation to simulate natural environmental conditions.
#'
#' @seealso [grow_macroalgae()] for using this data in growth simulations.
#'
#' @examples
#' # View structure of environmental data
#' str(env)
#' 
#' # Plot temperature over time
#' plot(env$t, env$temperature, type = "l", 
#'      xlab = "Day of year", ylab = "Temperature (°C)")
#' 
#' # Summary of all environmental variables
#' summary(env)
"env"

# ---- Source: C:/Users/treimer/Documents/R-temp-files/macrogrow/R/drag_functions.R ----
#' Relative water attenuation within canopy
#' 
#' @description
#' Calculates relative water attenuation within the canopy, based on Plew, D. R. (2011). "Depth-Averaged Drag Coefficient for Modeling Flow through Suspended Canopies". Journal of Hydraulic Engineering, 137(2), 234–247.
#'
#' @param U0 incoming incident water velocity (m/s)
#' @param macro_state vector of named numbers. Must include:
#'  * `biomass`, macroalgae wet weight (g)
#'  * `hm`, algae height (m)
#' @param site_params vector of named numbers. Must include:
#'  * `hz`, total water depth (m)
#'  * `hc`, vertical water column occupied by the canopy (m)
#'  * `d_top`, depth of the top of the canopy beneath the water surface (m)
#' @param spec_params vector of named numbers. Must include:
#'  * `SA_WW`, conversion of wet weight to surface area
#' @param constants vector of named numbers defining extra constants for the attenuation submodel. Must include all constants, which have the following default values:
#'  * `s` = 0.0045
#'  * `gam` = 1.13
#'  * `a2` = 0.2^2
#'  * `Cb` = 0.0025
#'
#' @return the relative water attenuation coefficient (u_c)
#' @export
#' @seealso [height()], [u_b()], [C_t()]
#' 
u_c <- function(U0, macro_state, site_params, spec_params, 
  constants = c(s = 0.0045, gam = 1.13, a2 = 0.2^2, Cb = 0.0025)){
  
  if (missing(spec_params) | is.na(spec_params["SA_WW"])) {
    SA_WW <- 0.5 * (0.0306/2) # default is based on Macrocystis pyrifera
  } else {
    SA_WW <- spec_params["SA_WW"]
  }
  
  D <- SA_WW * min(macro_state['hm']/abs(site_params['hc']), 1) * macro_state['biomass']
  Kd <- 0.5 * abs(site_params['hz']) * D * constants['s'] * U0^(constants['gam'] - 2)
  Hc <- (abs(site_params['d_top']) + abs(site_params['hc'])) / abs(site_params['hz'])
  
  drag_test <- unname(suppressWarnings(
    sqrt(Kd * (1 - Hc) * Hc * (constants['Cb'] * Hc + constants['a2']) - constants['a2'] * constants['Cb'] * Hc)
    ))
    
  u_c <- if (is.na(drag_test)) {1} else {
    (-constants['a2'] - constants['Cb'] * Hc ^ 2 + (1 - Hc) * drag_test) / (Kd * Hc * (1 - Hc) ^ 3 - constants['a2'] - constants['Cb'] * Hc ^ 3)
  }
  return(unname(u_c))
}

#' Relative water attenuation beneath canopy
#'
#' @description
#' Calculates the relative water attenuation beneath the canopy, based on Plew, D. R. (2011). "Depth-Averaged Drag Coefficient for Modeling Flow through Suspended Canopies". Journal of Hydraulic Engineering, 137(2), 234–247. This function is not actually used within the main function `grow_macroalgae()`.
#' 
#' @inheritParams u_c
#'
#' @return a scalar of relative water attenuation beneath canopy
#' @export
#' 
#' @examples examples
u_b <- function(U0, macro_state, SA_WW = 0.5 * (0.0306/2), site_params, 
  constants = c(s = 0.0045, gam = 1.13, a2 = 0.2^2, Cb = 0.0025)){
  uc <- uc(U0, macro_state, SA_WW, site_params, constants)
  Hc <- (site_params['d_top'] + site_params['hc']) / site_params['hz']
  u_b <- (1 - u_c * Hc) / (1 - Hc)
  return(unname(u_b))
}

#' Total drag coefficient
#'
#' @description
#' Calculates the total drag coefficient, based on Plew, D. R. (2011). "Depth-Averaged Drag Coefficient for Modeling Flow through Suspended Canopies". Journal of Hydraulic Engineering, 137(2), 234–247. This function is not actually used within the main function `grow_macroalgae()`.
#'
#' @inheritParams u_c
#' 
#' @return The total drag coefficient C_t
#' 
C_t <- function(U0, macro_state, site_params, spec_params, 
  constants = c(s = 0.0045, gam = 1.13, a2 = 0.2^2, Cb = 0.0025)) {
  
  if (missing(spec_params) | is.na(spec_params["SA_WW"])) {
    SA_WW <- 0.5 * (0.0306/2) # default is based on Macrocystis pyrifera
  } else {
    SA_WW <- spec_params["SA_WW"]
  }
  D <- SA_WW * min(macro_state['hm']/abs(site_params['hc']), 1) * macro_state['biomass']
  Kd <- 0.5 * abs(site_params['hz']) * D * constants['s'] * U0^(constants['gam'] - 2)
  Hc <- (abs(site_params['d_top']) + abs(site_params['hc'])) / abs(site_params['hz'])
  u_c <- u_c(U0, macro_state, site_params, spec_params, constants)

  C_t <- (Kd * Hc * u_c ^ 2 + constants['Cb'] * u_b^2)
  return(unname(C_t))
}

# ---- Source: C:/Users/treimer/Documents/R-temp-files/macrogrow/R/errors.R ----
#' @title Custom error for missing parameters
#'
#' @param param missing parameter
#' @param place place where parameter should be defined
#' 
#' @keywords internal
#'
abort_missing_parameter <- function(param, place) {
  if (is.null(place)) {
    msg <- glue::glue("`{param}` must be defined")
  } else {
    msg <- glue::glue("`{param}` must be defined in {place}")
  }
  
  rlang::abort(
    class = "error_missing_parameter", 
        message = msg, 
        param = param, 
        place = place
  )
}


# ---- Source: C:/Users/treimer/Documents/R-temp-files/macrogrow/R/get_limiting.R ----
#' Determine limiting factor in growth
#' 
#' @description
#' `r lifecycle::badge("experimental")`
#' For each timestep in the growth output matrix, determine the factor(s) most likely to be limiting growth. This function is unfinished and subject to change.
#'
#' @param output matrix, direct output of 'grow_macroalgae' function (can be sparse)
#' @param spec_params a named vector of site-specific parameters - same as input to 'grow_macroalgae' function
#'
#' @return vector of limiting factors at each timestep
#' @export 
#' 
#' @details
#' This function gives an idea of what factor(s) are most likely to be limiting growth at each timestep. The choices are listed in order below. For each limitation, it can be assumed that the factors above it are within "ideal" ranges 
#' 1. `T_lim`: The macroalgae is limited by temperature
#' 2. `I_lim`: The macroalgae is limited by light availability
#' 3. `S_lim`: The macroalgae is limited by salinity
#' 4. `Q_lim`: The macroalgae does not have enough stored nutrients to fix into biomass at maximum rate allowed by temperature, light and salinity
#' 5. `Ns_to_Nf`: The macroalgae does not have enough stored nutrients to fix into biomass at maximum rate allowed by temperature, light and salinity
#' 6. `Nf_loss`: More Nf was lost than was fixed
#' 7. `conc_nitrate`/`conc_ammonium`/`conc_other`: There were not enough ambient nutrients available to replace fixed Ns
#' 8. `Ns_fixed_not_replaced`: the Ns fixed last timestep was not replaced this timestep (Q_rel is declining)
#' 9. `Ns_loss_not_replaced`: the Ns lost last timestep was not replaced this timestep (Q_rel is declining)
#'
#' @seealso [grow_macroalgae()]
get_limiting <- function(output, spec_params) {
  
  limiting <- rep(NA, nrow(output))
  first_lim <- second_lim <- NA

  for (i in 1:nrow(output)) {
    # growth_rate[i]  <- unname(spec_params['mu'] * min(T_lim[i], I_lim[i], S_lim[i]) * Q_lim[i])
    growth_factors <- c(output[i,'T_lim'], output[i,'I_lim'], output[i,'S_lim'], output[i,'Q_lim'])
    
    if (output[i,'growth_rate']/spec_params['mu'] == 1) {
      limiting[i] <- "No_limit"
    } else if (any(growth_factors == 0)) {
      zero_lim <- which(growth_factors == 0)
      zero_lim <- names(zero_lim)
      limiting[i] <- paste(zero_lim, collapse = " & ")
    } else {
      # 1. Temperature, light, salinity, or nutrients
      first_lim <- min(output[i,'T_lim'], output[i,'I_lim'], output[i,'S_lim'], output[i,'Q_lim'])
      first_lim <- which(c(output[i,'T_lim'], output[i,'I_lim'], output[i,'S_lim'], output[i,'Q_lim']) == first_lim)
      first_lim <- names(first_lim)
      
      # 2. Internal nutrients (ie there's not enough Ns to fix to Nf)
      second_lim_1 <- ifelse(output[i,'Ns_to_Nf'] > output[i, 'Ns'], "Ns_to_Nf", NA)

      # 3. Loss (ie more Nf is lost than fixed)
      second_lim_2 <- ifelse(output[i,'Nf_loss'] > output[i, 'Ns_to_Nf'], "Nf_loss", NA)
      
      # 4. Nitrogen uptake (ie not enough nitrogen was available to be taken up last step)
      check_conc <- c(
        output[i-1,'conc_nitrate'] - output[i-1,'up_Ni'],
        output[i-1,'conc_ammonium'] - output[i-1,'up_Am'],
        output[i-1,'conc_other'] - output[i-1,'up_Ot']
      )
      second_lim_3 <- ifelse(any(check_conc < 0, na.rm = T), paste(names(check_conc)[which(check_conc < 0)], collapse = " & "), NA)
      
      # 5. Nitrogen uptake (ie not enough nitrogen was taken up last step to replace fixed Ns)
      uptake <- sum(c(output[i-1,'up_Ni'], output[i-1,'up_Am'], output[i-1,'up_Ot']), na.rm = T)
      second_lim_4 <- ifelse(output[i,'Ns_to_Nf'] > uptake, "Ns_fixed_not_replaced", NA)

      # 6. Nitrogen uptake (ie not enough nitrogen was taken up last step to replace lost Ns)
      second_lim_5 <- ifelse(output[i,'Ns_loss'] > uptake, "Ns_loss_not_replaced", NA)
      
      second_lim <- c(second_lim_1, second_lim_2, second_lim_3, second_lim_4, second_lim_5)
      second_lim <- second_lim[!is.na(second_lim)]
      second_lim <- paste(second_lim, collapse = " & ")
        
      lims <- c(first_lim, second_lim)
      lims <- lims[!is.na(lims)]
      limiting[i] <- paste(lims, collapse = " & ")
      second_lim <- second_lim_1 <- second_lim_2 <- second_lim_3 <- second_lim_4 <- second_lim_5 <- NA
    }
  }
  return(limiting)
}

# ---- Source: C:/Users/treimer/Documents/R-temp-files/macrogrow/R/grow_macroalgae.R ----
#' Grow macroalgae
#' 
#' @description
#' Initiate and grow macroalgae.
#'
#' @param start integer, start of the growth period (day of at-sea deployment). Defaults to 1.
#' @param grow_days integer, number of day in growing period - if missing will take the length of the temperature vector
#' @param temperature a vector of daily temperatures (\eqn{^{\circ}}C)
#' @param salinity a vector of daily salt concentrations (g L\eqn{^{-1}})
#' @param light a vector of surface light (umol m\eqn{^{-2}} s\eqn{^{-1}})
#' @param kW a vector of light attenuation coefficients for open water (m\eqn{^{-1}})
#' @param velocity a vector of water velocities (m s\eqn{^{-1}})
#' @param nitrate a vector of nitrate concentrations (mg m\eqn{^{-3}})
#' @param ammonium a vector of ammonium concentrations (mg m\eqn{^{-3}})
#' @param ni_uptake shape for nitrate uptake, either "MM" for Michaelis-Menton saturation kinetics or "linear" for linear uptake
#' @param am_uptake shape for ammonium uptake, either "MM" for Michaelis-Menton saturation kinetics or "linear" for linear uptake
#' @param site_params a named vector of species-specific parameters - see details
#' @param spec_params a named vector of site-specific parameters - see details
#' @param initials a named vector of the macroalgae starting conditions
#' @param sparse_output logical, whether to include input vectors and other non-essential information in final dataframe (default = TRUE)
#' @param other_constants a named vector of other constants for intermediate equations (see `u_c()`)
#'
#' @importFrom glue glue
#' @import rlang
#' 
#' @return matrix of outputs
#' @export 
#' 
#' @details
#' Example csv with all the spec_params & site_params required? 
#' - Note that the final growth dataframe is inclusive of the start and end date, so the environmental vectors must be the same
#'
#' @examples "see here" link?
#' @seealso [check_grow()]
grow_macroalgae <- function(
  t = 1:30,
  temperature,
  salinity,
  light,
  kW,
  velocity,
  nitrate,
  ammonium,
  ni_uptake = NA,
  am_uptake = NA,
  site_params,
  spec_params,
  initials,
  sparse_output = T,
  other_constants = c(s = 0.0045, gam = 1.13, a2 = 0.2^2, Cb = 0.0025)
) {
  
  # Placeholder vectors
  u_c <- I_top <- conc_nitrate <- conc_ammonium <- Nf <- Ns <- B_dw.mg <- B_ww.mg <- hm <- lambda <- lambda_0 <- Q_int <- Q_rel <- T_lim <- S_lim <- Q_lim <- I_lim <- growth_rate <- Ns_to_Nf <- Ns_loss <- Nf_loss <- up_Am <- up_Ni <- 
    # red_Am <- remin <- 
    as.numeric(rep(NA, length(t)))

  add_ammonium     <- ammonium
  add_nitrate      <- nitrate

  use_Slim  <- ifelse(any(is.na(salinity)), yes = F, no = T)
  use_Ilim  <- ifelse(any(is.na(c(light, kW))), yes = F, no = T)
  use_Uc    <- ifelse(any(is.na(velocity)), yes = F, no = T)

  # External starting state
  conc_ammonium[1] <- add_ammonium[1]
  conc_nitrate[1]  <- add_nitrate[1]
  
  # Macroalgae starting state
  if (is.na(initials['Q_int'])) {
    initials['Q_int'] <- Q_int(Nf = initials['Nf'], Q_rel = initials['Q_rel'], spec_params = spec_params)
  }
  Nf[1]              <- unname(initials['Nf'])
  Ns[1]              <- unname(Nf[1]*(initials['Q_int']/spec_params['Q_min'] - 1))

  # Main run, after starting state
  for (i in 1:length(t)) {
    # Macroalgae state at start of day
    Q_int[i]       <- Q_int(Nf = Nf[i], Ns = Ns[i], spec_params = spec_params)
    Q_rel[i]       <- Q_rel(Q_int = Q_int[i], spec_params = spec_params)
    B_dw.mg[i]     <- Nf_to_biomass(Nf = Nf[i], Ns = Ns[i], Q_rel = Q_rel[i], dry = T, spec_params = spec_params)
    B_ww.mg[i]     <- Nf_to_biomass(Nf = Nf[i], Ns = Ns[i], Q_rel = Q_rel[i], dry = F, spec_params = spec_params)
    hm[i]          <- height(Nf[i], spec_params)
    
    # Environmental state (incoming)
    if (use_Uc) {
      biom <- B_ww.mg[i] / 1000 # mg -> g
      u_c[i] <- u_c(
        U0 = velocity[i], # m/s
        macro_state = c(biomass = biom, hm = hm[i]),
        site_params = site_params,
        spec_params = spec_params,
        constants = other_constants
      )
      U_0 <- velocity[i] * 86400 # m/s -> m/d
      lambda[i]      <- (u_c[i] * U_0)/unname(site_params['farmA'] * site_params['hc']) 
      lambda_0[i]    <- U_0/unname(site_params['farmA'] * site_params['hc']) 
    } else {
      u_c[i] <- lambda[i] <- lambda_0[i] <- 1
    }
    
    # Nutrient delivery
    conc_ammonium[i]     <- add_ammonium[i] * lambda_0[i]/lambda[i]
    conc_nitrate[i]      <- add_nitrate[i]  * lambda_0[i]/lambda[i]
    
    # Nitrogen uptake
    up_Am[i]        <- min(
      # Concentration of ammonium available
      conc_ammonium[i],
      # Uptake rate in current state
      (1 - Q_rel[i]) * B_dw.mg[i]/1000 *
        suppressMessages(get_uptake(
          conc = conc_ammonium[i], 
          spec_params = spec_params,
          uptake_shape = am_uptake, 
          Nform_abbr = "am"))
    )
    
    up_Ni[i]        <- min(
      conc_nitrate[i],
      (1 - Q_rel[i]) * (B_dw.mg[i]/1000) * 
        suppressMessages(get_uptake(
          conc = conc_nitrate[i], 
          uptake_shape = ni_uptake, 
          Nform_abbr = "ni", 
          spec_params = spec_params))
    )
    
    # Environmental limitation on growth
    T_lim[i]       <- T_lim(Tc = temperature[i], spec_params = spec_params)
    I_top[i]       <- unname(light[i] * exp(-kW[i] * site_params['d_top']))
    I_lim[i]       <- ifelse(use_Ilim, I_lim(Nf[i], I_top[i], kW[i], spec_params, site_params), 1)
    S_lim[i]       <- ifelse(use_Slim, S_lim(salinity[i], spec_params), 1)
    Q_lim[i]       <- Q_lim(Nf[i], Ns[i], spec_params)
    growth_rate[i] <- unname(spec_params['mu'] * min(T_lim[i], I_lim[i], S_lim[i]) * Q_lim[i])
    
    # Nitrogen fixation
    Ns_to_Nf[i]     <- min(growth_rate[i] * Ns[i], Ns[i])
    
    # Biomass loss
    U_c            <- velocity[i] * u_c[i] # m/s
    D_m            <- loss(U0 = U_c, spec_params = spec_params)
    N_loss         <- biomass_to_Nf(biomass = (B_ww.mg[i] * D_m), Q_rel = Q_rel[i], spec_params = spec_params, dry = F)
    Ns_loss[i]     <- unname(N_loss["Nf"])
    Nf_loss[i]     <- unname(N_loss["Ns"])
    
    # If you're not at the final day, set up for next day
    if (i < length(t)) {
      # IF MACROALGAE DIES
      if (Nf[i] <= 0) break
      
      # Change in algae state
      Nf[i+1]      <- Nf[i] + Ns_to_Nf[i] - Nf_loss[i]
      Ns[i+1]      <- Ns[i] + up_Am[i] + up_Ni[i] - Ns_to_Nf[i] - Ns_loss[i]
    }
  # End of main model run
  }
  
  # Some quick renaming
  add_nitrate <- nitrate
  add_ammonium <- ammonium
  
  # Put all the data together for outputs 
  if (sparse_output == F) {
    df <- cbind(t, Nf, Ns, growth_rate, Ns_to_Nf, Ns_loss, Nf_loss, Q_int, Q_rel, Q_lim, B_dw.mg, B_ww.mg, hm, add_nitrate, conc_nitrate, up_Ni, add_ammonium, conc_ammonium, up_Am, temperature, T_lim, salinity, S_lim, light, I_top, I_lim, velocity, u_c, lambda)
  } else {
    df <- cbind(t, Nf, Ns, growth_rate, Ns_to_Nf, Ns_loss, Nf_loss, Q_int, Q_rel, Q_lim, B_dw.mg, B_ww.mg, hm, conc_nitrate, up_Ni, conc_ammonium, up_Am, T_lim, S_lim, I_top, I_lim, u_c)
  }
  return(df)
}



# ---- Source: C:/Users/treimer/Documents/R-temp-files/macrogrow/R/I_lim.R ----
#' Light limitation on growth
#' 
#' @description
#' Calculates the relative limitation on growth rate due to light availability, via:
#' \deqn{I_{lim} = \frac{e}{K \cdot d_{top}} \times \Biggl[e^{-\frac{I_z e^{-K \cdot d_{top}}}{I_o}} - e^{-\frac{I_z}{I_o}} \Biggr]}
#' where \eqn{I_{z}=I e^{-k_W \cdot d_{top}}} is the irradiance at the cultivation depth 
#' and \eqn{K=k_{m}+kW} is the total attenuation coefficient. 
#' 
#' \eqn{k_{m}} is the additional attenuation coefficient from macroalgae biomass, calculated as:
#' \deqn{k_{m} = a_{cs} \times N_f \times \text{max} \Biggl( \frac{h_{m}}{d_{top}}, 1 \Biggr) \times \frac{1}{\text{min}(h_{m}, d_{top})}}
#' where \eqn{h_m} is the macroalgae height. 
#' 
#' @inheritParams Q_int
#' @param I the surface irradiance, PAR (\eqn{\mu}mol photons m\eqn{^{-2}} s\eqn{^{-1}})
#' @param kW the light attenuation coefficient for open water (m\eqn{^{-1}})
#' @param spec_params a vector of named numbers. Must include:
#' * `a_cs`, the carbon-specific self-shading constant
#' * `I_o`, the light saturation parameter
#' * `h_a`, `h_b` and `h_c`, parameters governing height change with `N_f`
#' * `h_max`, maximum species height
#' @param site_params A vector of named numbers. Must include:
#' * `d_top`, the below-surface  depth (m) of the top of the macroalgae culture
#'
#' @return a scalar of relative light limitation on growth (between 0 and 1)
#' @export
#'
#' @examples 
#' my_species <- c(I_o = 200, a_cs = 0.001, h_a = 750, h_b = 0.5, h_c = 0.01, h_max = 1)
#' my_site <- c(d_top = 1)
#' I <- seq(0, 2000, 10)
#' kW <- rnorm(length(I), 0.05, 0.01)
#' 
#' # Use purrr::map2 to calculate light limitation over a range of light and attenuation values, but fixed Nf
#' \dontrun{
#' df <- purrr::map2_dfr(I, kW, function(Iz, kWz) {
#'   data.frame(
#'     I = Iz,
#'     kW = kWz,
#'     lim = I_lim(Nf = 1000, I = Iz, kW = kWz, spec_params = my_species, site_params = my_site)
#'   )
#' })
#' }
#' 
#' @seealso [height()]
#' 
I_lim <- function(Nf, I, kW, spec_params, site_params) {
  # Check that required parameters are supplied
  if (is.na(site_params['d_top'])) {abort_missing_parameter(param = "d_top", place = "site_params")}
  if (is.na(spec_params['a_cs'])) {abort_missing_parameter(param = "a_cs", place = "spec_params")}
  if (is.na(spec_params['I_o'])) {abort_missing_parameter(param = "I_o", place = "spec_params")}
  
  I_top <- I * exp(-(kW*site_params['d_top']))
  
  h_m <- height(Nf, spec_params)
  k_ma <-  Nf * h_m * spec_params['a_cs'] * pmax(h_m/site_params['d_top'], 1) * 1/(pmin(h_m, site_params['d_top']))
  K <- k_ma + kW
  
  Ilim <- exp(1)/(K*h_m) *
    (
      exp(-(I_top*exp(-K*h_m))/spec_params['I_o']) - exp(-I_top/spec_params['I_o'])
      )
  
  return(unname(Ilim))
}

# ---- Source: C:/Users/treimer/Documents/R-temp-files/macrogrow/R/loss.R ----
#' Loss through fragmentation
#' 
#' @description
#' The macroalgae fragmentation rate (\eqn{D}) is calculated as:
#' \deqn{D = D_{ve} \times U_0 + D_{tu} + D_{m}} 
#' where \eqn{D_{ve}} is the linear coefficient of loss with laminar velocity (\eqn{U_0}), 
#' \eqn{D_{tu}} is the loss rate at the specified level of turbulence (static, low, medium, or high), and
#' \eqn{D_{m}} is the base loss rate.
#'
#' @inheritParams u_c 
#' @param turbulence optional level of turbulence. Options are "none"/"static", "low", "medium"/"mid", and "high".
#' @param spec_params a vector of named numbers. All default to 0 if not supplied. Can include: 
#'  * `D_m`, a constant (base) loss rate
#'  * `D_ve`, the linear coefficient of loss with laminar velocity
#'  * `D_st`, the loss rate in static water (turbulence = 0, "none" or "static")
#'  * `D_lo`, the loss rate in turbulent water (turbulence = 1 or "low")
#'  * `D_mi`, the loss rate in turbulent water (turbulence = 2 or "medium")
#'  * `D_hi`, the loss rate in turbulent water (turbulence = 3 or "high")
#'
#' @details
#' Note that parameter `U0` has no default. 
#' If one is not supplied, loss due to velocity will default to 0 regardless of the value of `D_ve` supplied. 
#' 
#' @return a scalar for macroalgae loss (d-1)
#' @export
#'
#' @examples examples
#' 
loss <- function(U0, turbulence = NA, spec_params) {

  D_m <- ifelse(is.na(spec_params['D_m']), 0, spec_params['D_m'])  
  D_ve <- ifelse(is.na(spec_params['D_ve']), 0, spec_params['D_ve'])

  D_turbulence <- if(is.na(turbulence) | turbulence == "static" | turbulence == "none") {
    ifelse(is.na(spec_params['D_st']), 0, spec_params['D_st'])
  } else if (turbulence == "low") {
    ifelse(is.na(spec_params['D_lo']), 0, spec_params['D_lo'])
  } else if (turbulence == "medium" | turbulence == "mid") {
    ifelse(is.na(spec_params['D_mi']), 0, spec_params['D_mi'])
  } else if (turbulence == "high") {
    ifelse(is.na(spec_params['D_hi']), 0, spec_params['D_hi'])
  }

  # Actual loss calculation
  D_m <- U0 * D_ve + D_turbulence + D_m
  
  return(unname(D_m))
}

# ---- Source: C:/Users/treimer/Documents/R-temp-files/macrogrow/R/N_uptake_functions.R ----
#' Get the correct uptake rate
#' 
#' @description
#' Gives the substrate uptake rate based on the shape and species-specific parameters supplied. 
#' If no shape is supplied, attempts to infer the correct shape from supplied parameters. 
#'
#' @param conc substrate concentration
#' @param uptake_shape kinetic shape for substrate uptake. One of "linear" or "Michaelis-Menton" (or "MM"). Defaults to Michaelis-Menton. 
#' @param spec_params a vector of named numbers. Must include:
#' * `M` and `C` for linear uptake, OR
#' * `V` and `K` for Michaelis-Menton uptake
#' * `uptake_shape` of "linear" or "MM" (or "Michaelis-Menton"). Not required if only one set of the above parameters are provided.
#' @param Nform_abbr the abbreviation used in spec_params for the relevant nitrogen form. See details.
#'
#' @details
#' `Nform_abbr` allows the function to link species parameters with concentration. 
#' E.g. if `Nform_abbr` = "amm" (for ammonium) the function will look for `M_amm` and `C_amm` or `V_amm` and `K_amm` in spec_params and will ignore other uptake parameters which may be included for other substrates. 
#' 
#' @export
#' @seealso [lin_uptake()] [MM_uptake()]
#' 
#' @examples
#' Ni <- seq(0, 300, 10)
#' Am <- seq(0, 1000, 20)
#' my_species <- c(V_ni = 4.5, K_ni = 300, V_am = 5.5, K_am = 150, M_am = 1.5, C_am = 10)
#' 
#' # Get uptake for a vector of concentrations
#' sapply(Am, function(conc) {get_uptake(conc = conc, uptake_shape = "linear", Nform_abbr = "am", spec_params = my_species)})
#' 
#' # Uptake shape does not need to be specified if only one set of parameters is provided
#' sapply(Ni, function(conc) {get_uptake(conc = conc, uptake_shape = NA, Nform_abbr = "ni", spec_params = my_species)})
#' 
#' # Create a custom compound to uptake
#' my_species <- c(V_urea = 1.5, K_urea = 10)
#' Urea <- seq(0, 25, 0.5)
#' sapply(Urea, function(conc) {get_uptake(conc = conc, Nform_abbr = "urea", spec_params = my_species)})
#' 
get_uptake <- function(conc, uptake_shape = NA, Nform_abbr, spec_params) {
  
  if (is.na(conc)) {
    rlang::abort("Concentration passed to get_uptake cannot be NA", class = "error_bad_parameter")
  }
  
  M <- spec_params[paste0("M", "_", Nform_abbr)]
  C <- spec_params[paste0("C", "_", Nform_abbr)]
  V <- spec_params[paste0("V", "_", Nform_abbr)]
  K <- spec_params[paste0("K", "_", Nform_abbr)]
  
  if (is.na(uptake_shape)) {
    if (!is.na(V) & !is.na(K)) {
      up <- MM_uptake(conc = conc, V = V, K = K)
    } else if (!is.na(M) & !is.na(C)) {
      up <- lin_uptake(conc = conc, M = M, C = C)
    } else {
      rlang::abort(glue::glue("Error! You haven't provided the correct uptake parameters for '{form}'!", form = Nform_abbr))
    }
  } else if (uptake_shape == "linear") {
    if(is.na(M)) {
      abort_missing_parameter(param = paste0("M", "_", Nform_abbr), place = "spec_params for linear uptake")
    } else if(is.na(C)) {
      abort_missing_parameter(param = paste0("C", "_", Nform_abbr), place = "spec_params for linear uptake")
    } else {
      up <- lin_uptake(conc = conc, M = M, C = C)
    }
  } else if (uptake_shape == "MM" | uptake_shape == "Michaelis-Menton") {
    if(is.na(V)) {
      abort_missing_parameter(param = paste0("V", "_", Nform_abbr), place = "spec_params for Michaelis-Menton uptake")
    } else if(is.na(K)) {
      abort_missing_parameter(param = paste0("K", "_", Nform_abbr), place = "spec_params for linear uptake")
    } else {
      up <- MM_uptake(conc = conc, V = V, K = K)
    }
  } else {
    rlang::abort(glue::glue("Uptake shape `{uptake_shape}` is not recognised in FORT KICKASS"), class = "error_bad_parameter")
  }
  return(up)
}

#' Uptake rate (linear)
#' 
#' @inheritParams MM_uptake
#' @param M the slope of N uptake with increasing substrate concentration 
#' @param C the intercept
#'
#' @return the rate of uptake at the specified external concentration
#' @export
#'
#' @examples examples
lin_uptake <- function(conc, M, C) {
  uprate <- M * conc + C
  return(unname(uprate))
}

#' Uptake rate (Michaelis-Menton)
#'
#' @param conc external substrate concentration
#' @param V the maximum uptake rate \eqn{V_{max}}
#' @param K the half-saturation constant \eqn{K_{c}}
#'
#' @return the rate of uptake at the specified external concentration
#' @export
#'
MM_uptake <- function(conc, V, K) {
  if (missing(V) & missing(K)) {abort_missing_parameter(param = "V and K", place = "spec_params. Did you mean to use lin_uptake() instead?")}
  if (missing(V)) {abort_missing_parameter(param = "V", place = "spec_params and passed to function as 'V'")}
  if (missing(K)) {abort_missing_parameter(param = "K", place = "spec_params and passed to function as 'K'")}
  
  uprate <- (V * conc / (K + conc))
  return(unname(uprate))
}

# ---- Source: C:/Users/treimer/Documents/R-temp-files/macrogrow/R/Q_int.R ----
#' Non-dimensionalised internal nutrient quotient
#' 
#' @description Calculates the internal nutrient quotient \eqn{Q} via:
#' \deqn{\begin{array}[ccc] 
#' Q &=& Q_{min} \left(1 + \frac{N_s}{N_f}\right)
#' \end{array}}
#'
#' @param Nf Fixed nitrogen (mg m\eqn{^{-3}})
#' @param Ns Stored nitrogen (mg m\eqn{^{-3}})
#' @param spec_params A vector of named numbers. Must include:
#'  * `Q_min`, the minimum internal nutrient quotient
#'
#' @return the non-dimensionalised internal nutrient quotient
#' @export
#' @seealso [Q_rel()]
#'
#' @examples examples
Q_int <- function(Nf = NULL, Ns = NULL, Q_rel = NULL, spec_params) {
  if (!is.null(Q_rel)) {
    Q_int <- spec_params['Q_max'] - (1 - Q_rel)*(spec_params['Q_max'] - spec_params['Q_min'])
  } else if (!is.null(Nf) & !is.null(Ns)) {
    Q_int <- spec_params['Q_min'] * (1 + Ns/Nf)
  }
  if (Q_int > spec_params['Q_max']) {Q_int <- spec_params['Q_max']}
  if (Q_int < spec_params['Q_min']) {Q_int <- spec_params['Q_min']}

  return(unname(Q_int))
}



# ---- Source: C:/Users/treimer/Documents/R-temp-files/macrogrow/R/Q_lim.R ----
#' Q_int limitation on growth
#'
#' @description Calculates the internal nutrient quotient \eqn{Q_{int}} and its relative effect on growth rate \eqn{Q_{lim}} via:
#' \deqn{\begin{array}[ccc] 
#' Q_{int} &=& Q_{min} \left(1 + \frac{N_s}{N_f}\right) \\
#' Q_{lim} &=& \frac{Q_{int} - Q_{min}}{Q_{int} - K_c}
#' \end{array}}
#' 
#' @inheritParams Q_int
#' @param spec_params a vector of named numbers. Must include:
#'  * `Q_min`, the minimum nutrient quotient
#'  * `Q_max`, the maximum nutrient quotient
#'  * `K_c`, the half-saturation growth constant
#' 
#' @return a scalar of relative limitation from internal nutrient reserves on growth (between 0 and 1)
#' @export
#' 
#' @examples examples
#' @seealso [Q_int()], [Q_rel()]
#' 
Q_lim <- function(Nf, Ns, spec_params) {
  Q_int <- Q_int(Nf = Nf, Ns = Ns, spec_params = spec_params)
  if (spec_params['K_c'] >= spec_params['Q_min']) {rlang::abort("K_c must be less than Q_min", class = "error_bad_parameter")}
  if (Q_int < spec_params['Q_min']) {Q_int <- spec_params['Q_min']} 
  if (Q_int > spec_params['Q_max']) {Q_int <- spec_params['Q_max']}
  
  Q_lim <- (Q_int - spec_params['Q_min'])/(Q_int - spec_params['K_c'])
  return(unname(Q_lim))
}

# ---- Source: C:/Users/treimer/Documents/R-temp-files/macrogrow/R/Q_rel.R ----
#' Relative nutrient quotient
#'
#' @param Q_int Non-dimensionalised internal nutrient quotient
#' @param spec_params A vector of named numbers. Must include:
#'  * `Q_min`, the minimum internal nutrient quotient
#'  * `Q_max`, the maximum internal nutrient quotient
#'
#' @return the relative (0-1) internal nutrient quotient
#' @export
#' @seealso [Q_int()]
#'
Q_rel <- function(Nf = NULL, Ns = NULL, Q_int = NULL, spec_params) {
  if (!is.null(Nf) & !is.null(Ns)) {
    Q_int <- Q_int(Nf = Nf, Ns = Ns, spec_params = spec_params)
  }
  if (!is.null(Q_int)) {
    Q_rel <- 1 - (spec_params['Q_max'] - Q_int)/(spec_params['Q_max'] - spec_params['Q_min'])
  }
  return(unname(Q_rel))
}


# ---- Source: C:/Users/treimer/Documents/R-temp-files/macrogrow/R/S_lim.R ----
#' Salinity limitation on growth
#' 
#' Given species parameters, returns the relative limitation on growth rate according to a CTMI curve:
#' \deqn{
#'      S_{lim} &=& \frac{(Sal-S_{max})(Sal-S_{min})^2}{(S_{opt}-S_{min})[(S_{opt}-S_{min})(Sal-S_{opt})-(S_{opt}-S_{max})(S_{opt}+S_{min}-2Sal)]}
#' }
#' The CTMI curve was formulated for temperature.
#' 
#' @param Sal salinity to evaluate (g/L)
#' @param spec_params a vector of named numbers. Must include:
#'  * `S_opt` the optimum salinity for macroalgae growth
#'  * `S_min` the minimum salinity for macroalgae growth (when `Sal` <= `S_min`, growth = 0)
#'  * `S_max` the maximum salinity for macroalgae growth (when `Sal` >= `S_max`, growth = 0)
#'
#' @return a scalar of relative salinity limitation on growth (between 0 and 1)
#' @export
#'
#' @examples 
#' my_seaweed <- c(S_opt = 20, S_min = 5, S_max = 30)
#' 
#' S_lim(Sal = 22, spec_params = my_seaweed)
#' 
#' S_range <- 1:30
#' sapply(S_range, S_lim, spec_params = my_seaweed)
#' 
S_lim <- function(Sal, spec_params){
  
  # Check that required parameters are supplied
  if (is.na(spec_params['S_opt'])) {abort_missing_parameter(param = "S_opt", place = "spec_params")}
  if (is.na(spec_params['S_min'])) {abort_missing_parameter(param = "S_min", place = "spec_params")}
  if (is.na(spec_params['S_max'])) {abort_missing_parameter(param = "S_max", place = "spec_params")}
  
  if (spec_params['S_opt'] < spec_params['S_min']) {rlang::abort("error_bad_parameter", message = "Error: Minimum salinity is higher than optimum salinity")}
  if (spec_params['S_opt'] > spec_params['S_max']) {rlang::abort("error_bad_parameter", message = "Error: maximum salinity is lower than optimum salinity")} 
  if (spec_params['S_opt']-spec_params['S_min'] <= spec_params['S_max']-spec_params['S_opt']) {rlang::abort("error_bad_parameter", message = "Species CTMI function not valid! Must satisfy S_opt-S_min > S_max-S_opt")}
  
  Slim <- if (is.na(Sal)) {
    1
  } else if (Sal >= spec_params['S_max']) {
    0
  } else if (Sal <= spec_params['S_min']) {
    0
  } else {
    ((Sal - spec_params['S_max']) * (Sal - spec_params['S_min'])^2) / ((spec_params['S_opt'] - spec_params['S_min']) * ((spec_params['S_opt'] - spec_params['S_min']) * (Sal - spec_params['S_opt']) - (spec_params['S_opt'] - spec_params['S_max']) * (spec_params['S_opt'] + spec_params['S_min'] - 2 * Sal)))
  }
  return(unname(Slim))
}

# ---- Source: C:/Users/treimer/Documents/R-temp-files/macrogrow/R/Secchi_Kd.R ----
#' @title Get light attenuation parameter Kd from Secchi depth
#'
#' @description
#' Calculates the light attenuation coefficient in water via:
#' \deqn{\begin{array}[ccc] 
#' K_d(PAR) &=& 1.584 × SDD^{-0.894}
#' \end{array}}
#' where SDD is the Secchi disk depth. Equation is adapted from Zhang, Y., Liu, X., Yin, Y., Wang, M., & Qin, B. (2012). Predicting the light attenuation coefficient through Secchi disk depth and beam attenuation coefficient in a large, shallow, freshwater lake. Hydrobiologia, 693(1), 29–37. https://doi.org/10.1007/s10750-012-1084-2.
#' 
#' @param SDD Secchi disk depth (m)
#'
#' @return attenuation parameter for PAR
#' @export
#'
Secchi_to_Kd <- function(SDD) {
  1.584 * SDD^(-0.894)
}

#' @title Get Secchi depth from light attenuation parameter Kd
#'
#' @description
#' Calculates the light attenuation coefficient in water via:
#' \deqn{\begin{array}[ccc] 
#' K_d(PAR) &=& 1.584 × SDD^{-0.894}
#' \end{array}}
#' where SDD is the Secchi disk depth. Equation is adapted from Zhang, Y., Liu, X., Yin, Y., Wang, M., & Qin, B. (2012). Predicting the light attenuation coefficient through Secchi disk depth and beam attenuation coefficient in a large, shallow, freshwater lake. Hydrobiologia, 693(1), 29–37. https://doi.org/10.1007/s10750-012-1084-2.
#' 
#' @param Kd attenuation parameter for PAR
#'
#' @return Secchi disk depth (m)
#
# ---- Source: C:/Users/treimer/Documents/R-temp-files/macrogrow/R/solarMJ2ppfd.R ----
#' Solar Radiation to PPFD
#'
#' The following function and documentation was copied verbatim from file modules/data.atmosphere/R/metutils.R in https://github.com/PecanProject/pecan/.
#' 
#' There is no easy straight way to convert MJ/m2 to mu mol photons / m2 / s (PAR).
#' Note: 1 Watt = 1J/s
#' The above conversion is based on the following reasoning
#' 0.12 is about how much of the total radiation is expected to ocurr during the hour of maximum insolation (it is a guesstimate)
#' 2.07 is a coefficient which converts from MJ to mol photons (it is approximate and it is taken from ...
#' Campbell and Norman (1998). Introduction to Environmental Biophysics. pg 151 'the energy content of solar radiation in the PAR waveband is 2.35 x 10^5 J/mol'
#' See also the chapter radiation basics (10)
#' Here the input is the total solar radiation so to obtain in the PAR spectrum need to multiply by 0.486
#' This last value 0.486 is based on the approximation that PAR is 0.45-0.50 of the total radiation
#' This means that 1e6 / (2.35e6) * 0.486 = 2.07
#' 1e6 converts from mol to mu mol
#' 1/3600 divides the values in hours to seconds
#'
#' @title MJ to PPFD
#' @author Fernando Miguez
#' @author David LeBauer
#' @param solarMJ MJ per day
#' 
#' @return PPFD umol m-2 s-1
#' @export
#' 
solarMJ2ppfd <- function(solarMJ) {
  ppfd <- (0.12 * solarMJ) * 2.07 * 1e+06 / 3600
  return(ppfd)
} 

# ---- Source: C:/Users/treimer/Documents/R-temp-files/macrogrow/R/T_lim.R ----
#' Temperature limitation on growth
#' 
#' Given species parameters, returns the relative limitation on growth rate according to a CTMI curve:
#' \deqn{\begin{array}[ccc] 
#'      T_{lim} &=& \frac{(T_c-T_{max})(T_c-T_{min})^2}{(T_{opt}-T_{min})[(T_{opt}-T_{min})(T_c-T_{opt})-(T_{opt}-T_{max})(T_{opt}+T_{min}-2T_c)]}
#' \end{array}}
#' 
#' @param Tc temperature to evaluate
#' @param spec_params a vector of named numbers. Must include:
#'  * `T_opt` the optimum temperature for macroalgae growth
#'  * `T_min` the minimum temperature for macroalgae growth (when `T_c` < `T_min`, growth = 0)
#'  * `T_max` the maximum temperature for macroalgae growth (when `T_c` > `T_max`, growth = 0)
#'
#' @return a scalar of relative temperature limitation on growth (between 0 and 1)
#' @export
#'
#' @examples 
#' my_seaweed <- c(T_opt = 20, T_min = 5, T_max = 30)
#' 
#' T_lim(Tc = 22, spec_params = my_seaweed)
#' 
#' T_range <- 1:30
#' sapply(T_range, T_lim, spec_params = my_seaweed)
#' 
T_lim <- function(Tc, spec_params){
  
  # Check that required parameters are supplied
  if (is.na(spec_params['T_opt'])) {abort_missing_parameter(param = "T_opt", place = "spec_params")}
  if (is.na(spec_params['T_min'])) {abort_missing_parameter(param = "T_min", place = "spec_params")}
  if (is.na(spec_params['T_max'])) {abort_missing_parameter(param = "T_max", place = "spec_params")}
  
  if (spec_params['T_opt'] < spec_params['T_min']) {
    rlang::abort("error_bad_parameter", message = "Minimum temperature is higher than optimum temperature")
  }
  if (spec_params['T_opt'] > spec_params['T_max']) {
    rlang::abort("error_bad_parameter", message = "Error: maximum temperature is lower than optimum temperature")
  } 
  if (spec_params['T_opt']-spec_params['T_min'] <= spec_params['T_max']-spec_params['T_opt']) {
    rlang::abort("error_bad_parameter", message = "Species CTMI function not valid! Must satisfy T_opt-T_min > T_max-T_opt")
  }
  
  if (Tc >= spec_params['T_max']) {
    Tlim <- 0
  } else if (Tc <= spec_params['T_min']) {
    Tlim <- 0
  } else {
    Tlim <- ((Tc - spec_params['T_max'])*(Tc - spec_params['T_min'])^2)/((spec_params['T_opt'] - spec_params['T_min'])*((spec_params['T_opt'] - spec_params['T_min'])*(Tc - spec_params['T_opt']) - (spec_params['T_opt'] - spec_params['T_max'])*(spec_params['T_opt'] + spec_params['T_min'] - 2*Tc)))
  }
  return(unname(Tlim))
}
