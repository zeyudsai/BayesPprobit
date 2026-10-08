# prior_normal()      # validates Normal parameters + creates classed object
# prior_student_t()   # validates Student-t parameters + creates classed object

# validate_prior()    # only checks that it is a BayesPprobit prior object
# prior_type_code()   # translates family to 1L / 2L
# stan_p_probit()     # performs validation checks, creates the Stan data,
                      # locates the Stan file, compiles it, and samples
################################################################################

prior_normal <- function(location = 0, scale = 2.5) {
  # validation #
  if (!is.numeric(location) || length(location) != 1 || !is.finite(location)) {
    stop("'location' must be a finite real number.")
  }
  
  if (!is.numeric(scale) || length(scale) != 1 ||
      !is.finite(scale) || scale <= 0) {
    stop("'scale' must be a positive finite real number.")
  }
  
  # makes a BayesPprobit_prior object #
  structure(
    list(
      family = "normal",
      location = location,
      scale = scale,
      df = 1
    ),
    class = "BayesPprobit_prior"
  )
}

prior_student_t <- function(df = 3, location = 0, scale = 2.5) {
  # validation #
  if (!is.numeric(location) || length(location) != 1 || !is.finite(location)) {
    stop("'location' must be a finite real number.")
  }
  
  if (!is.numeric(scale) || length(scale) != 1 ||
      !is.finite(scale) || scale <= 0) {
    stop("'scale' must be a positive finite real number.")
  }
  
  if (!is.numeric(df) || length(df) != 1 || !is.finite(df) || df <= 0) {
    stop("'df' is the prior's number of degrees of freedom and must therefore be
         a positive finite real number.")
  }
  
  # makes a BayesPprobit_prior object #
  structure(
    list(
      family = "student_t",
      location = location,
      scale = scale,
      df = df
    ),
    class = "BayesPprobit_prior"
  )
}

prior_type_code <- function(prior, name = "prior") {
  switch(
    prior$family,
    normal = 1L,
    student_t = 2L,
    stop(sprintf(
      "'%s' uses an unsupported prior family: %s. You can only choose a normal or
      a Student's t-distribution.",
      name,
      prior$family
      ))
  )
}

validate_prior <- function(prior, name) {
  if (!inherits(prior, "BayesPprobit_prior")) {
    stop(sprintf(
      "'%s' must be created with 'prior_normal()' or 'prior_student_t()'.",
      name
    ))
  }
}

stan_p_probit <- function(
    X,
    y,
    prior_intercept = prior_normal(0, 5),
    prior_beta = prior_normal(0, 2.5),
    p_bounds = c(1, 6),
    chains = 4,
    parallel_chains = chains,
    iter_warmup = 1000,
    iter_sampling = 1000,
    seed = NULL,
    max_treedepth = 12,
    ...
) {
  
  ######################## validation ##########################################
  
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("'X' is either not numeric or not a matrix. Please fix that.")
  }
  
  if (anyNA(X)) {
    stop("I'm sorry to inform you that 'X' contains missing values. I can't handle
         that.")
  }
  
  if (length(y) != nrow(X)) {
    stop("Obviously, the length of 'y' should equal the number of rows in 'X'.
         Right now, it does not.")
  }
  
  if (!is.numeric(y)) {
    stop("Your target variable is not even numeric. How is that supposed to work?")
  }
  
  if (anyNA(y)) {
    stop("Your target variable is missing some values. That's a serious problem.")
  }
  
  if(!all(y %in% c(0, 1))) {
  stop("'y' needs to be a binary variable, fam.")
  }
  
  if (!is.numeric(p_bounds) ||
      length(p_bounds) != 2 ||
      p_bounds[1] <= 0 ||
      p_bounds[1] >= p_bounds[2]) {
    stop("'p_bounds' must contain two numbers with 0 < lower < upper. But I'm
         sure you know this and this was just a typo.")
  }
  
  validate_prior(prior_intercept, "prior_intercept")
  validate_prior(prior_beta, "prior_beta")
  
  ########################## create the data used in the Stan sampler ########## 
  
  stan_data <- list(
    N = nrow(X),
    d = ncol(X),
    X = X,
    y = as.integer(y),
    
    p_lower = p_bounds[1],
    p_upper = p_bounds[2],
    
    alpha_prior_type = prior_type_code(prior_intercept, "prior_intercept"),
    alpha_prior_loc = prior_intercept$location,
    alpha_prior_scale = prior_intercept$scale,
    alpha_prior_df = prior_intercept$df,
    
    beta_prior_type = prior_type_code(prior_beta, "prior_beta"),
    beta_prior_loc = prior_beta$location,
    beta_prior_scale = prior_beta$scale,
    beta_prior_df = prior_beta$df
  )
  
  ################################# locate the Stan file #######################
  
  stan_file <- system.file(
    "stan",
    "p_gen_probit.stan",
    package = "BayesPprobit"
  )
  
  if (stan_file == "") {
    stop("It seems that the Stan model file could not be found. That's certainly
         troubling.")
  }
  
  ############################## compile the Stan file #########################
  
  mod <- cmdstanr::cmdstan_model(stan_file)
  
  ############################## sample ########################################
  
  mod$sample(
    data = stan_data,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    max_treedepth = max_treedepth,
    ...
  )

}