functions {
  real Phi_p(real eta, real p) {
    real u;
    real g;

    u = abs(eta)^p / p;
    g = gamma_p(1.0 / p, u);

    if (eta >= 0)
      return 0.5 + 0.5 * g;
    else
      return 0.5 - 0.5 * g;
  }
}

data {
  int<lower=1> N;
  int<lower=1> d;
  matrix[N, d] X;
  array[N] int<lower=0, upper=1> y;
  
  real<lower=0> p_lower;
  real<lower=p_lower> p_upper;
  
  int<lower=1, upper=2> alpha_prior_type;
  real alpha_prior_loc;
  real<lower=0> alpha_prior_scale;
  real<lower=0> alpha_prior_df;
  
  int<lower=1, upper=2> beta_prior_type;
  real beta_prior_loc;
  real<lower=0> beta_prior_scale;
  real<lower=0> beta_prior_df;
}

parameters {
  real alpha; // intercept
  vector[d] beta; // coefficients
  real<lower=p_lower, upper=p_upper> p;
}

model {
  if (alpha_prior_type == 1)
    alpha ~ normal(alpha_prior_loc, alpha_prior_scale);
  else if (alpha_prior_type == 2)
    alpha ~ student_t(alpha_prior_df, alpha_prior_loc, alpha_prior_scale);
    
  if (beta_prior_type == 1)
    beta ~ normal(beta_prior_loc, beta_prior_scale);
  else if (beta_prior_type == 2)
    beta ~ student_t(beta_prior_df, beta_prior_loc, beta_prior_scale);
    
  vector[N] eta = alpha + X * beta;
  
  for (i in 1:N) {
    real prob = Phi_p(eta[i], p);
    prob = fmin(1 - 1e-12, fmax(1e-12, prob));
    target += bernoulli_lpmf(y[i] | prob);
  }
}
