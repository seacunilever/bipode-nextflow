// Implementation of the BIFROST model using a double skew logistic distribution.
// Assumes column biases are present in the data.
// Author: Joe Reynolds
functions {

  /*
     Helper function for the calculation of the beta logistic
     log density

     @param x vector value in support of distribution
     @param mu vector mean of the distribution
     @param sigma vector standard deviation of the distribtion
     @param a real shape parameter controlling the left tail
     @param b real shape parameter controlling the right tail
     @return vector log density for each component of x
   */
  real beta_logistic_lpdf(vector x, vector mu, real sigma, real a, real b) {

    int n = dims(x)[1];
    real u = digamma(a) - digamma(b);
    real v = sqrt(trigamma(a) + trigamma(b));
    vector[n] y = u + v .* (x - mu) / sigma;
    vector[n] log_dens = - a * log1p_exp(- y) - b * log1p_exp(y) - log(sigma) - lbeta(a, b) + log(v);
    return sum(log_dens);

  }

  // Quantile function defined by a logistic distribution on the logit scale
  real logistic_quantile(real logit_q, real lower_limit, real upper_limit, real slope, real loc) {
    real x = lower_limit + (upper_limit - lower_limit) / (1 + exp(-slope * (logit_q - loc)));
    return x;
  }

}
data {

  int n_sample; // Number of observations
  int n_treatment_batch; // Number of exposure plates/biological replicates
  array[n_sample] int count; // Counts for each sample
  vector[n_sample] total_count; // Maximum count for sample

  int n_batch; // Number of batches
  array[n_sample] int batch_index; // Batch index for each count

  int n_conc; // Number of unique treatment concentrations
  vector[n_conc] conc; // Treatment concentration
  array[n_sample] int conc_index; // Index of concentration group for each count

  int n_low_count; // Number of counts <= 100
  array[n_low_count] int low_count_index; // Sample index for low counts

  int n_high_count; // Number of counts > 100
  array[n_high_count] int high_count_index; // Sample index for low counts

}
transformed data {

  vector[n_low_count] log_total_count = log(total_count[low_count_index]);
  vector[n_high_count] n = to_vector(total_count[high_count_index]);

  real min_conc = min(conc);
  real max_conc = max(conc);

  real theta_l = min_conc - 1;
  real theta_u = max_conc + 1;
  real theta_a = -2 * log(theta_u - theta_l - 1) / logit(0.025);
  real theta_b = -log(theta_u - theta_l - 1) / theta_a;

  real l19 = log(19);

}
parameters {

  // Median shifted log-odds of reading probe on each plate
  real mu_loc;
  real<lower=0> mu_scale;
  vector[n_batch] mu_raw; // Mean log-odds on each plate

  real<lower=0> sigma;

  // Tail behaviour
  real<lower=0> a; // Parameter controlling excess left skew
  real<lower=0> b; // Parameter controlling excess right skew

  // Latent log-odds estimates
  vector[n_sample] log_odds;

  // Treatment effect parameters
  vector[n_conc] treatment_response;
  real theta_raw; // 5% effect concentration (depth)
  real<lower=0> delta; // Difference between 5% and 50% effect concentrations
  real<lower=0> gamma; // Maximum response
  real<lower=0> rho; // Lengthscale of treatment response

}
transformed parameters {

  vector[n_batch] mu;
  real theta;
  real beta;

  mu = mu_loc + mu_scale * mu_raw;
  theta = logistic_quantile(theta_raw, theta_l, theta_u, theta_a, theta_b);
  beta = theta + delta;

}
model {

  // Weakly informative priors for distribution of log-odds
  mu_loc ~ normal(0, 5);
  mu_scale ~ normal(0, 1);
  mu_raw ~ normal(0, 1);
  sigma ~ inv_gamma(2, 1);
  a ~ lognormal(0, 1);
  b ~ lognormal(0, 1);

  theta_raw ~ logistic(0, 1);
  delta ~ lognormal(-0.7, 0.5); // Median of 0.5, scale=0.5
  gamma ~ inv_gamma(2, 2);
  rho ~ inv_gamma(3, 10);

  // Calculate mixture probability over existence of treatment effect
  {
    matrix[n_conc, n_conc] Sigma;
    matrix[n_conc, n_conc] L_Sigma;
    for (i in 1:n_conc) {

      // Diagonal terms
      Sigma[i, i] = square(sigma) / n_treatment_batch
      + square(gamma / (1 + exp(- l19 * (conc[i] - beta) / delta)));

      // Off-diagonal terms
      for (j in 1:i-1) {
        Sigma[i, j] = square(gamma) / (1 + exp(- l19 * (conc[i] - beta) / delta))
        / (1 + exp(- l19 * (conc[j] - beta) / delta))
        * exp(-0.5 * square((conc[i] - conc[j]) / rho));

        Sigma[j, i] = Sigma[i, j];
      }
    }
    L_Sigma = cholesky_decompose(Sigma);
    treatment_response ~ multi_normal_cholesky(rep_vector(0, n_conc), L_Sigma);

  }

  // Likelihood of log-odds
  {
    vector[n_sample] z = rep_vector(0, n_sample);
    for (i in 1:n_sample) {
      if (conc_index[i] > 0) {
        z[i] = treatment_response[conc_index[i]];
      }
    }
    target += beta_logistic_lpdf(
      log_odds | mu[batch_index] + z - 10,
      sigma, a, b);
  }


  // Likelihood of counts (Poisson approximation for low counts)
  count[low_count_index] ~ poisson_log(log_odds[low_count_index] + log_total_count);

  // Likelihood of counts (normal approximation for high counts)
  {
    vector[n_high_count] p = inv_logit(log_odds[high_count_index]);
    to_vector(count[high_count_index]) ~ normal(n .* p, sqrt(n .* p .* (1 - p)));
  }

}

