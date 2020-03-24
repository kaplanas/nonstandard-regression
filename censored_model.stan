functions {
  // Get the indices of left-, right-, or non-censored observations.
  // type 1: left
  // type 2: right
  // type 3: non
  int[] inds_censored(vector y, int type) {
    int inds[rows(y)];
    int counter = 0;
    for(n in 1:rows(y)) {
      if((type == 1 && y[n] == 0) ||
         (type == 2 && y[n] == 100) ||
         (type == 3 && y[n] > 0 && y[n] < 100)) {
        counter += 1;
        inds[counter] = n;
      }
    }
    return(inds[1:counter]);
  }
}

data {
  int<lower=0> N; // number of observations
  int<lower=0> I; // number of predictors
  matrix[N,I] X; // matrix of predictors
  vector<lower=0,upper=100>[N] y; // final exam scores
}

transformed data {
  // indices of left-, right-, and non-censored observations
  int inds_left_censored[size(inds_censored(y, 1))] = inds_censored(y, 1);
  int inds_right_censored[size(inds_censored(y, 2))] = inds_censored(y, 2);
  int inds_non_censored[size(inds_censored(y, 3))] = inds_censored(y, 3);
  // add an intercept column to the matrix of predictors
  matrix[N,I+1] X_intercept = append_col(rep_vector(1, N), X);
  // scale outputs between 0 and 1
  // (underflow problems when we try to predict on the original scale)
  vector[N] y_scaled = y / 100;
}

parameters {
  vector[I+1] beta_scaled; // intercept and predictor coefficients
  real<lower=0> sigma_scaled; // scale parameter
}

transformed parameters {
  vector[I+1] beta = beta_scaled * 100;
  real<lower=0> sigma = sigma_scaled * 100;
}

model {
  // priors
  beta_scaled[1] ~ normal(0.5, 0.15);
  segment(beta_scaled, 2, I) ~ std_normal();
  sigma_scaled ~ exponential(5);
  // model
  y_scaled[inds_non_censored] ~ normal(X_intercept[inds_non_censored,] * beta_scaled, sigma_scaled);
  target += normal_lcdf(y_scaled[inds_left_censored] | X_intercept[inds_left_censored,] * beta_scaled,
                        sigma_scaled);
  target += normal_lccdf(y_scaled[inds_right_censored] | X_intercept[inds_right_censored,] * beta_scaled,
                         sigma_scaled);
}

generated quantities {
  real y_pred[N] = normal_rng(X_intercept * beta, sigma);
  for(n in 1:N) {
    if(y_pred[n] > 100) {
      y_pred[n] = 100;
    } else if(y_pred[n] < 0) {
      y_pred[n] = 0;
    }
  }
}
