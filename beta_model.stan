data {
  int<lower=0> N; // number of observations
  int<lower=0> I; // number of predictors
  matrix[N,I] X; // matrix of predictors
  vector<lower=0,upper=100>[N] y; // final exam scores
}

transformed data {
  // add an intercept column to the matrix of predictors
  matrix[N,I+1] X_intercept = append_col(rep_vector(1, N), X);
  vector[N] y_offset = y / 100; // final exam scores, scaled between 0 and 1 and offset
  for(n in 1:N) {
    if(y_offset[n] == 1) {
      y_offset[n] = 0.9999;
    } else if(y_offset[n] == 0) {
      y_offset[n] = 0.0001;
    }
  }
}

parameters {
  vector[I+1] beta; // intercept and predictor coefficients
  real<lower=0> phi; // scale parameter
}

model {
  // priors
  beta ~ std_normal();
  phi ~ exponential(1);
  // model
  y_offset ~ beta(inv_logit(X_intercept * beta) * phi,
                  (1 - inv_logit(X_intercept * beta)) * phi);
}

generated quantities {
  real y_pred[N] = beta_rng(inv_logit(X_intercept * beta) * phi,
                            (1 - inv_logit(X_intercept * beta)) * phi);
  for(n in 1:N) {
    y_pred[n] = y_pred[n] * 100;
  }
}
