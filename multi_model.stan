data {
  int<lower=0> N; // number of observations
  int<lower=0> I; // number of predictors
  matrix[N,I] X; // matrix of predictors
  int<lower=1> J; // number of possible meeting day patterns
  int<lower=1,upper=J> y[N]; // meeting days
}

transformed data {
  // add an intercept column to the matrix of predictors
  matrix[N,I+1] X_intercept = append_col(rep_vector(1, N), X);
}

parameters {
  matrix[I+1,J-1] beta; // predictor coefficients for each possible outcome except one
}

transformed parameters {
  matrix[I+1,J] beta_baseline; // add zeros for the first outcome category
  beta_baseline = append_col(rep_vector(0, I + 1), beta);
}

model {
  // priors
  to_vector(beta) ~ std_normal();
  // model
  for(n in 1:N) {
    y[n] ~ categorical_logit((X_intercept[n,] * beta_baseline)');
  }
}

generated quantities {
  int y_pred[N];
  for(n in 1:N) {
    y_pred[n] = categorical_logit_rng((X_intercept[n,] * beta_baseline)');
  }
}
