data {
  int<lower=0> N; // number of observations
  int<lower=0> I; // number of predictors
  matrix[N,I] X; // matrix of predictors
  int J; // number of possible grades
  real<lower=1,upper=J> y[N]; // final exam grades
}

transformed data {
  // add an intercept column to the matrix of predictors
  matrix[N,I+1] X_intercept = append_col(rep_vector(1, N), X);
}

parameters {
  vector[I+1] beta; // intercept and predictor coefficients
  real<lower=0> sigma; // scale parameter
}

model {
  // priors
  beta[1] ~ normal(3, 1);
  segment(beta, 2, I) ~ std_normal();
  sigma ~ exponential(0.5);
  // model
  y ~ normal(X_intercept * beta, sigma);
}

generated quantities {
  real y_pred[N] = round(normal_rng(X_intercept * beta, sigma));
}
