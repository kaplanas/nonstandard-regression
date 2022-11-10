data {
  int<lower=0> N; // number of observations
  int<lower=0> I; // number of predictors
  matrix[N,I] X; // matrix of predictors
  int J; // number of possible grades
  int<lower=1,upper=J> y[N]; // final exam grades
}

parameters {
  vector[I] beta; // predictor coefficients
  ordered[J-1] c; // cutpoints
}

model {
  // priors
  beta ~ std_normal();
  // model
  y ~ ordered_logistic(X * beta, c);
}

generated quantities {
  vector[N] y_latent;
  int y_pred[N];
  y_latent = X * beta;
  for(n in 1:N) {
    y_pred[n] = ordered_logistic_rng(y_latent[n], c);
  }
}
