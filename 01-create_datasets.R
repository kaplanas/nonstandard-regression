# Load packages.
library(dplyr)
library(ggplot2)
library(EnvStats)
library(betareg)
library(tidyr)
library(purrr)
library(broom)
library(censReg)
library(tibble)
library(gamlss)
library(logitnorm)

#########################################
# Parameters needed for all simulations #
#########################################

# Number of simulations.
n.sims = 100

# Distance to adjust scores away from exactly 0 or 1.
score.offset = 0.0001

# Table of location parameters, scale parameters, and coefficients.  Add target
# mean and standard deviation for each collection of parameters.
real.params.df = expand.grid(
  distribution = c("censored", "logit-normal", "beta"),
  location = c("mid", "high"),
  scale = c("narrow", "wide"),
  param = c("location", "scale", "continuous"),
  continuous = seq(0.01, 0.05, 0.01),
  true.value = NA_real_,
  stringsAsFactors = F
) %>%
  mutate(target.mean = case_when(location == "mid" ~ 0.5,
                                 location == "high" ~ 0.8),
         target.sd = case_when(scale == "narrow" ~ 0.05,
                               scale == "wide" ~ 0.2))

# We want the coefficient of the continuous predictor to be associated with
# roughly the same size of effect in the output scores, regardless of the
# distribution that generated the data.  The actual value of this coefficient
# will vary depending on the distribution _and_ on the mean of the output
# scores (for logit-normal and beta distributions).  Unfortunately, there isn't
# a closed-form solution for computing the coefficient based on the mean for
# these latter two (that I know of), so we'll have to do optimization.  Also
# unfortunately, we have to do this in a couple of separate places  below.  So
# let's write a function.
find.continuous.coef = function(distribution, location, target.coef) {
  coef.vals = c()
  for(i in 1:length(distribution)) {
    if(distribution[i] == "censored") {
      coef.vals = c(coef.vals, target.coef[i])
    } else if(distribution[i] %in% c("logit-normal", "beta")) {
      coef.vals = c(
        coef.vals,
        optim(
          target.coef[i],
          function(p, l, tc) {
            abs((invlogit(l + (0.1 * p)) - invlogit(l - (0.1 * p))) -
                  (0.2 * tc))
          },
          l = location[i], tc = target.coef[i], method = "Nelder-Mead"
        )$par
      )
    }
  }
  return(coef.vals)
}

# Parameter values.  Adjusted so that scores will have roughly the same mean
# and standard deviation regardless of the distribution that generated them.
# Also, normalize the coefficient for the continuous predictor so that a change
# of 0.2 around the middle of the distribution results in roughly the same
# change in score regardless of the distribution that generated the scores.
# Finding the right parameters is somewhat brittle, especially for the beta
# distribution, so let's give it multiple opportunities.
real.params.df = real.params.df %>%
  # The number of distinct sets of parameters we need to optimize.
  dplyr::select(distribution, target.mean, target.sd, continuous) %>%
  distinct() %>%
  # For each distinct distribution/location/scale/coefficient size, find the
  # location and scale parameters.  (The coefficient will be computed directly
  # from the mean.)
  mutate(params = pmap(
    list(d = distribution, tm = target.mean, ts = target.sd, tc = continuous),
    function(d, tm, ts, tc) {
      continue = T
      while(continue) {
        optim.result = optim(
          c(tm, log(ts)),
          function(p, d, tm, ts, tc) {
            # The coefficient of the continuous parameters is a function of the
            # target mean.
            c = find.continuous.coef(d, p[1], tc)
            # Set up some dummy continuous predictors.
            df = data.frame(x = rnorm(10000, 0, c))
            # Compute the outputs from the inputs and parameters.
            if(d == "censored") {
              df = df %>%
                mutate(y = rnorm(n(), p[1] + x, exp(p[2]))) %>%
                mutate(y = case_when(y >= 1 ~ 1,
                                     y <= 0 ~ 0,
                                     T ~ y))
            } else if(d == "logit-normal") {
              df = df %>%
                mutate(y = rlogitnorm(n(), p[1] + x, exp(p[2])))
            } else if(d == "beta") {
              df = df %>%
                mutate(y = rbeta(n(),
                                 invlogit(p[1] + x) * (exp(p[2]) * 10),
                                 (1 - invlogit(p[1] + x)) * (exp(p[2]) * 10)))
            }
            # Return the difference between the actual and target mean and
            # standard deviation.
            return(abs(mean(df$y) - tm) +
                     abs(sd(df$y) - ts))
          },
          d = d, tm = tm, ts = ts, tc = tc, method = "Nelder-Mead")
        # if(optim.result$value < 0.1) {
          continue = F
        # }
      }
      optim.result$par %>%
        as.list() %>% setNames(c("location", "scale")) %>%
        data.frame() %>% list()
    }
  )) %>%
  # Extract the parameters we computed.
  unnest(params) %>% unnest(params) %>%
  # Compute the coefficient of the continuous predictor as a function of the
  # location parameter.  Also fix the scale parameter (working with logs makes
  # optimization easier).
  mutate(scale = exp(scale) * ifelse(distribution == "beta", 10, 1),
         cont = find.continuous.coef(distribution, location, continuous)) %>%
  pivot_longer(cols = c(location, scale, cont),
               names_to = "param",
               values_to = "tv") %>%
  mutate(param = ifelse(param == "cont", "continuous", param)) %>%
  # Join back to the main table.
  right_join(real.params.df, by = c("distribution", "param", "continuous",
                               "target.mean", "target.sd")) %>%
  mutate(true.value = tv) %>%
  dplyr::select(distribution, location, scale, continuous, param, true.value,
                target.mean, target.sd)

# Parameters we're going to vary in each simulation.
sims.df = expand.grid(
  n.obs = c(100, 500),
  location = unique(real.params.df$location),
  scale = unique(real.params.df$scale),
  distribution = unique(real.params.df$distribution),
  continuous = unique(real.params.df$continuous),
  sub.sim.id = 1:n.sims,
  stringsAsFactors = F
) %>%
  rownames_to_column("sim.id")

# For each set of parameters, generate a dataset and fit models to it.
for(i in 1:nrow(sims.df)) {
  # Progress tracker.
  if(i %% 100 == 0) { print(paste(i, "/", nrow(sims.df))) }
  # Get parameters for this dataset.
  n.obs = sims.df$n.obs[i]
  temp.params.df = real.params.df %>%
    filter(distribution == sims.df$distribution[i] &
             location == sims.df$location[i] &
             scale == sims.df$scale[i] &
             continuous == sims.df$continuous[i]) %>%
    dplyr::select(param, true.value) %>%
    pivot_wider(names_from = "param", values_from = "true.value")
  # Generate the dataset.
  df = tibble(cs.gpa = rnorm(n.obs, 0, 1)) %>%
    mutate(mu = temp.params.df$location +
             (cs.gpa * temp.params.df$continuous),
           score = case_when(sims.df$distribution[i] == "censored" ~
                               rnorm(n.obs, mu, temp.params.df$scale),
                             sims.df$distribution[i] == "logit-normal" ~
                               rlogitnorm(n.obs, mu, temp.params.df$scale),
                             sims.df$distribution[i] == "beta" ~
                               rbeta(n.obs,
                                     invlogit(mu) * temp.params.df$scale,
                                     (1 - invlogit(mu)) * temp.params.df$scale)))
  # Censor score between 0 and 1.  Add an offset score that's never exactly 0
  # or 1 (necessary for some regression packages).
  df = df %>%
    mutate(score = case_when(score <= 0 ~ 0,
                             score >= 1 ~ 1,
                             T ~ score),
           score.o = case_when(score <= score.offset ~ score.offset,
                               score >= 1 - score.offset ~ 1 - score.offset,
                               T ~ score))
  # Add the dataset to the table of simulations.
  sims.df$data[i] = tibble(df)
  # Fit models to the dataset.  Add these models to the table of simulations.
  temp.fit = lm(score ~ cs.gpa, data = df)
  sims.df$normal.fit[i] = tibble(tidy(temp.fit))
  sims.df$normal.preds[i] = tibble(data.frame(temp.fit[c("fitted.values",
                                                         "residuals")]))
  if(min(df$score) == 0 | max(df$score) == 1) {
    temp.fit = censReg(score ~ cs.gpa,
                       left = 0, right = 1, data = df)
    temp = data.frame(coef(summary(temp.fit))) %>%
      rownames_to_column("term")
    colnames(temp) = c("term", "estimate", "std.error", "statistic", "p.value")
    sims.df$censored.fit[i] = tibble(temp)
    temp.predictors = c("(Intercept)", intersect(temp$term, names(df)))
    temp.preds = as.matrix(cbind(rep(1, nrow(df)), df[,temp.predictors[-1]])) %*%
      as.matrix(temp$estimate[temp$term %in% temp.predictors])
    sims.df$censored.preds[i] = data.frame(fitted.values = temp.preds) %>%
      mutate(fitted.values = case_when(fitted.values > 1 ~ 1,
                                       fitted.values < 0 ~ 0,
                                       T ~ fitted.values),
             residuals = df$score - fitted.values) %>%
      tibble()
  } else {
    sims.df$censored.fit[i] = sims.df$normal.fit[i]
    sims.df$censored.preds[i] = sims.df$normal.preds[i]
  }
  temp.fit = gamlss(score.o ~ cs.gpa,
                    data = df,
                    family = LOGITNO(),
                    control = gamlss.control(trace = F))
  sims.df$logit.fit[i] = tibble(tidy(temp.fit) %>%
                                  filter(parameter == "mu") %>%
                                  dplyr::select(-parameter))
  sims.df$logit.preds[i] = tibble(data.frame(fitted.values = fitted(temp.fit),
                                             residuals = df$score.o - fitted(temp.fit)))
  temp.fit = betareg(score.o ~ cs.gpa,
                     data = df,
                     link = "logit", link.phi = "log")
  sims.df$beta.fit[i] = tibble(tidy(temp.fit))
  sims.df$beta.preds[i] = tibble(data.frame(temp.fit[c("fitted.values",
                                                       "residuals")]))
}
rm(i, df, temp.fit, temp, temp.predictors, temp.preds, temp.params.df, n.obs)
