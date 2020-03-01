# Load packages.
library(dplyr)
library(ggplot2)
library(forcats)
library(tidyr)
library(logitnorm)

# Set the theme.
theme_set(theme_bw())

# Get distributional properties of the simulated datasets.
score.dist.df = sims.df %>%
  unnest(data) %>%
  group_by(sim.id, n.obs, location, scale, continuous, distribution,
           sub.sim.id) %>%
  summarize(mean.score = mean(score),
            mean.score.o = mean(score.o),
            sd.score = sd(score),
            sd.score.o = sd(score.o)) %>%
  ungroup() %>%
  left_join(params.df %>%
              dplyr::select(distribution, location, scale, continuous,
                            target.mean, target.sd) %>%
              distinct(),
            by = c("distribution", "location", "scale", "continuous")) %>%
  mutate(location = fct_relevel(location, "mid", "high"),
         scale = fct_relevel(scale, "wide", "narrow"),
         distribution = fct_relevel(distribution, "censored"))

# Check the mean and standard deviation of the scores in each simulated
# dataset.
score.dist.df %>%
  group_by(n.obs, location, scale, continuous, distribution, target.mean) %>%
  summarize(upper.95 = quantile(mean.score, 0.975),
            upper.50 = quantile(mean.score, 0.75),
            median = quantile(mean.score, 0.5),
            lower.50 = quantile(mean.score, 0.25),
            lower.95 = quantile(mean.score, 0.025)) %>%
  ggplot(aes(x = as.factor(continuous), y = median)) +
  geom_point(size = 3) +
  geom_linerange(aes(ymin = lower.50, ymax = upper.50), size = 2) +
  geom_linerange(aes(ymin = lower.95, ymax = upper.95)) +
  geom_hline(aes(yintercept = target.mean)) +
  facet_grid(scale + n.obs ~ location + distribution) +
  scale_x_discrete("Size of coefficient for continuous predictor") +
  scale_y_continuous("Mean score", breaks = seq(0, 1, 0.1)) +
  coord_flip()
score.dist.df %>%
  group_by(n.obs, location, scale, continuous, distribution, target.sd) %>%
  summarize(upper.95 = quantile(sd.score, 0.975),
            upper.50 = quantile(sd.score, 0.75),
            median = quantile(sd.score, 0.5),
            lower.50 = quantile(sd.score, 0.25),
            lower.95 = quantile(sd.score, 0.025)) %>%
  ggplot(aes(x = as.factor(continuous), y = median)) +
  geom_point(size = 3) +
  geom_linerange(aes(ymin = lower.50, ymax = upper.50), size = 2) +
  geom_linerange(aes(ymin = lower.95, ymax = upper.95)) +
  geom_hline(aes(yintercept = target.sd)) +
  facet_grid(location + n.obs ~ scale + distribution) +
  scale_x_discrete("Size of coefficients for continuous predictor") +
  scale_y_continuous("Standard deviation of score") +
  coord_flip()
sims.df %>%
  filter(continuous == 0.05 & sub.sim.id == 1) %>%
  mutate(distribution = fct_relevel(distribution, "censored", "logit-normal"),
         location = fct_relevel(location, "mid", "high"),
         scale = fct_relevel(scale, "wide", "narrow")) %>%
  unnest(data) %>%
  ggplot(aes(x = cs.gpa, y = score)) +
  geom_point(alpha = 0.1) +
  stat_smooth() +
  scale_x_continuous(breaks = seq(-5, 5, 1)) +
  facet_grid(location + distribution ~ scale + n.obs)

# Get fitted parameters.
fitted.params.df = bind_rows(
  normal = sims.df %>%
    unnest(normal.fit),
  censored = sims.df %>%
    unnest(censored.fit),
  logit_normal = sims.df %>%
    unnest(logit.fit),
  beta = sims.df %>%
    unnest(beta.fit) %>%
    filter(component == "mean") %>%
    dplyr::select(-component),
  .id = "fit"
) %>%
  dplyr::select(fit, n.obs, location, scale, continuous, distribution, sim.id,
                sub.sim.id, term, estimate, std.error) %>%
  mutate(distribution = gsub("_", "-", distribution),
         fit = gsub("_", "-", fit),
         param = case_when(term == "(Intercept)" ~ "location",
                           term == "cs.gpa" ~ "continuous",
                           T ~ term)) %>%
  dplyr::select(-term) %>%
  mutate(distribution = fct_relevel(distribution, "censored", "logit-normal"),
         fit = fct_relevel(fit, "normal", "censored", "logit-normal"),
         location = fct_relevel(location, "mid", "high"),
         scale = fct_relevel(scale, "wide", "narrow"))

# How often is the continuous coefficient significantly greater than 0?
fitted.params.df %>%
  filter(param == "continuous") %>%
  mutate(found.signif = as.numeric(estimate + (qnorm(0.025) * std.error) > 0)) %>%
  group_by(distribution, fit, n.obs, location, scale, continuous) %>%
  summarize(prop.found.signif = sum(found.signif) / n()) %>%
  ungroup() %>%
  mutate(fit.numeric = as.numeric(fit),
         match = as.character(fit) == as.character(distribution),
         line.width = ifelse(match, 3, 1),
         line.alpha = ifelse(match, 0.5, 1)) %>%
  ggplot(aes(x = as.factor(continuous), y = prop.found.signif, color = fit,
             size = line.width, alpha = line.alpha)) +
  geom_line(aes(group = fit)) +
  scale_x_discrete("Size of coefficient of GPA") +
  scale_y_continuous("Proportion found significant effect",
                     limits = c(0, 1)) +
  scale_color_discrete("Model fit to data") +
  scale_size_identity("Model matches true distribution",
                      guide = "legend",
                      breaks = c(1, 3),
                      labels = c("No", "Yes")) +
  scale_alpha_identity() +
  facet_grid(location + distribution ~ scale + n.obs)

# How often does the normal fit predict values greater than 1?
outside.bounds.df = sims.df %>%
  unnest(normal.preds) %>%
  dplyr::select(sim.id, n.obs, location, scale, distribution, continuous,
                sub.sim.id, fitted.values, residuals) %>%
  group_by(n.obs, location, scale, distribution, continuous) %>%
  summarize(n.over = sum(fitted.values > 1),
            n.under = sum(fitted.values < 0)) %>%
  ungroup() %>%
  mutate(distribution = fct_relevel(distribution, "censored", "logit-normal",
                                    "beta"),
         location = fct_relevel(location, "mid", "high"),
         scale = fct_relevel(scale, "wide", "narrow"))
outside.bounds.df %>%
  ggplot(aes(x = as.factor(continuous), y = n.over)) +
  geom_bar(stat = "identity") +
  scale_x_discrete("Size of coefficient of GPA") +
  scale_y_continuous("Number of predictions > 1") +
  facet_grid(location + distribution ~ scale + n.obs)

# How do the residuals from different fits compare?
preds.df = bind_rows(
  normal = sims.df %>%
    unnest(cols = c(data, normal.preds)),
  censored = sims.df %>%
    unnest(cols = c(data, censored.preds)),
  logit_normal = sims.df %>%
    unnest(cols = c(data, logit.preds)),
  beta = sims.df %>%
    unnest(cols = c(data, beta.preds)),
  .id = "fit"
) %>%
  dplyr::select(fit, sim.id, n.obs, location, scale, distribution, continuous,
                sub.sim.id, score, fitted.values, residuals) %>%
  mutate(distribution = gsub("_", "-", distribution),
         fit = gsub("_", "-", fit)) %>%
  mutate(distribution = fct_relevel(distribution, "censored", "logit-normal"),
         fit = fct_relevel(fit, "normal", "censored", "logit-normal"),
         location = fct_relevel(location, "mid", "high"),
         scale = fct_relevel(scale, "wide", "narrow"))
preds.df %>%
  filter(sub.sim.id == 1 & continuous == 0.05) %>%
  ggplot(aes(x = fitted.values, y = score, color = fit, fill = fit)) +
  geom_point(alpha = 0.05) +
  geom_abline(intercept = 0, slope = 1) +
  stat_smooth(alpha = 0.5) +
  scale_color_discrete("Model fit to data") +
  scale_fill_discrete("Model fit to data") +
  facet_grid(location + distribution ~ scale + n.obs)
preds.df %>%
  filter(continuous == 0.05) %>%
  ggplot(aes(x = residuals, fill = fit)) +
  geom_histogram(position = "identity", alpha = 0.5) +
  scale_fill_discrete("Model fit to data") +
  facet_grid(location + distribution ~ scale + n.obs, scales = "free")
preds.df %>%
  group_by(fit, n.obs, location, scale, distribution, continuous) %>%
  summarize(mean.resid = mean(abs(residuals))) %>%
  ggplot(aes(x = as.factor(continuous), y = mean.resid, fill = fit)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_discrete("Size of coefficient of GPA") +
  scale_y_continuous("Mean absolute residual") +
  scale_fill_discrete("Model fit to data") +
  facet_grid(location + distribution ~ scale + n.obs)
