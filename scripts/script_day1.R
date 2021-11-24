library(tidyverse)

y <- 50
beta_0 <- -20
beta_1 <- 4
varsigma <- 15
x <- 15

# the mean of the normal distribution over `dist` is
mu <- beta_0 + beta_1 * x

# the standard deviation of that normal distribution is varsigma

# the probability (density) of dist = 50 is
1/sqrt(2 * pi * varsigma^2) * exp(-(y - mu)^2/(2*varsigma^2))

# the easy way
dnorm(y, mean = mu, sd = varsigma)

# calculate the probability of y
# given x and assumed values of the parameters of 
# the simple linear model
prob_obs_lm <- function(y, x, beta_0, beta_1, sigma, log = F){
  mu <- beta_0 + beta_1 * x
  dnorm(y, mean = mu, sd = sigma, log = log)
}

# e.g. 
prob_obs_lm(y, x, beta_0, beta_1, varsigma)


# Log probability of *all* the observed data ------------------------------

y <- cars$dist
x <- cars$speed

# 
prob_obs_lm(y, x, beta_0, beta_1, varsigma, log = TRUE) %>% sum()

# algebraic equivalent calculation
log_prob_obs_lm <- function(y, x, beta_0, beta_1, sigma){
  n <- length(y)
  mu <- beta_0 + beta_1 * x
  -n/2 * log(2*pi) - n/2 * log(sigma^2) - sum((y - mu)^2)/(2*sigma^2)
}

log_prob_obs_lm(y, x, beta_0, beta_1, varsigma)



# Maximum likelihood estimates for the parameters of the normal li --------

M1 <- lm(dist ~ speed, data = cars)
coef(M1) # mle of beta_0 and beta_1

beta_0_mle <- coef(M1)['(Intercept)']
beta_1_mle <- coef(M1)['speed']

# mle of sigma 
sigma_mle <- sqrt(mean(residuals(M1)^2))

# not to be confused with sigma(M1)
sigma(M1)

# Probability of the observed values of the outcome 
# variable (i.e. `dist`, aka `y`), given the observed
# values of `speed` (aka `x`), and the given the 
# maximum likelihood values of the parameters
prob_obs_lm(y, x, beta_0_mle, beta_1_mle, sigma_mle, log = TRUE) %>% 
  sum()

logLik(M1)

# work it out ourselves
rss <- sum(residuals(M1)^2)
n <- length(y)

-(n/2) * (log(2*pi) - log(n) + log(rss) + 1)



# Deviance and GLMs -------------------------------------------------------
cars_df <- mutate(cars, z = dist > median(dist))

M2 <- glm(z ~ speed, data = cars_df,
          family = binomial(link = 'logit')
)

logLik(M2)
logLik(M2) * -2
deviance(M2)


theta_mle <- predict(M2, type = 'response')
# pull out the binary variable z
z <- cars_df$z

sum(z * log(theta_mle) + (1 - z) * log(1 - theta_mle))

# Deviance residuals

residuals(M2) # 



# Nested model comparison; normal linear models ---------------------------

M0 <- lm(Fertility ~ Agriculture + Education, data = swiss)
M1 <- lm(Fertility ~ Agriculture + Education + Catholic, data = swiss)

rss_0 <- sum(residuals(M0)^2)
rss_1 <- sum(residuals(M1)^2)

c(rss_0, rss_1)

# proportional increase in error (PIE)
(rss_0 - rss_1)/rss_1

# residual degrees of freedom
M0$df.residual
df.residual(M0)
nrow(swiss) - 2 - 1
nrow(swiss) - 3 - 1
M1$df.residual
df.residual(M1)

df_0 <- df.residual(M0)
df_1 <- df.residual(M1)

f_statistic <-  ((rss_0 - rss_1) / (df_0 - df_1)) / (rss_1/df_1)

c(f_statistic, df_0 - df_1, df_1)
pf(f_statistic, df_0 - df_1, df_1, lower.tail = F)

# or do it the easy way
anova(M0, M1)

# F tests via drop1
drop1(M1, scope = ~ Catholic, test = 'F')

drop1(M1, scope = ~ Education, test = 'F')


# what do values from an F dist (1, 43) look like
quantile(rf(10000, 1, 43), probs = c(0.025, 0.5, 0.975)) %>% round(3)


drop1(M1, scope = ~ Education + Agriculture, test = 'F')

# rss of M1 (intact)
sum(
  residuals(
  lm(Fertility ~ Agriculture + Education + Catholic, data = swiss)
)^2)

# rss of M1 after Education is dropped
sum(
  residuals(
    lm(Fertility ~ Agriculture + Catholic, data = swiss)
  )^2)

# rss of M1 after Agriculture is dropped
sum(
  residuals(
    lm(Fertility ~ Education + Catholic, data = swiss)
  )^2)

# another set of model comparisons 
anova(M0)


rss_null <- sum(
  residuals(
    lm(Fertility ~ 1, data = swiss)
  )^2)

rss_ag_only <- 
  sum(
    residuals(
      lm(Fertility ~ Agriculture, data = swiss)
    )^2)

rss_ag_educ <- 
  sum(
    residuals(
      lm(Fertility ~ Agriculture + Education, data = swiss)
    )^2)

# R squared ---------------------------------------------------------------

M_null <- lm(Fertility ~ 1, data = swiss)
rss_null <- sum(residuals(M_null)^2)
(rss_null - rss_0)/rss_null # R^2 for model M0
summary(M0)$r.squared

(rss_null - rss_1)/rss_null # R^2 for model M1
summary(M1)$r.squared



# Read in some utils ------------------------------------------------------

source("https://raw.githubusercontent.com/mark-andrews/msms02/main/utils/utils.R")



# Generate data and calculate R^2 -----------------------------------------

rsq_sample(N = 101, predictors = x_1:x_5)
# generate artificial data with 101 rows (100 reps)
# 5 predictors, only one non-null
quantile(map_dbl(seq(100), ~rsq_sample(N = 101, predictors = x_1:x_5)))

quantile(map_dbl(seq(100), ~rsq_sample(N = 101, predictors = x_1:x_10)))

quantile(map_dbl(seq(100), ~rsq_sample(N = 50, predictors = x_1:x_10)))

quantile(map_dbl(seq(100), ~rsq_sample(N = 50, predictors = x_1:x_10, rsq = 'adj.r.squared')))

# calculate adj R^2 
summary(M1)$adj.r.squared

rss_over_tss <-  1 - summary(M1)$r.squared
n <- nrow(M1$model)
K <- 3
penalty <- (n - 1) / (n - K - 1)
1 - rss_over_tss * penalty



# Generalized linear models and deviance ----------------------------------

swiss_df <- mutate(swiss, y = Fertility > median(Fertility))

M0 <- glm(y ~ Agriculture + Education, data = swiss_df, 
          family = binomial(link = 'logit'))
M1 <- glm(y ~ Agriculture + Education + Catholic, data = swiss_df, 
          family = binomial(link = 'logit'))

c(logLik(M0), logLik(M1))
c(deviance(M0), deviance(M1))

D_0 <- deviance(M0)
D_1 <- deviance(M1)


# proportional increase in error
(D_0 - D_1)/D_1

# proportional decrease in error
(D_0 - D_1)/D_0 # pseudo R^2 (assuming D_0 is the deviance of the null null)

c(D_0, D_1)

D_0 - D_1 # test statistic for a null hypothesis test

# prob of a value > D_0 - D_1 in a 1df chi ^ 2 dist
pchisq(D_0 - D_1, df = 1, lower.tail = F)

anova(M0, M1, test = 'Chisq')
