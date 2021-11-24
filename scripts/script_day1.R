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
