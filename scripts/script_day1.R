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
