library(tidyverse)
housing_df <- read_csv("https://raw.githubusercontent.com/mark-andrews/msms02/main/data/housing.csv")

hist(housing_df$price, 25)
hist(log(housing_df$price), 25)

M1 <- lm(log(price) ~ 1, data = housing_df)

i <- 42
# fit the model using all but row 42
M2 <- lm(log(price) ~ 1,
         data = slice(housing_df, -i)
)

# what is log probability of obs 42 given
# the model fitted with all but obs 42
mu <- coef(M2)
varsigma <- sigma(M2)

obs_i <- slice(housing_df, i) %>% unlist()

dnorm(log(obs_i), 
      mean  = mu, sd = varsigma, log = T)



# Non-nested model comparison using AIC -----------------------------------

# log normal model of house prices
M3 <- lm(log(price) ~ 1, data = housing_df)
# normal model of house prices
M4 <- lm(price ~ 1, data = housing_df)

D3 <- -2 * logLik(M3)
D4 <- -2 * logLik(M4)

K3 <- length(coef(M3)) + 1
K4 <- length(coef(M4)) + 1

# AIC of M3
2 * K3 + D3

# AIC of M4
2 * K4 + D4

# alternatively
AIC(M3)
AIC(M4)



# source our utils.R ------------------------------------------------------

source("https://raw.githubusercontent.com/mark-andrews/msms02/main/utils/utils.R")


M5 <- lm(Fertility ~ ., data = swiss) # full model of the swiss Fertility data
M6 <- lm(Fertility ~ Agriculture + Examination + Education, data = swiss)

all_loo_splits <- lm_loo(M5)
class(all_loo_splits)
all_loo_splits[[1]]
length(all_loo_splits)

lm_loo_cv(M5) # elpd 
loo_ic_M5 <- lm_loo_cv(M5) * -2 # put elpd on the deviance scale

loo_ic_M6 <- lm_loo_cv(M6) * -2

# Compare -2 x elpd from loo cv with AIC
aic_M5 <- AIC(M5)
aic_M6 <- AIC(M6)

aic_c_M5 <- AICc(M5)
aic_c_M6 <- AICc(M6)


# Using AIC etc for nonlinear regression complexity -----------------------

vocab_df <- read_csv("https://raw.githubusercontent.com/mark-andrews/msms02/main/data/GSSvocab.csv")


ggplot(vocab_df,
       aes(x = age, y = vocab)
) + geom_point()
