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

library(splines)

spline_reg_aic <- map_dbl(seq(3, 30),
        ~AICc(lm(vocab ~ ns(age, df = .), data = vocab_df))
)


min(spline_reg_aic)
which.min(spline_reg_aic)

# plot optimal model
plot(predict(lm(vocab ~ ns(age, df = 6), data = vocab_df)))
plot(predict(lm(vocab ~ ns(age, df = 5), data = vocab_df)))
plot(predict(lm(vocab ~ ns(age, df = 4), data = vocab_df)))

spline_reg_aic - min(spline_reg_aic)

plot(
  predict(lm(vocab ~ ns(age, df = 30), data = vocab_df)),
  type = 'l'
)

map_dbl(seq(3, 30),
        ~get_rsq(lm(vocab ~ ns(age, df = .), data = vocab_df))
)



# Variable selection ------------------------------------------------------

data_df <- make_data_set(N = 1000,
                         K = 10,
                         p = 3,
                         test_proportion = 0.1)

M7 <- lm(y ~ ., data = data_df$train)

M7_step_bw <- step(M7, direction = 'backward')

M8 <- lm(y ~ 1, data = data_df$train)

M8_step_fwd <- step(M8, direction = 'forward', scope = formula(M7))

M9_step_both <- step(M7, direction = 'both')

# all subsets regression# -------------------------------------------------

multimodel <- all_subsets_lm(data_df$train, y, x_1:x_10)

# all 2^10 models are in the `results` element of multimodel
length(multimodel$results)

multimodel_result <- all_subsets_lm_eval(multimodel)

multimodel_result %>% 
  select(-starts_with('x_')) %>% 
  mutate(cumsum_w = cumsum(w)) %>% 
  filter(cumsum_w < 0.9)

multimodel_result

relative_importance_weights(multimodel_result)

bootstrap_aic(multimodel)



# Model averaging ---------------------------------------------------------

library(MuMIn)
M5 <- lm(Fertility ~ ., data = swiss, na.action = 'na.fail')
mmM5 <- dredge(M5)

mmM5_avg <- model.avg(mmM5, subset = delta < 4)
mmM5_avg$coefArray[,,'Agriculture']

confint(mmM5_avg)

# Lasso etc ---------------------------------------------------------------

library(glmnet)

y <- data_df$train$y
X <- as.matrix(dplyr::select(data_df$train, starts_with('x_')))

M10_lasso <- glmnet(X, y, alpha = 1)

plot(M10_lasso, xvar='lambda', label = TRUE)

M10_lasso_lambda_cv <- cv.glmnet(X, y, alpha = 1)
plot(M10_lasso_lambda_cv)

M10_lasso_lambda_cv$lambda.1se %>% log()

M11_lasso <- glmnet(X, y, alpha = 1, lambda = 0.01)
coef(M11_lasso)

M11_lasso_b <- glmnet(X, y, alpha = 1, lambda = 0.2)
coef(M11_lasso_b)

# compare to the full model using mle
lm(y ~ ., data = data_df$train) %>% coef()


# Ridge regression

M12_ridge <- glmnet(X, y, alpha = 0)

plot(M12_ridge, xvar='lambda', label = TRUE)

M12_ridge_lambda_cv <- cv.glmnet(X, y, alpha = 0)
plot(M12_ridge_lambda_cv)

M12_ridge_lambda_cv$lambda.min
M12_ridge_lambda_cv$lambda.1se

M13_ridge <- glmnet(X, y, alpha = 0, lambda = 0.1)
coef(M13_ridge)

M14_ridge <- glmnet(X, y, alpha = 0, lambda = 0.8)
coef(M14_ridge)


# Bayesian methods --------------------------------------------------------

class_df <- read_csv("https://raw.githubusercontent.com/mark-andrews/msms02/main/data/classroom.csv")

library(lme4)

M15 <- lmer(mathscore ~ ses + (ses|schoolid) + (ses|classid), data=class_df)
M15a <- lmer(mathscore ~ ses + (ses|schoolid) + (ses||classid), data=class_df)
M15b <- lmer(mathscore ~ ses + (ses|schoolid) + (1|classid), data=class_df)

library(brms)

M16 <- brm(mathscore ~ ses + (ses|schoolid) + (ses|classid), 
           cores = 4,
           data=class_df)

M17 <- brm(mathscore ~ ses, 
           cores = 4,
           data=class_df)

loo(M16)
loo(M17)

waic(M16)
waic(M17)
