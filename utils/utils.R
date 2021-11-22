library(tidyverse)
library(magrittr)
library(modelr)

prob_obs_lm <- function(y, x, beta_0, beta_1, sigma, log = F){
  mu <- beta_0 + beta_1 * x
  dnorm(y, mean = mu, sd = sigma, log = log)
}

log_prob_obs_lm <- function(y, x, beta_0, beta_1, sigma){
  n <- length(y)
  mu <- beta_0 + beta_1 * x
  -n/2 * log(2*pi) - n/2 * log(sigma^2) - sum((y - mu)^2)/(2*sigma^2)
}

make_data_set <- function(N = 100, K = 10, p = 1, test_proportion = 0.25, lambda = 1){
  
  make_design_matrix <- function(N = 100, K = 25){
    x <- rnorm(N)
    
    X <- vector("list", K)
    for (i in seq(K)){
      x <- scale(lambda * x + rnorm(N))[,1]
      X[[i]] <- x
    }
    bind_cols(X, .name_repair = ~paste0('x_', seq(K))) 
  }
  
  X <- make_design_matrix(N = N, K = K)
  
  w <- c(rep(1/p, p),
         rep(0, K - p))
  
  y <- as.vector(as.matrix(X) %*% w) + rnorm(N)
  
  X <- mutate(X, y = y) %>%  relocate(y)
  
  # split X into train and test sets
  split_vector <- rep('train', N)
  split_vector[sample(seq(N))[1:ceiling(test_proportion * N)]] <- 'test'
  split(X, split_vector)
}

all_subsets_lm <- function(X, outcome, predictors){
  get_variable_inclusion_matrix <- function(X, outcome, predictors){
    
    select_loc <- function(data, ...) {
      tidyselect::eval_select(rlang::expr(c(...)), data)
    }
    
    outcome_ <- enquo(outcome)
    predictors_ <- enquo(predictors)
    
    outcome_str <- select_loc(X, !!outcome_) %>% names()
    
    outcome_str <- select_loc(X, !!predictors_) %>% names()
    
    variable_inclusion_matrix <- expand.grid(replicate(length(outcome_str), 
                                                       c(F, T), 
                                                       simplify = FALSE)) %>% 
      set_names(outcome_str) 
    
    row_names <- paste0('M', seq(0, nrow(variable_inclusion_matrix)-1))
    
    set_rownames(variable_inclusion_matrix, 
                 row_names)
    
  }
  
  make_formula <- function(vec){
    rhs <- unlist(vec) %>% 
      which() %>%
      names() %>% 
      c('1', .) %>% 
      paste(collapse = ' + ')
    
    paste(c(rlang::as_name(enquo(outcome)), rhs), collapse = ' ~ ') 
  }
  
  lm_model <- function(f){
    do.call("lm", list(as.formula(f), data = quote(X)))
  }
  
  outcome_ <- enquo(outcome)
  predictors_ <- enquo(predictors)
  
  variable_inclusion_matrix <- get_variable_inclusion_matrix(X, !!outcome_, !!predictors_)
  
  results <- variable_inclusion_matrix %>% 
    split(rownames(variable_inclusion_matrix)) %>% 
    map(make_formula) %>% 
    map(lm_model)
  
  list(variable_inclusion_matrix = variable_inclusion_matrix %>%
         as_tibble(rownames = 'model'),
       results = results)
}

akaike_weights <- function(aic){
  d_aic = aic - min(aic)
  f <- exp(-d_aic/2)
  f/sum(f)
}

relative_importance_weights <- function(.data){
  .data %>% summarise(across(starts_with('x_'), ~(sum(. * w))))
}

AICc <- function(model){
  
  LL <- logLik(model)
  k <- LL %>% attr('df')
  LL <- as.numeric(LL)
  N <- nrow(model$model)
  
  (-2 * LL + 2 * k) + (2 * k * (k + 1))/(N - k - 1)
  
}

all_subsets_lm_eval <- function(multimodel){
  
  rsq <- function(model) summary(model)$r.squared
  adjrsq <- function(model) summary(model)$adj.r.squared
  sig <- function(model) mean(residuals(model)^2)
  ll <- function(model) logLik(model)
  
  model_eval <- map_dfr(multimodel$results, 
                        ~c(rsq = rsq(.), 
                           adjrsq = adjrsq(.),
                           aic = AIC(.), 
                           aic_c = AICc(.)),
                        .id = 'model') %>% 
    mutate(daic = aic_c - min(aic_c),
           w = akaike_weights(aic_c))
  
  left_join(multimodel$variable_inclusion_matrix,
            model_eval,
            by = 'model') %>% 
    arrange(desc(w))
}


bootstrap_aic <- function(multimodel, top_n = 10, n = 1000, aic = AICc, probs = c(0.025, 0.5, 0.975)){
  
  bootstrap_lm <- function(lm_model, n = 10000){
    boot <- bootstrap(lm_model$model, n)
    map(boot$strap, ~lm(formula(lm_model), data = .))
  }
  
  all_subsets_lm_eval(multimodel) %>% 
    slice_max(w, n = top_n) %>% 
    pull(model) %>% 
    set_names(.,.) %>% 
    map_dfr(~bootstrap_lm(multimodel$results[[.]], n = n) %>%
              map_dbl(aic) %>%
              quantile(probs = probs),
            .id = 'model')
}

elpd_testing_set <- function(multimodel, test_df){
  
  elpd_lm <- function(lm_model, test_data){
    s <- sigma(lm_model)
    add_predictions(test_data, lm_model, var = 'y_hat') %>% 
      mutate(lpd = dnorm(test_data$y, mean = y_hat, sd = s, log = T)) %>% 
      summarise(elpd = sum(lpd)) %>% 
      unlist() %>% 
      multiply_by(-2)
  }
  
  map_dbl(multimodel$results, ~elpd_lm(., test_data = test_df)) %>% 
    enframe(name = 'model', value = 'elpd') %>% 
    arrange(elpd)
  
}



looic <- function(multimodel, top_n = 10){
  
  loo_lm <- function(lm_model){
    
    data_df <- lm_model$model
    
    n <- nrow(data_df)
    
    loocv <- function(i){
      
      train_df <- slice(data_df, -i)
      
      test_df <- slice(data_df, i)
      
      M <- lm(formula(lm_model), data = train_df) 
      predictions <- add_predictions(data = test_df, model = M, var = 'y_hat')
      
      mutate(predictions, 
             lpd = dnorm(y, mean = y_hat, sd = sigma(M), log = T)
      ) %>% pull(lpd) %>% unlist()
      
    }
    
    map_dbl(seq(n), loocv) %>% sum() %>% multiply_by(-2)
    
  }
  
  
  all_subsets_lm_eval(multimodel) %>% 
    slice_max(w, n = top_n) %>% 
    pull(model) %>% 
    set_names(.,.) %>% 
    map_dbl(~loo_lm(multimodel$results[[.]])) %>% 
    enframe(name = 'model', value = 'elpd')
}


rsq_sample <- function(N = 101, predictors, rsq = 'r.squared'){
  data_set <- make_data_set(N, test_proportion = 0.0)
  lm(y ~ ., data = select(data_set$train, y, !!enquo(predictors))) %>% 
    summary() %>% 
    extract2(rsq)
}

logprediction_housing_m1 <- function(i){
  m1_not_i <- lm(log(price) ~ 1, 
                 data = slice(housing_df, -i)
  )
  mu <- coef(m1_not_i)
  stdev <- sigma(m1_not_i)
  y_i <- slice(housing_df, i) %>% pull(price)
  dnorm(log(y_i), mean = mu, sd = stdev, log = T)
  
}

logprediction_housing_m0 <- function(i){
  m0_not_i <- lm(price ~ 1, 
                 data = slice(housing_df, -i)
  )
  mu <- coef(m0_not_i)
  stdev <- sigma(m0_not_i)
  y_i <- slice(housing_df, i) %>% pull(price)
  dnorm(y_i, mean = mu, sd = stdev, log = T)
  
}


lm_loo <- function(m){
  
  data_df <- m$model
  n <- nrow(data_df)
  formula <- formula(m)
  
  lm_drop_i <- function(i){
    m_not_i <- lm(formula, data = slice(data_df, -i))
  }
  
  map(seq(n), lm_drop_i)
  
}

lm_loo_cv <- function(m){
  
  data_df <- m$model
  n <- nrow(data_df)
  lm_formula <- formula(m)
  outcome_var <- all.vars(lm_formula)[1]
  
  lm_drop_i <- function(i){
    m_not_i <- lm(lm_formula, data = slice(data_df, -i))
    slice(data_df, i) %>% 
      add_predictions(m_not_i) %>% 
      transmute(lpd = dnorm(x = .[[outcome_var]], mean = pred, sd = sigma(m_not_i), log = TRUE)) %>% 
      unlist()
  }
  
  map_dbl(seq(n), lm_drop_i) %>% sum() # elpd 
  
}

get_rsq <- function(m) summary(m)$r.sq
get_rss <- function(m) sum(residuals(m)^2)
get_adjrsq <- function(m) summary(m)$adj.r.sq
