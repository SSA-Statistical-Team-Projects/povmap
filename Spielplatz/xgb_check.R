# data
library(AmesHousing)

# data cleaning
library(janitor)

# data prep
library(dplyr)

# tidymodels
library(rsample)
library(recipes)
library(parsnip)
library(tune)
library(dials)
library(workflows)
library(yardstick)
library(dplyr)
library(rlang)
library(furrr)
library(purrr)

# speed up computation with parrallel processing (optional)
library(doParallel)
all_cores <- parallel::detectCores(logical = FALSE)
registerDoParallel(cores = all_cores)

# set the random seed so we can reproduce any simulated results.
set.seed(1234)

# load the housing data and clean names
ames_data <- make_ames() %>%
  janitor::clean_names()


# split into training and testing datasets. Stratify by Sale price
ames_split <- rsample::initial_split(
  ames_data,
  prop = 0.2,
  strata = sale_price
)

# preprocessing "recipe"
preprocessing_recipe <-
  recipes::recipe(sale_price ~ ., data = training(ames_split)) %>%
  # convert categorical variables to factors
  recipes::step_string2factor(all_nominal()) %>%
  # combine low frequency factor levels
  recipes::step_other(all_nominal(), threshold = 0.01) %>%
  # remove no variance predictors which provide no predictive information
  recipes::step_nzv(all_nominal()) %>%
  prep()


ames_cv_folds <-
  recipes::bake(
    preprocessing_recipe,
    new_data = training(ames_split)
  ) %>%
  rsample::vfold_cv(v = 5)


# XGBoost model specification
xgboost_model <-
  parsnip::boost_tree(
    mode = "regression",
    trees = 1000,
    min_n = tune(),
    tree_depth = tune(),
    learn_rate = tune(),
    loss_reduction = tune()
  ) %>%
  set_engine("xgboost", objective = "reg:squarederror")


# grid specification
xgboost_params <-
  dials::parameters(
    min_n(),
    tree_depth(),
    learn_rate(),
    loss_reduction()
  )


xgboost_grid <-
  dials::grid_max_entropy(
    xgboost_params,
    size = 60
  )

knitr::kable(head(xgboost_grid))

### defining the workflow
xgboost_wf <-
  workflows::workflow() %>%
  workflows::add_model(xgboost_model) %>%
  workflows::add_formula(sale_price ~ .)

# hyperparameter tuning
xgboost_tuned <- tune::tune_grid(
  object = xgboost_wf,
  resamples = ames_cv_folds,
  grid = xgboost_grid,
  metrics = yardstick::metric_set(rmse, rsq, mae),
  control = tune::control_grid(verbose = TRUE)
)

### best hyperparameters at minimizing RMSE
xgboost_tuned %>%
  tune::show_best(metric = "rmse") %>%
  knitr::kable()


### next isolating the very best hyperparameter values
xgboost_best_params <-
  xgboost_tuned %>%
  tune::select_best(metric = "rmse")

knitr::kable(xgboost_best_params)

### finalize the xgboost model using the best tuning parameters
xgboost_model_final <-
  xgboost_model %>%
  finalize_model(xgboost_best_params)


### lets look at the performance of the model on training data
train_processed <- bake(preprocessing_recipe,  new_data = training(ames_split))

train_prediction <-
  xgboost_model_final %>%
  # fit the model on all the training data
  fit(
    formula = sale_price ~ .,
    data    = train_processed
  ) %>%
  # predict the sale prices for the training data
  predict(new_data = train_processed) %>%
  bind_cols(training(ames_split))

xgboost_score_train <-
  train_prediction %>%
  yardstick::metrics(sale_price, .pred) %>%
  mutate(.estimate = format(round(.estimate, 2), big.mark = ","))

knitr::kable(xgboost_score_train)

### and then on the testing data
test_processed  <- bake(preprocessing_recipe, new_data = testing(ames_split))
test_prediction <- xgboost_model_final %>%
  # fit the model on all the training data
  fit(
    formula = sale_price ~ .,
    data    = train_processed
  ) %>%
  # use the training model fit to predict the test data
  predict(new_data = test_processed) %>%
  bind_cols(testing(ames_split))
# measure the accuracy of our model using `yardstick`
xgboost_score <-
  test_prediction %>%
  yardstick::metrics(sale_price, .pred) %>%
  mutate(.estimate = format(round(.estimate, 2), big.mark = ","))
knitr::kable(xgboost_score)


### xgboost for poverty mapping

survey_dt <- eusilcA_smp
census_dt <- eusilcA_pop


traindt_split <-
  survey_dt |>
  rsample::initial_split(prop = 0.9,
                         strata = "eqIncome")


xvars_list <- colnames(survey_dt)[!colnames(survey_dt) %in%
                                    c("weight", "district",
                                    "state", "eqIncome")]

formula_obj <- paste0("eqIncome ~ ",
                      paste(xvars_list,
                            collapse = " + ")) |>
  as.formula()


preprocess_recipe <-
  recipes::recipe(formula_obj, data = training(traindt_split)) %>%
  # convert categorical variables to factors
  recipes::step_string2factor(all_nominal()) %>%
  # combine low frequency factor levels
  recipes::step_other(all_nominal(), threshold = 0.01) %>%
  # remove no variance predictors which provide no predictive information
  recipes::step_nzv(all_nominal()) %>%
  prep()


survey_cvfolds_obj <-
  bake(preprocess_recipe,
       new_data = training(traindt_split)) |>
  rsample::vfold_cv(v = 5)

# -----------------------------
# 1️⃣ XGBoost model specification
# -----------------------------
# xgboost_model <-
#   boost_tree(mode = "regression",
#              trees = tune(),
#              mtry = NULL,
#              min_n = tune(),
#              tree_depth = tune(),
#              learn_rate = tune(),
#              loss_reduction = tune(),
#              sample_size = tune(),
#              stop_iter = NULL) %>%
#   set_engine("xgboost",
#              objective = "reg:squarederror",
#              alpha = tune(),
#              lambda = tune(),
#              colsample_bytree = tune(),
#              colsample_bylevel = tune(),
#              colsample_bynode = tune(),
#              max_delta_step = tune(),
#              counts = FALSE)

# -----------------------------
# 2️⃣ Define tuning parameters
# -----------------------------
xgb_params <- dials::parameters(
  min_n(),
  tree_depth(),
  learn_rate(),
  loss_reduction(),
  sample_prop(),               # maps to sample_size
  alpha = penalty_L1(),
  lambda = penalty_L2(),
  colsample_bytree = new_quant_param(
    type = "double",
    range = c(0.5, 1),
    trans = NULL,
    inclusive = c(T, T)
  ),
  max_delta_step = new_quant_param(
    type = "double",
    range = c(0, 10),
    trans = NULL,
    inclusive = c(T, T)
  )
)

# -----------------------------
# 3️⃣ Create tuning grid
# -----------------------------
xgboost_grid <- grid_space_filling(
  xgb_params,
  size = 60
)

# -----------------------------
# 4️⃣ Define workflow
# -----------------------------

# xgboost_model <-
#   parsnip::boost_tree(
#     mode = "regression",
#     trees = 1000,
#     min_n = tune(),
#     tree_depth = tune(),
#     learn_rate = tune(),
#     loss_reduction = tune()
#   ) %>%
#   set_engine("xgboost", objective = "reg:squarederror")

xgboost_wf <- workflow() %>%
  add_model(xgboost_model) %>%
  add_formula(formula_obj)  # your regression formula

# -----------------------------
# 5️⃣ Hyperparameter tuning
# -----------------------------
xgboost_tuned <- tune_grid(
  object = xgboost_wf,
  resamples = survey_cvfolds_obj,     # your CV object
  grid = xgboost_grid,
  metrics = metric_set(rmse, rsq, mae),
  control = control_grid(verbose = TRUE)
)













### best hyperparameters at minimizing RMSE
xgboost_tuned %>%
  tune::show_best(metric = "rmse") %>%
  knitr::kable()


### next isolating the very best hyperparameter values
xgboost_best_params <-
  xgboost_tuned %>%
  tune::select_best(metric = "rmse")

knitr::kable(xgboost_best_params)

### finalize the xgboost model using the best tuning parameters
xgboost_model_final <-
  xgboost_model %>%
  finalize_model(xgboost_best_params)


### lets look at the performance of the model on training data
train_processed <- bake(preprocess_recipe,  new_data = training(traindt_split))

train_prediction <-
  xgboost_model_final %>%
  # fit the model on all the training data
  fit(
    formula = formula_obj,
    data    = train_processed
  ) %>%
  # predict the sale prices for the training data
  predict(new_data = train_processed) %>%
  bind_cols(training(traindt_split))

xgboost_score_train <-
  train_prediction %>%
  yardstick::metrics(eqIncome, .pred) %>%
  mutate(.estimate = format(round(.estimate, 2), big.mark = ","))

knitr::kable(xgboost_score_train)

### and then on the testing data
test_processed  <- bake(preprocess_recipe, new_data = testing(traindt_split))
test_prediction <- xgboost_model_final %>%
  # fit the model on all the training data
  fit(
    formula = formula_obj,
    data    = train_processed
  ) %>%
  # use the training model fit to predict the test data
  predict(new_data = test_processed) %>%
  bind_cols(testing(traindt_split))
# measure the accuracy of our model using `yardstick`
xgboost_score <-
  test_prediction %>%
  yardstick::metrics(eqIncome, .pred) %>%
  mutate(.estimate = format(round(.estimate, 2), big.mark = ","))
knitr::kable(xgboost_score)



### and then on the census data
test_processed <- bake(preprocess_recipe,
                       new_data = census_dt |> as_tibble() |> dplyr::select(-eqIncome))

test_prediction <- xgboost_model_final %>%
  # fit the model on all the training data
  fit(
    formula = formula_obj,
    data    = train_processed
  ) %>%
  # use the training model fit to predict the test data
  predict(new_data = test_processed) %>%
  bind_cols(census_dt |> as_tibble() |> dplyr::select(-eqIncome))















