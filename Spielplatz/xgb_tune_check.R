# Loading data - population and sample data
data("eusilcA_pop")
data("eusilcA_smp")

xgb_tune_model <- tidy_xgb_tune(fixed = eqIncome ~ eqsize + cash + self_empl + unempl_ben + age_ben +
                                  surv_ben + sick_ben + dis_ben + rent + fam_allow + house_allow +
                                  cap_inv + tax_adj,
                                smp_data = eusilcA_smp,
                                smp_weights = "weight",
                                domains = "district",
                                cluster = "district",
                                tune_size = 30,
                                parallel_over = "resamples",
                                folds = 3)


library(parsnip)
library(tune)

# 1️⃣ Fake tune parameter list
tune_params <- list(
  trees = c(100L, 200L),  # should be tuned
  min_n = 5L              # fixed
)

# 2️⃣ Define helper function
define_param <- function(param, tune_params) {
  val <- tune_params[[param]]
  if (!is.null(val) && length(val) > 1) return(tune())
  if (!is.null(val) && length(val) == 1) return(val[1])
  return(NULL)
}

# 3️⃣ Case A: pass function call directly → ends up unevaluated
xgb_wrong <- boost_tree(
  trees = define_param("trees", tune_params),
  min_n = define_param("min_n", tune_params)
)

xgb_wrong$args
# $trees
# <quosure>
# expr: ^define_param("trees", tune_params)   <--- UNEVALUATED
# $min_n
# <quosure>
# expr: ^define_param("min_n", tune_params)   <--- UNEVALUATED

# 4️⃣ Case B: evaluate first → pass values (tune() or numeric)
trees_val <- define_param("trees", tune_params)
min_n_val  <- define_param("min_n", tune_params)

xgb_fixed <- boost_tree(
  trees = trees_val,
  min_n = min_n_val
)

xgb_fixed$args
# $trees
# <quosure>
# expr: ^tune()          <--- correctly detected as tunable
# $min_n
# <quosure>
# expr: ^5               <--- correctly detected as fixed
