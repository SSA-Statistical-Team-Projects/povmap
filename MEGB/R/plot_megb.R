

#' Plot Random Effects, Residuals, Train vs. Validation Loss, and Feature Importance
#'
#' @param megb_object object of class MEGB which contains the boosting model
#' @param ... additional arguments
#' @import ggplot2
#' @import stats
#' @importFrom gridExtra grid.arrange
#' @importFrom xgboost xgb.plot.importance
#' @return series of four plots, i.e. a QQ-Plot of the Random Effects, a QQPlot of the Residuals, plot containing the loss for train and validation data, a feature importance plot
#' @export
#' @method plot MEGB
#' @name plot.MEGB


plot.MEGB <- function(megb_object) {

  # QQ Plot for Random Effects
  qqplot_random_effects <- function(megb_object) {
    qqnorm(megb_object$megb_model$dom_name_effects[[1]][, 1], main = "QQPlot for Random Effects")
    qqline(megb_object$megb_model$dom_name_effects[[1]][, 1])
  }
  
  # QQ Plot for Residuals
  qqplot_residuals <- function(megb_object) {
    qqnorm(megb_object$megb_model$residuals, main = "QQPlot for Residuals")
    qqline(megb_object$megb_model$residuals)
  }
  
  # Train vs. Validation Loss Plot
  plot_theme <- theme(
    plot.title = element_text(
      size = 13,
      hjust = 0.5,
      face = "bold"
    ),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )
  
  # Train vs. Validation Loss Plot
  train_validation_loss_plot <- function(megb_object) {
    plot <- megb_object$megb_model$eval_log %>%
      ggplot(aes(x = iter)) +
      geom_line(aes(y = train_rmse_mean, color = "Train"), lwd = 1) +
      geom_line(aes(y = test_rmse_mean, color = "Validation"), lwd = 1) +
      labs(
        title = "Train vs. Validation Loss
 XGBoost (inside EM-Algorithm)",
        x = "Iteration",
        y = "RMSE",
        color = "Dataset"
      ) +
      theme_minimal() +
      plot_theme + # Einheitliches Styling anwenden
      geom_vline(xintercept = megb_object$megb_model$boosting$niter,
                 col = "black") +
      annotate(
        "text",
        x = megb_object$megb_model$boosting$niter,
        y = max(megb_object$megb_model$eval_log$train_rmse_mean),
        label = "Used nrounds",
        vjust = -0.5,
        hjust = 1,
        color = "black",
        size = 3
      )
    
    # print plot
    print(plot)
  }
  
  # Feature Importance Plot
  importance_plot <- function(megb_object, ...) {
    # Plot feature importance
    importance <- megb_object$megb_model$importance_matrix
    xgb.plot.importance(importance,
                        main = "Feature Importance",
                        xlab = "Importance",
                        ylab = "Features",
                        ...)
  }
  
  # List of plot functions
  plot_functions <- list(
    qqplot_random_effects,
    qqplot_residuals,
    train_validation_loss_plot,
    importance_plot
  )
  if (megb_object$gbm_engine != "xgboost") {
    warning(
      "Boosting-related plots are only available for xgboost. ",
      "For catboost or lightgbm, please access boosting information via ",
      "`megb_object$megb_model$boosting`."
    )
    plot_functions <- plot_functions[1:2]
  }
  # plot each function:
  for (plot_fun in plot_functions) {
    # Plot aufrufen
    plot_fun(megb_object)
    
    # wait for Enter
    readline(prompt = "Press ENTER to continue...")
  }
}
