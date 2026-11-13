\name{xgb_tune_repeated}
\alias{xgb_tune_repeated}
\title{Tune on repeated cross-validation, the shared estimation convention}
\usage{
xgb_tune_repeated(fixed, smp_data, domains, fold_seeds = c(123, 124, 125),
  draw_seed = 20261010, n_iter = 150, bounds = NULL, fit_seed = NULL, ...)
}
\arguments{
\item{fixed, smp_data, domains, ...}{passed to \code{\link{xgb_tune}}. \code{grid}, \code{search}, \code{seed}, \code{n_iter} and \code{bounds} are set here.}

\item{fold_seeds}{integer vector of fold seeds. Default \code{c(123, 124, 125)}.}

\item{draw_seed}{seed of the one draw of configurations. Default 20261010.}

\item{n_iter}{number of configurations. Default 150.}

\item{bounds}{named list of ranges as in \code{xgb_tune}; \code{NULL} takes the convention's.}

\item{fit_seed}{seed passed to xgboost for every fit; \code{NULL} uses \code{draw_seed}.}
}
\value{
a list of \code{xgb_tune} results named \code{seed_<fold seed>}, with attributes \code{draw_seed}, \code{fold_seeds}, \code{fit_seed} and \code{bounds}.
}
\description{
Draws \code{n_iter} configurations once at \code{draw_seed} from \code{bounds} (default: the convention's bounds), then scores the same configurations with \code{\link{xgb_tune}} on each fold seed in \code{fold_seeds} (ten domain-grouped folds). \code{\link{xgb_top_configs}} averages losses and paired differences over the seeds; the list is passed to it as it is.
}
\seealso{
\code{\link{xgb_tune}}, \code{\link{xgb_top_configs}}
}
