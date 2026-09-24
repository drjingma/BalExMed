#' BalExMed: Balance-based Mediation Analysis for Compositional Exposures
#'
#' Bayesian mediation analysis for a compositional exposure, such as the
#' relative abundances of microbial taxa, a continuous mediator and a
#' continuous outcome. The composition acts through a latent balance
#' \eqn{B(z, x)}, the normalized log-ratio between the geometric means of two
#' unknown groups of parts, in the structural equation model
#' \deqn{m_i = \alpha B(z, x_i) + \psi_m^\top c_i + e_i,}
#' \deqn{y_i = \gamma B(z, x_i) + \beta m_i + \psi_y^\top c_i + \varepsilon_i.}
#' For a unit increase in the balance the direct effect is \eqn{\gamma}, the
#' indirect effect is \eqn{\alpha\beta} and the total effect is
#' \eqn{\gamma + \alpha\beta}. The parts that define the balance are selected
#' jointly with the effects.
#'
#' Fit the model with [balexmed()] and inspect it with
#' [summary.balexmed()], [taxa_table()], [effect_draws()] and
#' [plot.balexmed()].
#'
#' @keywords internal
"_PACKAGE"
