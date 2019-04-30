#' Causal Effect Decomposition Based on a Fitted \code{rwrmed} Model
#'
#' \code{decomp} is a function that implements causal effect decomposition based on a fitted
#'   \code{rwrmed} model. It returns a two-component decomposition of the total effect into
#'   the randomized interventional analogues of the natural direct effect (rNDE) and the natural
#'   indirect effect (rNIE). It also returns a four-component decomposition of the total effect into
#'   the controlled direct effect (CDE) and the randomized analogues of the reference interaction
#'   effect (rREFINT), the mediated interaction effect (rMEDINT), and the pure indirect effect (rPIE).
#'
#' @param object An object of class \code{rwrmed}.
#' @param a0 The baseline level of treatment.
#' @param a1 The level of treatment to be contrasted with the baseline.
#' @param m The level of the mediator at which the CDE is evaluated.
#'
#' @return A list of two elements.
#'  \item{twocomp}{Two component decomposition of the rATE into rNDE and rNIE.}
#'  \item{fourcomp}{Four component decomposition of the rATE into CDE, rINTREF, rPIE, and rINTMED.}
#' @import stats
#' @export
#' @seealso \code{\link{rwrmed}} for implementing the regression-with-residuals (RWR)
#'   approach to causal mediation.
#'
decomp <- function(object, a0 = 0, a1 = 1, m = 0){

  var_names <- object$var_names

  # mediator and outcome model coefs
  coefs_y <- coef(object$y_model)
  coefs_m <- coef(object$m_model)

  theta0 <- coefs_m["(Intercept)"]
  theta2 <- coefs_m[var_names$treatment]

  beta2 <- coefs_y[var_names$treatment]
  beta4 <- coefs_y[var_names$mediator]
  beta5 <- (coefs_y[paste0(var_names$treatment,":",var_names$mediator)] %||%
              coefs_y[paste0(var_names$mediator,":",var_names$treatment)]) %||% 0

  # effect decomposition
  CDE <- (a1 - a0) * (beta2 + beta5 * m)
  rINTREF <- (a1 - a0) * beta5 * (theta0 + theta2 * a0  - m)
  rNDE <- CDE + rINTREF

  rPIE <- (a1 - a0) * theta2 * (beta4 + beta5 * a0)
  rINTMED <- theta2 * beta5 * (a1 - a0)^2
  rNIE <- rPIE + rINTMED

  rATE <- rNDE + rNIE

  # output
  out <- list(NULL)
  out$twocomp <- setNames(c(rNDE, rNIE, rATE), c("rNDE", "rNIE", "rATE"))
  out$fourcomp = setNames(c(CDE, rINTREF, rPIE, rINTMED, rATE),
                          c("CDE", "rINTREF", "rPIE", "rINTMED", "rATE"))
  out
}
