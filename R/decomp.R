#' Causal Effect Decomposition Based on a Fitted \code{rwrmed} Model
#'
#' \code{decomp} is a function that implements causal effect decomposition based on a fitted
#'   \code{rwrmed} model. It returns a two-component decomposition of the total effect into
#'   the randomized interventional analogues of the natural direct effect (rNDE) and the natural
#'   indirect effect (rNIE). It also returns a four-component decomposition of the total effect into
#'   the controlled direct effect (CDE) and the randomized analogues of the reference interaction
#'   effect (rINTREF), the mediated interaction effect (rINTMED), and the pure indirect effect (rPIE).
#'
#' @param object An object of class \code{rwrmed}.
#' @param a0 The baseline level of treatment.
#' @param a1 The level of treatment to be contrasted with the baseline.
#' @param m The level of the mediator at which the CDE is evaluated.
#' @param bootstrap Whether to compute standard errors and 95% confidence intervals using the
#'   nonparametric bootstrap.
#' @param rep Number of bootstrap replications if \code{bootstrap = TRUE}. Default is 250.
#'
#' @return A list of two elements.
#'  \item{twocomp}{Two component decomposition of the rATE into rNDE and rNIE.}
#'  \item{fourcomp}{Four component decomposition of the rATE into CDE, rINTREF, rPIE, and rINTMED.}
#' @import stats
#' @importFrom pryr partial
#' @export
#' @seealso \code{\link{rwrmed}} for implementing the regression-with-residuals (RWR)
#'   approach to causal mediation.
#'
#' a0 = 0; a1 = 1; m = 0; bootstrap = TRUE; rep = 250
decomp <- function(object, a0 = 0, a1 = 1, m = 0, bootstrap = TRUE, rep = 250){

  var_names <- object$var_names

  # outcome model coefs

  coefs_y <- coef(object$y_model)

  beta2 <- coefs_y[var_names$treatment]
  beta4 <- coefs_y[var_names$mediator]
  beta5 <- (coefs_y[paste0(var_names$treatment,":",var_names$mediator)] %||%
              coefs_y[paste0(var_names$mediator,":",var_names$treatment)]) %||% 0

  # pred1 and pred0 from the mediator model
  if (object$m_model$family$link != "identity"){

    newdata0 <- newdata1 <- object$data_ed
    newdata0[[var_names$treatment]] <- a0
    newdata1[[var_names$treatment]] <- a1

    pred0 <- weighted.mean(predict.glm(object$m_model, newdata0, type = "response"), w = newdata0$weights)
    pred1 <- weighted.mean(predict.glm(object$m_model, newdata1, type = "response"), w = newdata1$weights)

  } else {

    coefs_m <- coef(object$m_model)

    theta0 <- coefs_m["(Intercept)"]
    theta2 <- coefs_m[var_names$treatment]

    pred0 <- theta0 + theta2 * a0
    pred1 <- theta0 + theta2 * a1
  }

  # effect decomposition
  CDE <- (a1 - a0) * (beta2 + beta5 * m)
  rINTREF <- (a1 - a0) * beta5 * (pred0  - m)
  rNDE <- CDE + rINTREF
  rPIE <- (pred1 - pred0) * (beta4 + beta5 * a0)
  rINTMED <- (pred1 - pred0) * beta5 * (a1 - a0)
  rNIE <- rPIE + rINTMED
  rATE <- rNDE + rNIE

  # prepare output
  est <- se <- lower <- upper <-
    setNames(rep(NA, 7), c("CDE", "rINTREF", "rNDE", "rPIE", "rINTMED", "rNIE", "rATE"))
  est[] <- c(CDE, rINTREF, rNDE, rPIE, rINTMED, rNIE, rATE)

  # boostrap standard errors
  if (bootstrap == TRUE){

    holder <- `colnames<-`(matrix(NA, rep, 7), c("CDE", "rINTREF", "rNDE", "rPIE", "rINTMED", "rNIE", "rATE"))

    y_form <- formula(object$y_model)
    m_form <- formula(object$m_model)
    m_family <- object$m_model$family
    z_forms <- lapply(object$zmodels, formula)
    z_families <- lapply(object$zmodels, family)

    nz <- length(object$zmodels)

    for (k in seq(1, rep))
    {
      if(k %% 10 == 0) cat(".")

      # resample from the original data
      indices <- sample.int(nrow(object$data), replace = TRUE)
      data <- object$data[indices, , drop = FALSE]

      # refit the post-treatment confounder models
      glm_partial <- partial(glm, data = data, weights = weights)
      zmodels <- Map(glm_partial, z_forms, z_families)

      # take pre- and post-treatment confounders
      x <- data[, var_names$pre_cov, drop = FALSE]
      z <- data[, var_names$post_cov, drop = FALSE]

      # copy data to date_ed
      data_ed <- data

      # demean x and residualize z
      for(i in seq_along(x)) data_ed[[names(x)[i]]] <- demean(x[[i]], data_ed$weights)
      for(i in seq_along(z)) data_ed[[names(z)[i]]] <- z[[i]] - zmodels[[i]][["fitted.values"]]

      # mediator and outcome models
      m_model <- glm(formula = m_form, family = m_family, data = data_ed, weights = weights)
      y_model <- lm(formula = y_form, data = data_ed, weights = weights)

      # outcome model coefs
      coefs_y <- coef(y_model)

      beta2 <- coefs_y[var_names$treatment]
      beta4 <- coefs_y[var_names$mediator]
      beta5 <- (coefs_y[paste0(var_names$treatment,":",var_names$mediator)] %||%
                  coefs_y[paste0(var_names$mediator,":",var_names$treatment)]) %||% 0

      # pred1 and pred0 from the mediator model

      if (m_family$link != "identity"){

        newdata0 <- newdata1 <- data_ed
        newdata0[[var_names$treatment]] <- a0
        newdata1[[var_names$treatment]] <- a1

        pred0 <- weighted.mean(predict.glm(m_model, newdata0, type = "response"), w = newdata0$weights)
        pred1 <- weighted.mean(predict.glm(m_model, newdata1, type = "response"), w = newdata1$weights)

      } else {

        coefs_m <- coef(m_model)

        theta0 <- coefs_m["(Intercept)"]
        theta2 <- coefs_m[var_names$treatment]

        pred0 <- theta0 + theta2 * a0
        pred1 <- theta0 + theta2 * a1
      }

      # effect decomposition
      holder[k, "CDE"] <- (a1 - a0) * (beta2 + beta5 * m)
      holder[k, "rINTREF"] <- (a1 - a0) * beta5 * (pred0  - m)

      holder[k, "rPIE"] <- (pred1 - pred0) * (beta4 + beta5 * a0)
      holder[k, "rINTMED"] <- (pred1 - pred0) * beta5 * (a1 - a0)
    }
    holder[, "rNDE"] <- holder[, "CDE"] + holder[, "rINTREF"]
    holder[, "rNIE"] <- holder[, "rPIE"] + holder[, "rINTMED"]
    holder[, "rATE"] <- holder[, "rNDE"] + holder[, "rNIE"]

    se <- apply(holder, 2, sd)
    lower <- apply(holder, 2, quantile, 0.025)
    upper <- apply(holder, 2, quantile, 0.975)
  }

  fullmat <- cbind(est, se, lower, upper)
  colnames(fullmat) <- c("Estimate", "SE", "2.5% Perc", "97.5% Perc")

  out <- NULL
  out$twocomp <- fullmat[c("rNDE", "rNIE", "rATE"), , drop = FALSE]
  out$fourcomp <- fullmat[c("CDE", "rINTREF", "rPIE", "rINTMED", "rATE"), , drop = FALSE]

  out
}
