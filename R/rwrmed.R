#' Regression-with-Residuals Analysis of Causal Mediation
#'
#' \code{rwrmed} is a function that implements the regression-with-residuals (RWR)
#'   approach to causal mediation, allowing for post-treatment confounding of the mediator-outcome
#'   relationship. Specifically, it fits user-specified mediator and outcome models with a set of
#'   residualized post-treatment confounders. It returns an object of class \code{rwrmed}, which can be
#'   used for effect decomposition via the \code{\link{decomp}} function.
#'
#' @param treatment A character string indicating the name of the treatment variable, which
#'   can be either binary or continuous.
#' @param pre_cov A character string vector indicating the names of the pretreatment covariates
#'   that are to be centered at their means. Thus factor-valued covariates should be pre-coded as
#'   dummy variables.
#' @param zmodels A list of fitted \code{lm} or \code{glm} objects for the
#'   post-treatment confounders of the mediator-outcome relationship.
#' @param y_form Formula for the outcome model.
#' @param m_form Formula for the mediator model.
#' @param m_family The family of the mediator model to be specified in \code{\link[stats]{glm}}.
#' @param weights An optional vector of weights to be used in fitting the outcome and mediator models.
#' @param data A data frame containing the variables in the model.
#'
#' @return An object of class \code{rwrmed}.
#'  \item{y_model}{The fitted outcome model.}
#'  \item{m_model}{The fitted mediator model.}
#'  \item{var_names}{Names of the treatment variable, the mediator, the outcome, and
#'     the pretreatment and posttreatment covariates}
#'  \item{data_ed}{The data frame purged of observations with missing variables.}
#'  \item{call}{The matched call.}
#' @import stats
#' @export
#' @seealso \code{\link{decomp}} for effect decomposition based on a fitted \code{rwrmed}
#'  object.
#'
#' @examples
#'
#' ## an example with a continuous mediator ##
#'
#' # treatment and pre-treatment covariates
#' treatment <- "tone_eth"
#' pre_cov <- c("ppage", "female", "hs", "sc", "ba", "ppincimp")
#'
#' # outcome and mediator models
#' y_form <- immigr ~ ppage + female + hs + sc + ba + ppincimp +
#'  tone_eth + emo + tone_eth * emo + p_harm
#' m_form <- emo ~ ppage + female + hs + sc + ba + ppincimp + tone_eth
#'
#' # model for a post-treatment covariate
#' m1 <- lm(p_harm ~ ppage + female + hs + sc + ba + ppincimp + tone_eth,
#'  data = immigration)
#'
#' # effect decomposition using RWR
#' fit1 <- rwrmed(treatment = treatment, pre_cov = pre_cov, zmodels = list(m1),
#'   y_form = y_form, m_form = m_form, data = immigration)
#'
#' out1 <- decomp(fit1)
#'
#' ## an example with a binary mediator ##
#'
#' # treatment and pre-treatment covariates
#' treatment <- "democ"
#' pre_cov <- c("ally", "trade", "h1", "i1", "p1", "e1", "r1", "male", "white", "age", "ed4")
#'
#' # outcome and mediator models
#' y_form <- strike ~ ally + trade + h1 + i1 + p1 + e1 + r1 + male + white + age + ed4 + democ +
#'   immoral + democ * immoral + threatc + cost + successc
#' m_form <- immoral ~ ally + trade + h1 + i1 + p1 + e1 + r1 + male + white + age + ed4 + democ
#'
#' # models for post-treatment covariates
#' m1 <- lm(threatc ~ ally + trade + h1 + i1 + p1 + e1 + r1 + male + white + age + ed4 + democ,
#'   data = peace)
#' m2 <- lm(cost ~ ally + trade + h1 + i1 + p1 + e1 + r1 + male + white + age + ed4 + democ,
#'   data = peace)
#' m3 <- lm(successc ~ ally + trade + h1 + i1 + p1 + e1 + r1 + male + white + age + ed4 + democ,
#'   data = peace)
#'
#' # effect decomposition using RWR
#'
#' fit2 <- rwrmed(treatment = treatment, pre_cov = pre_cov, zmodels = list(m1, m2, m3),
#'   y_form = y_form, m_form = m_form, m_family = binomial("logit"), data = peace)
#'
#' summary(fit2$m_model)
#'
#' out2 <- decomp(fit2)
#'
#'
#' @references Wodtke, Geoffrey T. and Xiang Zhou. 2019. "Effect Decomposition in the Presence of
#'   Treatment-induced Confounding: A Regression-with-Residuals Approach."
#' @references  Zhou, Xiang and Geoffrey T. Wodtke. 2019. "A Regression-with-Residuals Method for Estimating
#'   Controlled Direct Effects." Political Analysis.
#' @references  VanderWeele, Tyler J, Stijn Vansteelandt and James M Robins. 2014. "Effect Decomposition in
#'   the Presence of an Exposure-induced Mediator-outcome Confounder." Epidemiology 25:300-306.
#' @references  VanderWeele, Tyler J. 2014. "A Unification of Mediation and Interaction: a Four-Way
#'   Decomposition." Epidemiology 25(5):749-761.
#'
rwrmed <- function(treatment, pre_cov, zmodels, y_form, m_form,
                   m_family = gaussian, weights, data){

  # match call
  cl <- match.call()

  # check data
  if(!is.data.frame(data)) stop("data must be a data.frame.")
  n <- nrow(data)

  # check treatment
  if(missing(treatment)) stop("treatment must be provided.")
  if(!is.character(treatment)) stop("treatment must be a character string.")

  # check pretreatment covariates
  if(missing(pre_cov)) pre_cov <- character(0L)
  if(!is.character(pre_cov)) stop("pre_cov must be a character string vector.")

  # check zmodels
  if(missing(zmodels)) zmodels <- list(NULL)
  if(!is.list(zmodels)) stop("zmodels must be a list.")
  if(!all(unlist(lapply(zmodels, inherits, "lm")))){
    stop("Each element of zmodels must be an object of class `glm` or `lm`")
  }

  # check missing formulas
  if(missing(y_form)) stop("y_form must be provided.")
  if(missing(m_form)) stop("m_form must be provided.")

  # check weights
  if(missing(weights)) weights <- rep(1, n)
  if((!is.numeric(weights)) || (length(weights) != n)) stop("weights must be a numeric vector of length nrow(data)")

  # get mediator, outcome, and post_cov names
  mediator <- all.vars(m_form)[1L]
  outcome <- all.vars(y_form)[1L]
  post_cov <- vapply(zmodels, function(x) names(x[["model"]])[1], character(1L))

  # get treatment, mediator, outcome, and pre_cov and post_cov
  a <- data[[treatment]]
  m <- data[[mediator]]
  y <- data[[outcome]]
  x <- data[, pre_cov, drop = FALSE]
  z <- data[, post_cov, drop = FALSE]

  # delete rows with missing data
  badRow <- is.na(a) | is.infinite(a) | is.na(m) | is.infinite(m) | is.na(y) | is.infinite(y)
  badRow <- badRow | apply(x, 1, function(v) any(is.na(v) | is.infinite(v)))
  badRow <- badRow | apply(z, 1, function(v) any(is.na(v) | is.infinite(v)))

  # cleaned variables
  data_ed <- data[!badRow, , drop = FALSE]
  n <- nrow(data_ed)

  a <- data_ed[[treatment]]
  m <- data_ed[[mediator]]
  y <- data_ed[[outcome]]
  x <- data_ed[, pre_cov, drop = FALSE]
  z <- data_ed[, post_cov, drop = FALSE]
  data_ed$weights <- weights[!badRow];

  # demean x and residualize z

  for(i in seq_along(x)) data_ed[[names(x)[i]]] <- demean(x[[i]], data_ed$weights)
  for(i in seq_along(z)) data_ed[[names(z)[i]]] <- z[[i]] - zmodels[[i]][["fitted.values"]]

  # fit the mediator and outcome models
  m_model <- glm(formula = m_form, family = m_family, data = data_ed, weights = weights)
  y_model <- lm(formula = y_form, data = data_ed, weights = weights)

  out <- list(y_model = y_model, m_model = m_model,
              var_names = list(treatment = treatment,
                               mediator = mediator,
                               outcome = outcome,
                               pre_cov = pre_cov,
                               post_cov = post_cov),
              data_ed = data_ed, call = cl)
  class(out) = c("rwrmed", "list")
  out
}


