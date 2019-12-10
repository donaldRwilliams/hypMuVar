#' @title Confirmatory Hypothesis Testing for \code{melsm} Objects
#' @name confirm.melsm
#' @description Confirmatory hypotheses testing of scientific expectations.
#'
#' @param object object of class \code{melsm}
#' @param hyp hypothesis to test
#' @param ... currently ignored
#' @return object of class \code{confirm}
#' \itemize{
#' \item \code{BF} Bayes factor
#' \item \code{post_probs} posterior probabilities
#' \item \code{prior_probs} prior probabilities
#' \item \code{model} model type (\code{melsm})
#' \item \code{mixture} mixture type
#' }
#' @export
#'
#' @examples
#'
#' \dontrun{
#' dat <- nlme::Orthodont
#' fit <- melsm(fixed_location = distance ~  age,
#'             random_location = ~ age | Subject,
#'             fixed_scale = sigma ~  1, k = 3,
#'             random_scale = ~ 1 | Subject,
#'             adapt = 5000,
#'             iter = 10000,
#'             rho_test = "all",
#'             mixture = "SSVS",
#'             data = dat)
#'
#'
#'# hypothesis (rho_01 is null and rho_02 is negative )
#' # vs compliment (not h1)
#' hyp <- list(h1 = c("rho_01 = k1", "rho_02 = k3"))
#' confirm(fit, hyp)
#'
#' # hypothesis  (both positive) vs both null
#' hyp <- list(h1 = c("rho_01 = k2", "rho_02 = k2"),
#' h2 = c("rho_01 = k1", "rho_02 = k1"))
#' confirm(fit, hyp)
#'
#'}
confirm.melsm <- function(object, hyp, ...){

  if(class(object) == "melsm"){

    samps <- do.call(rbind.data.frame, object$melsm_fit)
    n_samps <- nrow(samps)
    n_cors <- length(grep("k", colnames(samps)))
    cor_names <- colnames(model_prob(object)$model_mat)[3:(2 + n_cors)]

    rho_string <- unlist(strsplit(unlist(hyp), split = "="))
    rho_string <- rho_string[stringr::str_starts(rho_string, pattern = "rho")]

  if(!all( gsub(" ", "", rho_string) %in% cor_names)){
    stop("correlations not found")
  }
  colnames(samps)[grep("k", colnames(samps))] <- cor_names
  if (length(hyp) > 2) {
    stop("currently only two hypotheses can be tested")
  }
  mixture <- object$mixture
  n_hyps <- length(hyp)
  k <- object$k
  # hyp <- unlist(hyp)
  post_probs <- unlist(lapply(1:n_hyps, function(x) hyp_helper(hyp[[x]],
                                                               mixture = mixture,
                                                               cor_names = cor_names,
                                                               samps = samps,
                                                               n_samps = n_samps)))

  if(length(hyp) == 1){
    # message("testing against the compliment")
    # compliment post prob
    c_post_prob <- 1 - post_probs

    # h1 prior prob
    H1_prior_prob <- (1/k)^length(unlist(hyp))

    # compliment prior prob
    c_prior_prob <- 1 - H1_prior_prob

    # posterior odds
    post_odds <- post_probs / c_post_prob

    # prior odds
    prior_odds <- H1_prior_prob / c_prior_prob

    # Bayes factor
    BF <- post_odds / prior_odds

    prior_probs = data.frame(h1 = H1_prior_prob, hc = 1 - H1_prior_prob)

  } else {

    prior_odds <- ((1/k)^length(hyp[[1]]))/((1/k)^length(hyp[[2]]))

    post_odds <- post_probs[1] / post_probs[2]

    BF <- post_odds / prior_odds

    prior_probs = data.frame(h1 = ((1/k)^length(hyp[[1]])), h2 = ((1/k)^length(hyp[[2]])))
  }

  returned_object <- list(BF = BF,
                          post_probs = post_probs,
                          prior_probs = prior_probs,
                          model = "melsm",
                          hyp = hyp,
                          mixture = mixture)

  }

  class(returned_object) <- "confirm"
  returned_object

}

#' @title S3 \code{confirm} method
#'
#' @description S3 \code{confirm} method
#' @param object \code{melsm} object
#' @param hyp hypothesis to test
#' @param ... currently ignored
#' @seealso \code{\link{confirm}}
#' @export
confirm <- function(object, hyp, ...){
UseMethod("confirm")
}

#' Summary Method for \code{confirm} Object
#'
#' @param object \code{confirm} Object
#' @param ... currently ignored
#'
#' @return object of class \code{summary.confirm}
#' \itemize{
#' \item \code{dat_prob} data frame with probabilities
#' \item \code{BF_12} Bayes factor in favor of H1
#' \item \code{object} \code{melsm} object
#' }
#' @export
#' @seealso \code{\link{confirm}}
summary.confirm <- function(object,...){

  if(length(object$hyp) == 1){
    dat_prob <- data.frame(Hyp = c("H1", "H2"),
                           Post.prob =  round(c(object$post_probs, 1 - object$post_probs),3),
                           Prior.Prob = round(t(object$prior_probs),3))
    dat_summary <- list(dat_prob = dat_prob,
                        BF_12 = object$BF,
                        object = object)

  } else {

    dat_prob <- data.frame(Hyp = c("H1", "H2"),
                           Post.prob =  round(as.matrix(object$post_probs),3),
                           Prior.Prob = round(t(as.matrix(object$prior_probs)),3))

    dat_summary <- list(dat_prob = dat_prob,
                        BF_12 = object$BF,
                        object = object)

  }

  class(dat_summary) <- "summary.confirm"
  dat_summary
}
#' Print Method for \code{summary.confirm} Objects
#'
#' @param x \code{summary.confirm} object
#' @param ... currently ignored
#' @seealso \code{\link{confirm}}
#' @export
print.summary.confirm <- function(x,...){
  cat(paste("hypMuVar: Bayesian Hypothesis Testing",
            "\n          of Mean--Variance Relations\n"))
  cat("Model: MELSM\n")
  cat("Mixture:", x$object$mixture, "\n")
  cat("----\n")
  cat("Hypotheses\n")
  if(length(hyp) == 1){
    cat("H1:", x$object$hyp[[1]], "\n")
    cat("H2:", "compliment\n")
  } else {
    cat("H1:",  x$object$hyp[[1]], "\n")
    cat("H2:",  x$object$hyp[[2]] , "\n")
  }
  cat("\nBF_12 =", round(x$BF_12,3), "\n")
  cat("----\n")
  print(x$dat_prob, row.names = F, digits = 3)
}



#' Print Method for \code{confirm} Objects
#'
#' @param x \code{confirm} object
#' @param ... currently ignored
#' @seealso \code{\link{confirm}}
#' @export
#'
print.confirm <- function(x,...){
 print(summary(x))
}
