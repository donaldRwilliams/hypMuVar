#' Bayesian Model Averaging
#'
#' @param object \code{melsm} object
#' @param top number of models to average
#' @param ... currently ignored
#'
#' @return object of class \code{melsm}
#' @export
#'
#' @examples
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
#'# Bayesian model averaging
#'bma(fit, top = 2)
#'}
bma <- function(object, top = 2, ...){

  # model probability
  prob_mat <- model_prob(object)$model_mat

  # check at leas
  if(nrow(prob_mat) < top){

    stop("only one model was sampled (no need to model average)")
  }

  # model weights
  model_weights <-  prob_mat[1:top,1] / sum(prob_mat[1:top,1])

  # samples
  samps <- do.call(rbind.data.frame, object$melsm_fit)

  # number of correlations
  n_cors <- length(grep("k", colnames(samps)))

  # correlation names
  cor_names <- colnames(prob_mat)[3:(2 + n_cors)]

  # which were testd
  which_k <- grep("k", colnames(samps))

  # original names
  orig_k <- colnames(samps)[grep("k", colnames(samps))]

  # change names
  colnames(samps)[grep("k", colnames(samps))] <- cor_names

  # samples from each model
  samps_mod <- list()

  for(i in 1:top){

    mod_string <- paste("samps$`",
                        colnames(prob_mat)[3:(2+n_cors)],
                        "` == ",   prob_mat[i,3:(2+n_cors)],
                        sep = "", collapse = " & ")
    mod_string <- paste("samps[", mod_string, ",]")
    samps_mod[[i]] <- eval(parse(text=mod_string))
    colnames(samps_mod[[i]])[which_k] <- orig_k

    }

  # rename for compatability
  samps_mod <- lapply(1:top, function(x) {colnames(samps_mod[[x]])[which_k] <-  orig_k;
                              samps_mod[[x]]})

  # overwrite old melsm fit object
  object$melsm_fit <- samps_mod

  # bma is TRUE
  object$bma <- TRUE

  # bma samples
  object$bma_samps <- sum(prob_mat[1:top,1] * nrow(samps))

  # number of models
  object$top <- top

  object
}



