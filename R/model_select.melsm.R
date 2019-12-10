#' @title Select Model for \code{melsm} Objects
#' @description Model selection for mixed effects location scale models.
#' @name model_select.melsm
#' @param object \code{melsm} object
#' @param type model to select (currently the only option is the highest posterior model \code{hpm})
#' @param ... currently ignored
#'
#' @return object of class \code{model_select.melsm}
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
#'# highest posterior model
#' hpm <- model_select(fit, type = "hpm")
#' summary(hpm, cred = 0.95)
#'}
model_select.melsm <- function(object, type = "hpm", ...){

  if(object$mixture == "none"){
    stop("mixture prior distribution required")
  }

  mod_prob <- model_prob(object)
  top_mods <- mod_prob$model_mat[1:2,]
  samps <- as.data.frame(do.call(rbind.data.frame, object$melsm_fit))

  cors_tested <- length( colnames(samps)[grep("^k", colnames(samps))] )
  if(type == "hpm"){

    colnames(samps)[grep("^k", colnames(samps))] <- colnames(top_mods)[3:(2 + cors_tested)]
    mod_prob <- model_prob(object)
    top_model <- mod_prob$model_mat[1:2,]
    str <- paste("samps[", paste("samps$`", colnames(top_model)[3:(2 + cors_tested)], "`", " == ",
                                 top_model[1, 3:(2 + cors_tested)],
                                 collapse = " & ", sep = ""), ",]")
    hpm_posterior <- eval(parse(text = str))

  } else{
    stop("currently only type = 'hpm' (highest posterior model) is allowed")
  }
  returned_object <- list(hpm_posterior = hpm_posterior,
                          type = type,
                          object = object)
  class(returned_object) <- "model_select.melsm"
  returned_object
}



#' @title S3 \code{model_select} method
#' @param ... currently ignored
#' @seealso  \code{\link{model_select.melsm}}
#' @export
model_select <- function(...){
  UseMethod("model_select")
}

#' Summary Method for \code{model_select.melsm} Objects
#'
#' @param object \code{model_select.melsm} object
#' @param cred credible interval
#' @param ... currently ignored
#' @seealso  \code{\link{model_select.melsm}}
#' @export

summary.model_select.melsm <- function(object, cred = 0.95,...){
  melsm <- object$object
  type <- object$type
  samps <- do.call(rbind.data.frame, melsm$melsm_fit)
  melsm$melsm_fit <- list(object$hpm_posterior)
  colnames(melsm$melsm_fit[[1]]) <- colnames(samps)
  object <- melsm
  cat(paste("hypMuVar: Bayesian Hypothesis Testing",
            "\n          of Mean--Variance Relations\n"))
  cat("Model: MELSM\n")
  cat("Chains:", object$chains, "\n")
  cat("Samples:", (nrow(melsm$melsm_fit[[1]])) , "\n")
  cat("Mixture:", object$mixture, "\n")
  if(object$mixture == "none"){
    cat("Rho Test:", "none\n")
  } else {
    cat("Rho Test:", object$rho_test, "\n")
  }
  cat("Credible Interval:", cred, "\n")

  if(type == "hpm"){
    cat("Selected: Highest Posterior Model\n")
  }  else {
    cat("Selected: Median Posterior Model\n")

  }
  cat("----\n")
  cat("Random Effects Correlations\n")
  print(rho_helper(object, cred = cred), right = F, row.names= F)
  cat("\n----\n")
  cat("\nRandom Effects Standard Deviations\n")
  print(tau_helper(object, cred=cred ), right = F, row.names = F)
  cat("\n----\n")
  cat("\nFixed Effects\n")
  temp <- rbind.data.frame(beta_helper(object, cred = cred), eta_helper(object, cred = cred))
  print(temp, row.names = F, right = F)
  cat("\n----\n")
  if(object$mixture == "SSVS"){
    if(object$k == 2){
      cat("note \nk1: 'spike' \nk2: slab")
    } else {
      cat("note \nk1: 'spike' \nk2: positive slab \nk3: negative slab")
    }
  }
}

#' @title Print Method for \code{model_select.melsm} objects
#'
#' @param x \code{model_select.melsm} object
#' @param ... currently ignored
#' @seealso  \code{\link{model_select.melsm}}
#' @export
print.model_select.melsm <- function(x, ...){
  object <- x
  melsm <- object$object
  type <- object$type
  samps <- do.call(rbind.data.frame, melsm$melsm_fit)
  melsm$melsm_fit <- list(object$hpm_posterior)
  colnames(melsm$melsm_fit[[1]]) <- colnames(samps)
  object <- melsm
  cat(paste("hypMuVar: Bayesian Hypothesis Testing",
            "\n          of Mean--Variance Relations\n"))
  cat("Model: MELSM\n")
  cat("Chains:", object$chains, "\n")
  cat("Samples:", (nrow(melsm$melsm_fit[[1]])) , "\n")
  cat("Mixture:", object$mixture, "\n")
  if(object$mixture == "none"){
    cat("Rho Test:", "none\n")
  } else {
    cat("Rho Test:", object$rho_test, "\n")
  }
 if(type == "hpm"){
    cat("Selected: Highest Posterior Model\n")
  }  else {
    cat("Selected: Median Posterior Model\n")

  }
  cat("----\n")
}





