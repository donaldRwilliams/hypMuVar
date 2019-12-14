#' Compute Model Probabilities for \code{melsm} Objects
#'
#' @param object \code{melsm} object
#' @param ... currently ignored
#'
#' @return list of class \code{model_prob}
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
#' model_prob(fit)
#' }
model_prob <- function(object, ...){

  if(object$mixture == "none"){

    stop("mixture is required to compute model probabilities")

  }
  # samples
  samps <- do.call(rbind.data.frame, object$melsm_fit)

  # indicators
  k_indicators <- as.matrix(samps[,grep(pattern = "^k", x = colnames(samps))])

  # number of samples
  n_samps <- nrow(samps)

  # models visited
  mods_visited <- unlist(lapply(1:n_samps, function(x) paste(k_indicators[x,], collapse = "")))

  # probs in tbale
  tab_prob <- table(mods_visited)

  # post probs
  post_probs <- tab_prob / n_samps
  post_probs <- post_probs[order(post_probs, decreasing = T)]

  if(!class(post_probs) == "table"){
    post_probs <- as.table(post_probs)
}

  # number of models visited
  n_mods_visited <- dim(post_probs)

  # number of tested correlations
  n_cors <- ncol(k_indicators)
  rho_sum <- hypMuVar:::rho_helper(object, cred = 0.95)
  rho_sum[rho_sum==""]<-NA
  rho_sum <- na.omit(rho_sum)
  rho_split <- strsplit(as.character(rho_sum[,1]), split = ":")
  model_mat <- matrix(0, nrow = dim(post_probs) , ncol = 2 + n_cors)
  colnames(model_mat) <-  c("Post.prob", "BF.1i", unlist(lapply(1:n_cors, function(x) rho_split[[x]][1])))

  for(i in 1:n_mods_visited){
    model_mat[i,3:(2+n_cors)] <- as.numeric(strsplit(names(post_probs), split = "")[[i]])
    model_mat[i,1] <- round(post_probs[i],3)
    model_mat[i,2] <- round((post_probs[1] / post_probs)[i], 3)
  }

  row.names(model_mat) <- paste("Model", 1:n_mods_visited, sep = " ")

  returned_object <- list(model_mat = model_mat,
                          mixture = object$mixture,
                          rho_test = object$rho_test,
                          k = object$k,
                          rho_summary = rho_sum)

  class(returned_object) <- "model_prob"
  returned_object
}


#' Print Method for \code{model_prob} Objects
#'
#' @param x \code{model_prob} object
#' @param ... currently ignored
#' @export
print.model_prob <- function(x, ...){
  object <- x
  cat(paste("hypMuVar: Bayesian Hypothesis Testing",
            "\n          of Mean--Variance Relations\n"))
  cat("Model: MELSM\n")
  cat("Mixture:", object$mixture, "\n")
  cat("Rho Test:", object$rho_test, "\n")
  cat("Components:", object$k, "\n")
  cat("----\n\n")
  print(object$model_mat,right = T, ...)
  cat("----\n\n")
  cat(paste(object$rho_summary[,1],collapse =   "\n"))
}

