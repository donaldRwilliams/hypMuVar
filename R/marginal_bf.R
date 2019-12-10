#' Compute Marginal Bayes Factors for \code{melsm} Objects
#'
#' @param object \code{melsm} object
#' @param H1 hypothesis one (only for \code{mixture = "SSVS"})
#' @param H2 hypothesis two (only for \code{mixture = "SSVS"})
#' @param ... currently ignored
#'
#' @return list of class \code{marginal bf}
#' @export
marginal_bf <- function(object,
                        H1 = "k1",
                        H2 = "compliment", ...){

  # check class
  if (class(object) != "melsm") {

    stop("object must be of class melsm")

  }

  # rho summary
  rho_summary <- rho_helper(object, cred = 0.95)

  # mixture SSVS
  if(object$mixture == "SSVS"){

    # two component
    if(object$k == 2){

      # post probs
      post_probs <- rho_summary[,6:7]

      # H1 post prob
      H1_probs <- post_probs[,grep(H1, x = colnames(post_probs))]

      # compliment
      if(H2 == "compliment"){

        BF_12 <- as.numeric(H1_probs) / as.numeric(post_probs[,-grep(H1, x = colnames(post_probs))])

      } else {

        BF_12 <- as.numeric(H1_probs) /  as.numeric(post_probs[,grep(H2, x = colnames(post_probs))])
      }

      BF_dat <- data.frame(Correlation = rho_summary[,1],
                           BF_12  = BF_12)
    } else if (object$k == 3){

      # post probs
      post_probs <- rho_summary[,6:8]

      # h1 post probs
      H1_probs <- post_probs[,grep(H1, x = colnames(post_probs))]

      # compliment
      if(H2 == "compliment"){

        BF_12 <- (H1_probs / rowSums(post_probs[,-grep(H1, x = colnames(post_probs))])) / 0.5

      } else {

        BF_12 <- H1_probs /  post_probs[,grep(H2, x = colnames(post_probs))]
      }

      BF_dat <- data.frame(Correlation = rho_summary[,1],
                           BF_12  = round(BF_12,3))

    }

    # mixture KM
  } else if (object$mixture == "KM") {

    # post probs
    BF_10 <- as.numeric(rho_summary[,6]) / as.numeric(rho_summary[,7])

    BF_dat <- data.frame(Correlation = rho_summary[,1],
                         BF_10  = round(BF_10, 3))


  } else {

    stop("no mixture for the random effects correlations")

  }

  returned_object <- list(BF_dat = BF_dat,
                          mixture = object$mixture,
                          k = object$k,
                          H1 = H1,
                          H2 = H2,
                          rho_test = object$rho_test)

  class(returned_object) <- "marginal_bf"

  returned_object

}


#' Print Method for \code{marginal_bf} objects
#'
#' @param x \code{marginal_bf} object
#' @param ... currently ignored
#' @export
print.marginal_bf <- function(x, ...){
  object <- x
  cat(paste("hypMuVar: Bayesian Hypothesis Testing",
            "\n          of Mean--Variance Relations\n"))
  cat("Model: MELSM\n")
  cat("Mixture:", object$mixture, "\n")
  cat("Rho Test:", object$rho_test, "\n")
  cat("Hypotheses\n")

  if(object$mixture == "SSVS"){
    cat("H1: rho =", object$H1, "\n")
    cat("H2: rho =", object$H2, "\n")

  } else {

    cat("H0: rho = 0", "\n")
    cat("H1: rho != 0", "\n")

  }
  cat("-----\n")
  print(object$BF_dat, row.names = F, right = F)
}
