globalVariables(c('re_slc_int_form',
                  'ordered_index',
                  'Estimate',
                  'cred.ub',
                  'cred.lb',
                  'sig',
                  'hyp'))

#' @importFrom utils combn
#' @importFrom stringr str_extract
numextract <- function(string){
  stringr::str_extract(string, "\\-*\\d+\\.*\\d*")
}

rho_helper <- function(object, cred = cred,....){

  # cred lower bound
  lb <- (1 - cred) / 2

  # cred upper bound
  ub <- 1 - lb

  # posterior samples
  post_samps <- do.call(rbind.data.frame, object$melsm_fit)

  # random effects cors
  rhos <- as.matrix(post_samps[,grep(pattern = "^rho", x = colnames(post_samps))])

  # posterior means
  rhos_mu <- colMeans(rhos)

  # posterior sds
  rhos_sd <- apply(rhos, 2, sd)

  # credible intervals
  rhos_cred <- apply(rhos, 2, quantile, c(lb, ub))

  # location model matrix names
  re_loc_names <- paste("location_", colnames(object$model_matrices$re_loc_mm),sep = "")

  # scale model matrix names
  re_scl_names <- paste("scale_", colnames(object$model_matrices$re_scl_mm),sep = "")

  # random effects names
  re_names <- c(re_loc_names, re_scl_names)

  # number of random effects
  total_res <- length(re_loc_names) + length(re_scl_names)

  # rho index pairwise combinatinos
  rho_index <- as.matrix(t(combn(1:total_res, 2)))

  # store names
  cor_names <- list()

  # correlation names
  for(i in 1:nrow(rho_index)){

    cor_names[[i]] <- paste("rho_",
                            rho_index[i,1] - 1,
                            rho_index[i,2] - 1, ": ",
                            re_names[rho_index[i,1]], "_" ,
                            re_names[rho_index[i,2]], sep = "")
  }

  # rho summary
  rho_summary <- data.frame(Post.mean = rhos_mu,
                            Post.sd = rhos_sd,
                            Cred.lb = rhos_cred[1,],
                            Cred.ub = rhos_cred[2,])

  # add rho names
  rho_summary <- data.frame(Correlations = unlist(cor_names),
                            round(rho_summary, 3))

  # remove first column names
  colnames(rho_summary)[1] <- ""

  # start: KM mixture
  if(object$mixture == "KM"){

    # indicators
    rhos_k <- post_samps[,grep(pattern = "^k", x = colnames(post_samps))]

    # one correlation
    if(is.null(dim(rhos_k))){

      rhos_prob <- round(colMeans(as.matrix(rhos_k)), 3)
      rho_summary[,6] <-  rhos_prob
      rho_summary[,7] <-  1 - rhos_prob
      colnames(rho_summary)[6] <- "Pr.slab"
      colnames(rho_summary)[7] <- "Pr.spike"

      # more than one correlations
      } else {

      rhos_prob <- round(colMeans(rhos_k), 3)
      rho_summary[as.numeric(numextract(names(rhos_prob))),6] <-  rhos_prob
      rho_summary[as.numeric(numextract(names(rhos_prob))),7] <-  1 - rhos_prob
      colnames(rho_summary)[6] <- "Pr.slab"
      colnames(rho_summary)[7] <- "Pr.spike"
      }
    } # end KM mixture

  # start SSVS mixture
  if (object$mixture == "SSVS") {

    # indicators
    rhos_k <- as.matrix(post_samps[,grep(pattern = "^k", x = colnames(post_samps))])

    # one correlations
    if (dim(rhos_k)[2] == 1) {

      # two components
      if (object$k == 2) {

        rhos_prob_k1 <- round(colMeans(rhos_k == 1), 3)
        rhos_prob_k2 <- round(colMeans(rhos_k == 2), 3)
        rho_summary[,6] <-  rhos_prob_k1
        rho_summary[,7] <-  rhos_prob_k2
        colnames(rho_summary)[6] <- "Pr.k1"
        colnames(rho_summary)[7] <- "Pr.k2"

        # three components
        } else {

          rhos_prob_k1 <- round(colMeans(rhos_k == 1), 3)
          rhos_prob_k2 <- round(colMeans(rhos_k == 2), 3)
          rhos_prob_k3 <- round(colMeans(rhos_k == 3), 3)
          rho_summary[,6] <-  rhos_prob_k1
          rho_summary[,7] <-  rhos_prob_k2
          rho_summary[,8] <-  rhos_prob_k3
          colnames(rho_summary)[6] <- "Pr.k1"
          colnames(rho_summary)[7] <- "Pr.k2"
          colnames(rho_summary)[8] <- "Pr.k3"
          }
      # more than one correlation
      } else {

        # two components
        if(object$k == 2){

          rhos_prob_k1 <- round(colMeans(rhos_k == 1), 3)
          rhos_prob_k2 <- round(colMeans(rhos_k == 2), 3)
          rho_summary[as.numeric(numextract(names(rhos_prob_k1))),6] <-  rhos_prob_k1
          rho_summary[as.numeric(numextract(names(rhos_prob_k1))),7] <-  rhos_prob_k2
          colnames(rho_summary)[6] <- "Pr.k1"
          colnames(rho_summary)[7] <- "Pr.k2"

          # three components
          } else {

            rhos_prob_k1 <- round(colMeans(rhos_k == 1), 3)
            rhos_prob_k2 <- round(colMeans(rhos_k == 2), 3)
            rhos_prob_k3 <- round(colMeans(rhos_k == 3), 3)
            rho_summary[as.numeric(numextract(names(rhos_prob_k1))),6] <-  rhos_prob_k1
            rho_summary[as.numeric(numextract(names(rhos_prob_k1))),7] <-  rhos_prob_k2
            rho_summary[as.numeric(numextract(names(rhos_prob_k1))),8] <-  rhos_prob_k3
            colnames(rho_summary)[6] <- "Pr.k1"
            colnames(rho_summary)[7] <- "Pr.k2"
            colnames(rho_summary)[8] <- "Pr.k3"

          }
      }
    }

  # remove first column name
  rho_summary[is.na(rho_summary)] <- ""

  # rho summary
  rho_summary

}


# extract and summarise random effects SDs
tau_helper <- function(object, cred = cred,...){

  # cred lower bound
  lb <- (1 - cred) / 2

  # cred upper bound
  ub <- 1 - lb

  # posterior samples
  post_samps <- do.call(rbind.data.frame, object$melsm_fit)

  # random effects SDs
  taus <- post_samps[,grep(pattern = "^Tau", x = colnames(post_samps))]

  # set 0s to NA (off-diagonal elements)
  taus[taus == 0] <- NA

  # remove NA
  taus <- t(na.omit(t(taus)))

  # posterior means
  taus_mu <- colMeans(taus)

  # posterior sds
  taus_sd <- apply(taus, 2, sd)

  # credible intervals
  taus_cred <- apply(taus, 2, quantile, c(lb, ub))

  # location RE SD names
  re_loc_names <- paste("location_", colnames(object$model_matrices$re_loc_mm),sep = "")

  # scale RE SD names
  re_scl_names <- paste("scale_", colnames(object$model_matrices$re_scl_mm),sep = "")

  # re names
  re_names <- c(re_loc_names, re_scl_names)

  # RE sd names
  tau_names <- paste("tau_", 0:(length(re_names)-1), ": ", re_names, sep = "")

  # tau summary
  tau_summary <- round(data.frame(
    Post.mean = taus_mu,
    Post.sd = taus_sd,
    Cred.lb = taus_cred[1,],
    Cred.ub = taus_cred[2,]), 3)

  # add RE SD names
  tau_summary <- cbind.data.frame(tau_names, tau_summary)

  # remove first column name
  colnames(tau_summary)[1] <- ""

  # RE sd summary
  tau_summary
}

# extract and summarise location fixed effects
beta_helper <- function(object, cred = cred, ...){

  # cred lower bound
  lb <- (1 - cred) / 2

  # cred upper bound
  ub <- 1 - lb

  # posterior samples
  post_samps <- do.call(rbind.data.frame, object$melsm_fit)

  # location fixed effects
  betas <- as.matrix(post_samps[,grep(pattern = "^beta", x = colnames(post_samps))])

  # posterior means
  betas_mu <- colMeans(betas)

  # posterior sds
  betas_sd <- apply(betas, 2, sd)

  # credible intervals
  betas_cred <- apply(betas, 2, quantile, c(lb, ub))

  # model matrix names
  re_loc_names <- paste("beta_", 0:(ncol(betas)-1), ": ", "location_",
                        colnames(object$model_matrices$fe_loc_mm),sep = "")

  # beta summary
  beta_summary <- round(data.frame(Post.mean = betas_mu,
                                   Post.sd = betas_sd,
                                   Cred.lb = betas_cred[1,],
                                   Cred.ub = betas_cred[2,]), 3)

  # add coefficient names
  beta_summary <- cbind.data.frame(re_loc_names, beta_summary)

  # remove first column name
  colnames(beta_summary)[1] <- ""

  # beta summary
  beta_summary
}

# extract and summarise scale fixed effects
eta_helper <- function(object, cred = cred, ...){
  # cred lower bound
  lb <- (1 - cred) / 2

  # cred upper bound
  ub <- 1 - lb

  # posterior samples
  post_samps <- do.call(rbind.data.frame, object$melsm_fit)

  # scale fixed effects
  etas <- as.matrix(post_samps[,grep(pattern = "^eta", x = colnames(post_samps))])

  # posterior means
  etas_mu <- colMeans(etas)

  # posterior sds
  etas_sd <- apply(etas, 2, sd)

  # credible interval
  etas_cred <- apply(etas, 2, quantile, c(lb, ub))

  # model matrix names
  re_scl_names <- paste("eta_", 0:(ncol(etas)-1), ": ",
                        " scale_", colnames(object$model_matrices$fe_scl_mm),sep = "")

  # eta summary
  eta_summary <- round(data.frame(Post.mean = etas_mu,
                                  Post.sd = etas_sd,
                                  Cred.lb = etas_cred[1,],
                                  Cred.ub = etas_cred[2,]), 3)

  # add coefficient names
  eta_summary <- cbind.data.frame(re_scl_names, eta_summary)

  # remove first column name
  colnames(eta_summary)[1] <- ""

  # eta summary
  eta_summary
}

hyp_helper <- function(x, mixture,
                       cor_names, samps,
                       n_samps){

  temp <- list()
 if(mixture == "SSVS"){

    for(i in 1:length(x)){

      if(length( grep("k1", x = x[i]) ) == 1 ){

        tested_rho <- unlist(stringr::str_extract_all( x[i], cor_names))
        temp[[i]] <- paste("samps$`", tested_rho, "` == 1", sep = "")

      } else if ( length( grep("k2", x = x[i]) ) == 1 ){

        tested_rho <- unlist(stringr::str_extract_all(x[i], cor_names))
        temp[[i]] <- paste("samps$`", tested_rho, "` == 2", sep = "")

      } else if (length( grep("k3", x = x[i]) ) == 1 ){

        tested_rho <- unlist(stringr::str_extract_all( x[i], cor_names))
        temp[[i]] <- paste("samps$`", tested_rho, "` == 3", sep = "")

      } else {

        stop("error in hypothesis formulation")

      }
    }
  } else {

    for(i in 1:length(x)){

      if(length( grep("k1", x = x[i]) ) == 1 ){

        tested_rho <- unlist(stringr::str_extract_all( x[i], cor_names))
        temp[[i]] <- paste("samps$`", tested_rho, "` == 0", sep = "")

      } else if ( length( grep("k2", x = x[i]) ) == 1 ){

        tested_rho <- unlist(stringr::str_extract_all(x[i], cor_names))
        temp[[i]] <- paste("samps$`", tested_rho, "` == 1", sep = "")

      }  else {

        stop("error in hypothesis formulation")

      }
    }
  }

  hyp_string <- paste("samps[", paste(temp, collapse = " & "), ",]", sep = "")
  H1_post_prob <- nrow(eval(parse(text = hyp_string))) / n_samps
  H1_post_prob

}


