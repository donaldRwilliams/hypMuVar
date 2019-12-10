#' Prior Names
#'
#' @param fixed_location a \code{formula} object for the fixed-effects part of the location sub-model, with the response (on the left) and
#'                       the terms (on the right) sepearted by \code{~}.
#' @param random_location a \code{formula} object for the random effects part of the location sub-model, with the random effect (on the left)
#'                        and the cluster (on the right) seperated by \code{~}.
#' @param fixed_scale a \code{formula} object for the fixed-effects part of the scale sub-model, with the response (on the left) and
#'                       the terms (on the right) sepearted by \code{~}.
#' @param random_scale a \code{formula} object for the random effects part of the location sub-model, with the random effect (on the left)
#'                      and the cluster (on the right) seperated by \code{~}.
#' @param data a data frame  containing the variables named in \code{fixed_location}, \code{random_location},
#'             \code{fixed_scale}, and \code{random_scale}
#'
#' @export
#'
#' @examples
#'
#' dat <- nlme::Orthodont
#' prior_names(fixed_location = distance ~ Sex + age,
#'            fixed_scale = sigma ~ age,
#'            random_location = ~ 1 | Subject,
#'            random_scale = ~ 1 | Subject,
#'            data = dat)
prior_names <- function(fixed_location,
                        random_location,
                        fixed_scale,
                        random_scale,
                        data){


  # all.vars(fixed_location) %in% colnames(dat)
  fe_loc_mm <- model.matrix(fixed_location, data)
  fe_scl_mm <- model.matrix(as.formula(paste("~",
                                             all.vars(fixed_scale)[-1])),
                            data)

  fe_loc_prior <- paste("beta_",
                        0:(ncol(fe_loc_mm)-1), sep = "")
  fe_loc_mat <- matrix(0, nrow = ncol(fe_loc_mm), ncol = 2)
  fe_loc_mat[1:ncol(fe_loc_mm),1] <- fe_loc_prior
  fe_loc_mat[1:ncol(fe_loc_mm),2] <- paste("location_", colnames(fe_loc_mm), sep = "")
  fe_loc_prior_df <- as.data.frame(fe_loc_mat)
  colnames(fe_loc_prior_df) <- c("prior_name", "parameter")



  fe_scl_prior <- paste("eta_",
                        0:(ncol(fe_scl_mm)-1), sep = "")
  fe_scl_mat <- matrix(0, nrow = ncol(fe_scl_mm), ncol = 2)
  fe_scl_mat[1:ncol(fe_scl_mm),1] <- fe_scl_prior
  fe_scl_mat[1:ncol(fe_scl_mm),2] <- paste("scale_", colnames(fe_scl_mm), sep = "")
  fe_scl_prior_df <- as.data.frame(fe_scl_mat)
  colnames(fe_scl_prior_df) <- c("prior_name", "parameter")


  # split RE formulas
  split_re_loc <- all.vars(random_location)
  split_re_scl <- all.vars(random_scale)

  if(length(split_re_loc) == 1){
    re_loc_mm <- model.matrix(~ 1, data)
  } else {

    re_loc_mm <- model.matrix(as.formula(paste("~", split_re_loc[1])), data)

  }
  if(length(split_re_scl) == 1){

    re_scl_mm <- model.matrix(~ 1, data)

  } else {

    re_scl_mm <- model.matrix(as.formula(paste("~", split_re_scl[1])), data)

  }


  # location model matrix names
  re_loc_names <- paste("location_",
                        colnames(re_loc_mm),
                        sep = "")

  # scale model matrix names
  re_scl_names <- paste("scale_",
                        colnames(re_scl_mm),
                        sep = "")

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


  rho_prior_df <- do.call(rbind.data.frame, strsplit( unlist(cor_names), split = ":"))
  colnames(rho_prior_df) <- c("prior_name", "parameter")


  loc_re_mat <- matrix(0, nrow = ncol(re_loc_mm), ncol = 2)
  loc_re_mat[1:ncol(re_loc_mm),1] <- paste("tau_", 0:(ncol(re_loc_mm)-1), sep = "")
  loc_re_mat[1:ncol(re_loc_mm),2] <- re_loc_names

  loc_re_prior_df <- as.data.frame(loc_re_mat)
  colnames(loc_re_prior_df) <-  c("prior_name", "parameter")





  scl_re_mat <- matrix(0, nrow = ncol(re_scl_mm), ncol = 2)
  scl_re_mat[1:ncol(re_scl_mm),1] <- paste("tau_", ncol(re_loc_mm):(total_res-1), sep = "")
  scl_re_mat[1:ncol(re_scl_mm),2] <- re_scl_names

  scl_re_prior_df <- as.data.frame(scl_re_mat)
  colnames(scl_re_prior_df) <-  c("prior_name", "parameter")

  cat("Correlation\n")
  print(rho_prior_df, right = F, row.names = F )
  cat("\nnote: currently each correlation is assinged same spike and slab prior\n")
  cat("\ne.g., rho = list(slab = 0.5, spike = 0.01)\n")
  cat("----\n")
  cat("Fixed Effects\n\n")
  print(rbind.data.frame(fe_loc_prior_df, fe_scl_prior_df),
        row.names = F, right = F)
  cat("\ne.g., beta = list(beta_1 = 'dnorm(0, 1)')\n")
  cat("----\n")

  cat("Random Effects\n\n")
  print(rbind.data.frame(loc_re_prior_df, scl_re_prior_df),
        row.names = F, right = F)
  cat("\ne.g., tau = list(tau_0 = 'dunif(0, 10)')\n")
  cat("----\n")
}


