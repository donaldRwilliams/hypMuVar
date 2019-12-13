#' @title Extract \code{melsm} Random Effects
#' @name ranef.melsm
#' @param object \code{melsm} object
#' @param cred credible interval
#' @param ... currently ignored
#'
#' @return 3D array of class \code{ranef.melsm}
#' @note this function differs from \code{nlme},
#' in that additional information (not only the random effects) is returned in a list. This is used
#' to plot the random effects (\code{\link{plot.ranef.melsm}})
#' @export
#' @importFrom nlme ranef
#' @export ranef
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
#' ranef(fit, cred = 0.90)
#'}
ranef.melsm <- function(object, cred = 0.95,...){

  # cred lower bound
  lb <- (1 - cred) / 2

  # cred upper bound
  ub <- 1 - lb

  # posterior samples
  post_samps <- do.call(rbind.data.frame, object$melsm_fit)

  # random effects
  res <- post_samps[,grep(pattern = "^u", x = colnames(post_samps))]

  # number of clusters
  J <- object$melsm_initialize$data()$J

  # number of location random effects
  n_loc_res <- ncol(object$model_matrices$re_loc_mm)

  # number of scale random effects
  n_scl_res  <- ncol(object$model_matrices$re_scl_mm)

  # total number of random effects
  total_res <- n_loc_res + n_scl_res

  # random effects array
  re_array <- array(0, c(J,  4, total_res))


  if(is.null(object$dat[, object$cluster])){

    dimnames(re_array)[[1]] <- 1:length(unique(object$dat[, object$cluster]))

    } else {

  # row names
  dimnames(re_array)[[1]] <- levels((object$dat[, object$cluster]))

  }
  # colnames
  dimnames(re_array)[[2]] <- c("Estimate", "Est.Error", "cred.lb", "cred.ub")

  # dimension names
  dimnames(re_array)[[3]] <- c(paste("location_",
                                     colnames(object$model_matrices$re_loc_mm),
                                     sep = ""),
                               paste("scale_",
                                     colnames(object$model_matrices$re_scl_mm),
                                     sep = ""))

  for(i in 1:n_loc_res){

    re_temp <- res[,  paste("u[",1:J,",", i, "]", sep = "")]
    re_array[,1,i] <- apply(re_temp, 2, mean)
    re_array[,2,i] <- apply(re_temp, 2, sd)
    re_array[,3:4,i] <- t(apply(re_temp, 2, quantile, c(lb, ub)))

    }

  for(i in (n_loc_res+1):total_res){

    re_temp <- res[,  paste("u[",1:J,",", i, "]", sep = "")]
    re_array[,1,i] <- apply(re_temp, 2, mean)
    re_array[,2,i] <- apply(re_temp, 2, sd)
    re_array[,3:4,i] <- t(apply(re_temp, 2, quantile, c(lb, ub)))

    }

  returned_object <- list(re_array = re_array)
  class(returned_object) <- "ranef.melsm"
  returned_object

}




#' Plot \code{ranef.melsm} Objects
#' @name plot.ranef.melsm
#' @param x \code{ranef.melsm} object
#' @param ... currently ignored
#' @return list of \code{ggplot} objects
#' @import ggplot2
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
#' res <- ranef(fit, cred = 0.90)
#'
#' plot(res)
#'}

plot.ranef.melsm <- function(x,...){

  coefs <- x$re_array

  dim_names <- dimnames(coefs)

  re_names <- dim_names[[3]]

  # fe_est <- colMeans(cbind(x$fe_loc, x$fe_scl))

  plot_list <- list()

  for(i in 1:length(re_names)) {

    dat_plot <- coefs[,,re_names[i]]

    dat_plot <- as.data.frame( dat_plot[order(dat_plot[,1]),] )

    dat_plot$ordered_index <- as.factor(1:nrow(dat_plot))

    dat_plot$sig <- as.factor(with(dat_plot, ifelse(cred.lb < 0 &
                                                      cred.ub > 0, 0, 1)))

    plt <- ggplot(dat_plot, aes(x = ordered_index,
                                y = Estimate)) +
      geom_hline(yintercept = 0,
                 linetype = "dotted") +
      geom_errorbar(aes(ymax = cred.ub,
                        ymin = cred.lb, color = sig),
                    width = 0) +
      geom_point() +
      ylab(re_names[i]) +
      theme(legend.position = "none")

    plot_list[[i]] <- plt
  }

  plot_list

}

#' Print Method for \code{ranef.melsm} Objects
#'
#' @param x \code{ranef.melsm} object
#' @param ... currently ignored
#' @export
#'
print.ranef.melsm <- function(x,...){
  print(x$re_array)

}


