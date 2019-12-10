#' @title Extract \code{melsm} Coefficients
#'
#' @description Extract coefficients for each cluster (e.g., person, etc.), which obtained by adding the fixed effect
#' to each random effect.
#'
#' @inheritParams ranef.melsm
#'
#' @note this function differs from \code{nlme},
#' in that additional information (not only the coefficients) is returned in a list. This is used
#' to plot the coefficients (\code{\link{plot.coef.melsm}})
#' @return 3D array of class \code{coef.melsm}
#' @export
#' @importFrom stats coef
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
#' coef(fit, cred = 0.90)
#'}
coef.melsm <- function(object, cred = 0.95, ...){

  # cred lower bound
  lb <- (1 - cred) / 2

  # cred upper bound
  ub <- 1 - lb

  post_samps <- do.call(rbind.data.frame, object$melsm_fit)


  fe_scl <- as.matrix(post_samps[,grep(pattern = "^eta", x = colnames(post_samps))])


  fe_loc <- as.matrix(post_samps[,grep(pattern = "^beta", x = colnames(post_samps))])

  res <- post_samps[,grep(pattern = "^u", x = colnames(post_samps))]
  J <- object$melsm_initialize$data()$J
  loc_res <- ncol(object$model_matrices$re_loc_mm)
  scl_res  <- ncol(object$model_matrices$re_scl_mm)
  total_res <- loc_res + scl_res

  fe_loc_ind <- which(colnames(object$model_matrices$fe_loc_mm) ==
                        colnames(object$model_matrices$re_loc_mm))
  fe_scl_ind <- which(colnames(object$model_matrices$fe_scl_mm) ==
                        colnames(object$model_matrices$re_scl_mm))

  fe_ind <- c(fe_loc_ind, fe_scl_ind)



  coef_array <- array(0, c(J,  4, total_res))


  if(is.null(object$dat[, object$cluster])){

    dimnames(coef_array)[[1]] <-  1:length(unique(object$dat[, object$cluster]))

    } else {

      dimnames(coef_array)[[1]] <- levels((object$dat[, object$cluster]))

      }

  dimnames(coef_array)[[2]] <- c("Estimate", "Est.Error", "cred.lb", "cred.ub")

  dimnames(coef_array)[[3]] <- c(paste("location_",
                                       colnames(object$model_matrices$re_loc_mm),
                                       sep = ""),
                                 paste("scale_",
                                       colnames(object$model_matrices$re_scl_mm),
                                       sep = ""))

  for(i in 1:loc_res){
    coef_temp <- res[,  paste("u[",1:J,",", i, "]", sep = "")] + fe_loc[, fe_ind[i]]
    coef_array[,1,i] <- apply(coef_temp, 2, mean)
    coef_array[,2,i] <- apply(coef_temp, 2, sd)
    coef_array[,3:4,i] <- t(apply(coef_temp, 2, quantile, c(lb, ub)))
  }

  for(i in (loc_res+1):total_res){
    coef_temp <- res[,  paste("u[",1:J,",", i, "]", sep = "")] + fe_scl[, fe_ind[i] ]
    coef_array[,1,i] <- apply(coef_temp, 2, mean)
    coef_array[,2,i] <- apply(coef_temp, 2, sd)
    coef_array[,3:4,i] <- t(apply(coef_temp, 2, quantile, c(lb, ub)))

  }

  returned_object <- list(coef_array = coef_array,
                          fe_scl = fe_scl,
                          fe_loc = fe_loc)
  class(returned_object) <- "coef.melsm"
  returned_object

}

#' Plot \code{coef.melsm} Objects
#' @name plot.coef.melsm
#' @param x \code{coef.melsm} object
#' @param ... currently ignored
#' @return list of \code{ggplot} objects
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
#' coefs <- coef(fit, cred = 0.90)
#'
#' plot(coefs)
#'}
plot.coef.melsm <- function(x,...){

  coefs <- x$coef_array

  dim_names <- dimnames(coefs)

  re_names <- dim_names[[3]]

  fe_est <- colMeans(cbind(x$fe_loc, x$fe_scl))

  plot_list <- list()

  for(i in 1:length(re_names)) {

    dat_plot <- coefs[,,re_names[i]]

    dat_plot <- as.data.frame( dat_plot[order(dat_plot[,1]),] )

    dat_plot$ordered_index <- as.factor(1:nrow(dat_plot))

    dat_plot$sig <- as.factor(with(dat_plot, ifelse(cred.lb < fe_est[i] &
                                                    cred.ub > fe_est[i], 0, 1)))

    plt <- ggplot(dat_plot, aes(x = ordered_index,
                                y = Estimate)) +
      geom_hline(yintercept = fe_est[i],
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

#' Print Method for \code{coef.melsm} Objects
#'
#' @param x \code{coef.melsm} object
#' @param ... currently ignored
#' @export
#'
print.coef.melsm <- function(x, ...){
  print(x$coef_array)

}




