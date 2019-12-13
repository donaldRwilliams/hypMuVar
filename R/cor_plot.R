#' Plot Coefficient Scatteplots
#'
#' @param object \code{melsm} object
#' @param palette paletter for 2d density heat map
#' @param direction direction of the gradient
#' @param ... currently ignored
#'
#' @return list of \code{ggplot} objects
#' @export
#'
#' @examples
#'
#'
#'
#' \dontrun{
#' fit_int <- melsm(fixed_location = rt ~ 1,
#'                  random_location = ~ 1| ID,
#'                  fixed_scale = sigma ~ 1,
#'                  random_scale = ~ 1 | ID, k = 3,
#'                  mixture = "SSVS", adapt = 5000,
#'                  iter = 10000, data = stroop)
#'
#' cor_plot(fit_int)
#' }
cor_plot <- function(object, palette = "YlOrRd", direction = -1,...){


  coefs <- coef(object)
  coefs_names <- dimnames( coefs$coef_array)[[3]]
  index <- t(combn(1:dim(coefs$coef_array)[3], 2))

  plts <- list()
  for(i in 1:nrow(index)){
    dat_plot <-  data.frame(x = coefs$coef_array[,1,index[i,1]], y =  coefs$coef_array[,1,index[i,2]])
    plts[[i]] <- ggplot(dat_plot, aes(x = x, y= y)) +
      xlab(coefs_names[index[i,1]]) +
      ylab(coefs_names[index[i,2]]) +

      stat_density_2d(aes(fill = ..density..),
                      geom = "raster",
                      contour = FALSE,
                      alpha = .75,
                      show.legend = F) +
      # spectral gradient
      scale_fill_distiller(palette= palette,
                           direction=direction)  +
      geom_point() +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0,0)) +
      geom_smooth(method = "lm",
                  color = "white",
                  se = FALSE)

  }

  plts
}
