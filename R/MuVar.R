#' Title
#'
#' @param y
#' @param cluster
#' @param predictor
#' @param data
#' @param mixture
#' @param type
#' @param tau_sd
#' @param tau_sd_df
#' @param z_slab_sd
#' @param z_spike_sd
#' @param pr_prob_rho
#' @param chains
#' @param iter
#' @param warmup
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
muvar <- function(y,
                  cluster,
                  predictor,
                  data,
                  mixture = FALSE,
                  type = NULL,
                  k = 2,
                  tau_sd = NULL,
                  tau_sd_df = NULL,
                  z_slab_sd = NULL,
                  z_spike_sd = NULL,
                  pr_prob_rho = NULL,
                  chains = 4,
                  iter = 5000,
                  warmup = 1000,
                  ...) {
  # arguments
  arg <- match.call()

  # data frame
  data <- data.frame(na.omit(data))

  # check outcome
  check_y <- eval(arg[[match("y", names(arg))]], envir = data)
  if (is.function(check_y)) {
    stop("y not found in the data")
  }

  # check cluster
  check_cluster <-
    eval(arg[[match("cluster", names(arg))]], envir = data)
  if (is.function(check_cluster)) {
    stop("cluster not found in the data")
  }

  # check location formula
  if (isFALSE(class(predictor) == "formula")) {
    stop("location must be a formula")
  }

  # location terms
  mu_terms <- labels(terms(predictor))

  # scale terms
  sd_terms <- labels(terms(predictor))

  # check # of location terms
  if (length(mu_terms) > 1) {
    stop("currently only one predictor can be included in the model")
  }

  # check # of scale terms
  if (length(sd_terms) > 1) {
    stop("currently only one predictor can be included in the model")
  }

  # number of random effects
  if (length(mu_terms) == 0) {
    mu_REs <- 1
  } else{
    mu_REs <- 2
  }
  if (length(sd_terms) == 0) {
    sd_REs <- 1
  } else{
    sd_REs <- 2
  }

  # total # of REs
  total_REs <- mu_REs + sd_REs

  # total # of RE correlations
  total_cors <- 0.5 * (total_REs * (total_REs - 1))

  # outcome
  y <- as.numeric(eval(arg[[match("y", names(arg))]], envir = data))

  y_scale <- as.numeric(scale(y))

  # cluster
  ID <-
    as.numeric(eval(arg[[match("cluster", names(arg))]], envir = data))

  # mean structure
  mu_formula <- model.matrix(predictor, data)

  # check number of weights
  if (ncol(mu_formula) > 2) {
    stop("currently only an intercept and one (regression) weight can be estimated")
  }

  # location weight names
  mu_weights <- colnames(mu_formula)

  # variance structure
  sd_formula <- model.matrix(predictor, data)

  # check number of weights
  if (ncol(sd_formula) > 2) {
    stop("currently only an intercept and one (regression) weight can be estimated.")
  }

  # scale weight names
  sd_weights <- colnames(sd_formula)

  if (mu_REs == 1) {
    # random intercept
    mu_structure <- paste("yhat[i] <-  u[ID[i], 1]\n")
    var_structure <- paste("shat[i] <- u[ID[i], 2]\n}\n")



    # default tau sd
    if (is.null(tau_sd)) {
      tau11_sd <- sd(y)
      tau22_sd <- sd(y)
    } else {
      # check tau_df length
      if (length(tau_sd) != 2 | !is.list(tau_sd)) {
        stop("two values for the random effects sd required in a list.")
      } else {
        # check names
        if (isFALSE(all(c("tau11_sd", "tau22_sd") %in% names(tau_sd)))) {
          stop("incorrect names in tau_sd")
        } else {
          # user defined prior sd
          tau11_sd <- tau_sd$tau11_sd
          tau22_sd <- tau_sd$tau22_sd
        }
      }
    }

    # default df
    if (is.null(tau_sd_df)) {
      tau11_df <- 10
      tau22_df <- 10
    } else {
      # check tau_df length
      if (length(tau_sd_df) != 2 | !is.list(tau_sd_df)) {
        stop("two values for the random effects sd required in a list.")
      } else{
        # check names
        if (isFALSE(all(c("tau11_sd_df", "tau22_sd_df") %in% names(tau_sd_df)))) {
          stop("incorrect names in tau_sd")
        } else {
          # user defined prior df
          tau11_df <- tau_sd_df$tau11_sd_df
          tau22_df <- tau_sd_df$tau11_sd_df
        }
      }
    }
    if (isFALSE(mixture)) {
      warning("mixture set to FALSE. model estimated with no spike and slab formulation.")

      # default rho
      if (is.null(z_slab_sd)) {
        z_slab_sd <- 0.5
      } else {
        # check tau_df length
        if (!is.list(z_slab_sd)) {
          stop("z sd required in a list.")
        } else{
          # check names
          if (isFALSE(all(c("z_slab_sd") %in% names(z_slab_sd)))) {
            stop("incorrect names in tau_sd")
          } else {
            # user defined prior df
            z_slab_sd <- z_slab_sd$z_slab_sd
          }
        }
      }

      beta0_mean <- mean(y)
      beta0_sd   <- sd(y)

      eta0_mean  <- 0
      eta0_sd    <- 10


      model <- paste(
        likelihood,
        mu_structure,
        var_structure,
        re_structure,
        bhat,
        beta,
        re_sigma,
        rho,
        cov_mat
      )

      model_initialize <- jags.model(
        textConnection(model),
        list(
          y = y_scale,
          D = total_REs,
          total_cors = total_cors,
          J = length(unique(ID)),
          n = length(y),
          ID = ID,
          beta0_mean = beta0_mean,
          beta0_sd = beta0_sd,
          eta0_mean = eta0_mean,
          eta0_sd = eta0_sd,
          tau_11_sd = tau11_sd,
          tau_11_df = tau11_df,
          tau_22_sd = tau22_sd,
          tau_22_df = tau22_df,
          z_sd = z_slab_sd
        ),
        n.chains = chains,
        n.adapt = warmup
      )

      samps <- coda.samples(
        model_initialize,
        variable.names = c("beta0", "eta0",
                           "tau", "rho", "u"),
        n.iter = iter
      )

      res <-  list(
        model = model,
        y = y,
        ID = ID,
        samps = samps,
        model = model,
        model_initialize = model_initialize
      )

    } else {
      # check mixture type
      if (isFALSE(any(c("KM", "SVSS") %in% type))) {
        stop("type must be KM or SVSS")
      } else {
        # Kuo & Mallick
        if (type == "KM") {
          if (is.null(pr_prob_rho)) {
            rho_pr_prob <- 0.5
          } else {
            # check pr prob length
            if (length(pr_prob_rho) != 1 |
                !is.list(pr_prob_rho)) {
              stop("two values for the random effects sd required in a list.")
            } else{
              # check names
              if (isFALSE(all(c("pr_prob") %in% names(pr_prob_rho)))) {
                stop("incorrect names in pr_prob_rho")
              } else {
                # user defined prior prob
                rho_pr_prob <- pr_prob_rho$pr_prob
              }
            }
          } # end of is.null

          # default rho
          if (is.null(z_slab_sd)) {
            z12_slab_sd <- 0.5
          } else {
            # check tau_df length
            if (!is.list(z_slab_sd)) {
              stop("z sd required in a list.")
            } else{
              # check names
              if (isFALSE(all(c("z12_slab_sd") %in% names(z_slab_sd)))) {
                stop("incorrect names in tau_sd")
              } else {
                # user defined prior df
                z12_slab_sd <- z_slab_sd$z12_slab_sd
              }
            }
          }

          model <- paste(
            likelihood,
            mu_structure,
            var_structure,
            re_structure,
            bhat,
            beta,
            re_sigma,
            rho_bern,
            cov_mat
          )

          model_initialize <- jags.model(
            textConnection(model),
            list(
              y = y_scale,
              D = total_REs,
              total_cors = total_cors,
              J = length(unique(ID)),
              n = length(y),
              ID = ID,
              beta0_mean = beta0_mean,
              beta0_sd = beta0_sd,
              eta0_mean = eta0_mean,
              eta0_sd = eta0_sd,
              tau_11_sd = tau11_sd,
              tau_11_df = tau11_df,
              tau_22_sd = tau22_sd,
              tau_22_df = tau22_df,
              z_sd = z12_slab_sd,
              pr_prob_rho = rho_pr_prob
            ),
            n.chains = chains,
            n.adapt = warmup
          )

          samps <- coda.samples(
            model_initialize,
            variable.names = c("beta0", "eta0",
                               "pi", "tau", "rho", "u"),
            n.iter = iter
          )

          res <-  list(
            model = model,
            y = y,
            ID = ID,
            samps = samps,
            model = model,
            model_initialize = model_initialize
          )





        } # end of KM

      } # end of

    } # end mixture


  } # end re = 1

  class(res) <- "muvar"
  return(res)

}
