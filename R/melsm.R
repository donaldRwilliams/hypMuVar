#' @title Mixed-Effects Location Scale Model
#' @name melsm.default
#' @param fixed_location a \code{formula} object for the fixed-effects part of the location sub-model, with the response (on the left) and
#'                       the terms (on the right) sepearted by \code{~}.
#' @param random_location a \code{formula} object for the random effects part of the location sub-model, with the random effect (on the left)
#'                        and the cluster (on the right) seperated by \code{~}.
#' @param fixed_scale a \code{formula} object for the fixed-effects part of the scale sub-model, with the response (on the left) and
#'                       the terms (on the right) sepearted by \code{~}.
#' @param random_scale a \code{formula} object for the random effects part of the location sub-model, with the random effect (on the left)
#'                      and the cluster (on the right) seperated by \code{~}.
#' @param prior a list specifying non-default prior distribution. Set to \code{NULL} for the defaults.
#' @param mixture type of mixture prior distribution for the random effects correlations. Options are \code{mixture  = "KM"} for Kuo and
#'                Mallick (1998) or \code{mixture = SSVS} for stochastic search variable section (X). Default is set to \code{NULL}, wherein
#'                no covariance selection is performed. See notes for futher details.
#' @param k number of mixture components. options are \code{k = 2} or \code{k = 3} for \code{mixture = "SSVS"}.
#'          When \code{mixture = "KM"}. only two components are allowed.
#'
#' @param rho_test \code{rho_test = "all"} tests all of the random effects correlations, whereas \code{rho_test = "muvar"} tests
#'                 the random effects correlations that capture mean--variance relations across the location and scale sub-models.
#' @param adapt adaptive phase. see X.
#' @param chains number of chains.
#' @param iter number of posterior samples for each chain
#' @param thin thinning interval (e.g., with \code{thin = 5} every 5 iteraions are saved)
#' @param data a data frame  containing the variables named in \code{fixed_location}, \code{random_location},
#'             \code{fixed_scale}, and \code{random_scale}
#' @param ... currently ignored
#'
#' @return object of class \code{melsm}
#' @importFrom stats as.formula model.matrix na.omit quantile sd update
#' @import rjags
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
#'}
melsm.default <- function(fixed_location,
                  random_location,
                  fixed_scale,
                  random_scale,
                  prior = NULL,
                  mixture = NULL,
                  k = 2,
                  rho_test = "muvar",
                  adapt = 1000,
                  chains = 4,
                  iter = 5000,
                  thin = 1,
                  data,...) {

  # give null mixture name
  if(is.null(mixture)){
    mixture <- "none"
  }

  # give null prior name
  if(is.null(prior)){
    prior <- "defaults"
  }

  # check k
  if(!any(k %in% c(2, 3))){
    stop("number of components must be 2 or 3")
  }

  if (mixture == "KM") {
    k <- 2
  }

  if(!all(all.vars(fixed_location) %in% colnames(data))){
    stop("fixed_location variables not found in data")
  }

  if(!all(all.vars(fixed_scale)[-1] %in% colnames(data))){
    stop("fixed_scale variables not found in data")
  }

  # split RE formulas
  split_re_loc <- all.vars(random_location)
  split_re_scl <- all.vars(random_scale)

  # outcome
  y <- data[, all.vars(fixed_location)[1]]

  # fixed effects location model matrix
  fe_loc_mm <- model.matrix(fixed_location, data = data)

  # fixed effects scale model matrix
  fe_scl_mm <- model.matrix(as.formula(gsub("sigma.*","",fixed_scale)[-2]), data = data)

  # check for model type
  if (length(split_re_loc) == 1) {
    re_loc_type <- "random_intercept"
  } else {
    re_loc_type <- "random_slope"
  }

  if (length(split_re_scl) == 1) {
    re_scl_type <- "random_intercept"
  } else {
    re_scl_type <- "random_slope"
  }

  # number of location FEs
  betas_loc <- ncol(fe_loc_mm)

  # number of scale FEs
  etas_scl <- ncol(fe_scl_mm)

  # prior names
  prior_names <- names(prior)

  # sd y
  sd_y <- sd(y)

  # mean y
  mean_y <- mean(y)

  # start: random location and scale intercepts
  if (re_loc_type == "random_intercept" &
      re_scl_type == "random_intercept") {

    # tau priors (any user defined)
    if (any(prior_names %in% "tau")) {

      # tau priors
      tau_prior <- prior$tau

      # tau 0 (location intercepts)
      if (any(names(tau_prior) == "tau_0")) {

        # user defined
        tau_0 <- paste("Tau[1,1] ~ ",  tau_prior$tau_0)

        } else {

          # defaults
          tau_0 <- paste("Tau[1,1] ~ dt(0, pow(", sd_y, ",-2), 10)T(0,)")

        }

      # tau 1 (scale intercepts)
      if(any(names(tau_prior) == "tau_1")){

        # user defined
        tau_1 <- paste("Tau[2,2] ~ ",  tau_prior$tau_1)

        } else {

          # defaults
          tau_1 <- paste("Tau[2,2] ~ dt(0, pow(", sd_y, ",-2), 10)T(0,)")
        }

      # tau prior
      tau_prior <-  paste("", "# RE sd priors\n",
                          tau_0, "\n",
                          tau_1, "\n",
                          "Tau[1,2] <- 0\n",
                          "Tau[2,1] <- 0\n")

      # all defaults
      } else {

        tau_prior <- paste("", "# RE sd priors\n",
                         "Tau[1,1] ~ dt(0, pow(", sd_y, ",-2), 10)T(0,)\n",
                         "Tau[2,2] ~ dt(0, pow(", sd_y, ",-2), 10)T(0,)\n",
                         "Tau[1,2] <- 0\n",
                         "Tau[2,1] <- 0\n")

      } # end of tau prior

    # beta priors (any user defined)
    if (any(prior_names %in% "beta")) {

      beta_prior <- prior$beta

      user_prior <- unlist(lapply(1:length(beta_prior), function(x) numextract(names(beta_prior)[x]))) + 1

      # check prior exists
      if(any(user_prior > betas_loc)){

        stop("beta not found")

        }

      beta_prior <- paste("beta[", user_prior, "]~", beta_prior, collapse = "\n", sep = "")

      suppressWarnings(default_beta <- which(user_prior != 1:betas_loc))

      if(length(default_beta) > 0){

      beta_prior <- paste("", "# beta priors\n", beta_prior, "\n",
                          paste("beta[", default_beta, "]~ dnorm(0, 0.001)",
                                collapse = "\n",
                                sep = ""))
       }

      # default priors
      } else {

        beta_prior <-paste("beta[", 1:betas_loc,
                         "]~ dnorm(",c(mean_y, rep(0,betas_loc-1)) , ", 0.001)",
                         sep = "",
                         collapse = "\n")
    } # end beta prior

    # eta priors (any user defined)
    if (any(prior_names %in% "eta")) {

      eta_prior <- prior$eta

      user_prior <- unlist(lapply(1:length(eta_prior), function(x) numextract(names(eta_prior)[x]))) + 1

      # check prior exists
      if (any(user_prior > etas_scl)) {

        stop("eta not found")

        }

      eta_prior <- paste("eta[", user_prior, "]~", eta_prior, collapse = "\n", sep = "")

      suppressWarnings(default_eta <- which(user_prior != 1:etas_scl))

      if (length(default_eta) > 0) {

        eta_prior <- paste("", "# eta priors\n", eta_prior, "\n",
                          paste("eta[", default_eta, "]~ dnorm(0, 0.001)",
                               collapse = "\n",
                               sep = ""))
        }

    # default priors
    } else {

      eta_prior <-paste("eta[", 1:etas_scl,
                        "]~ dnorm(0, 0.001)",
                        sep = "", collapse = " \n")
    } # end eta prior

    # start: rho priors (user defined)
    if(any(prior_names %in% "rho")){

      rho_prior <- prior$rho


      if(!any(c('slab') %in% names(rho_prior))) {

        stop("spike and slab standard deviation must be specified")

      }

      # start: SSVS
      if(mixture == "SSVS"){

        # check spike and slab
        if(!all(c("spike", 'slab') %in% names(rho_prior))) {
          stop("spike and slab standard deviation must be specified")
        }

        # check components
        if(k == 0 | is.null(k)){
          stop("k (number of mixture components must be specified)")
        }

        # start: two components (k = 2)
        if (k == 2) {

          rho_prior <- paste("for(k in 1:cors){\n",
                             "k_rho[k] ~ dcat(pr_k[])\n",
                             "z[k] ~ dnorm(0, rho_sd[k_rho[k]]))\n",
                             "rho[1] <- tanh(z[k])}\n",
                             paste("rho_sd[1] <- pow(", rho_prior$spike , ", -2)\n"),
                             paste("rho_sd[2] <- pow(", rho_prior$slab , ", -2)\n"),
                             "pr_k[1] <- 1/2\n",
                             "pr_k[2] <- 1/2\n")

          } # end: two components (k = 2)

        # start: three components (k = 3)
        if (k == 3) {

          rho_prior <-  paste("for(k in 1:cors){\n",
                              "k_rho[k] ~ dcat(pr_k[])\n",
                              "z[k] ~ dnorm(0,  rho_sd[k_rho[k]])T(T1[k_rho[k]], T2[k_rho[k]])\n",
                              "rho[k] <- tanh(z[k]) \n}\n\n",
                              paste("rho_sd[1] <- pow(", rho_prior$spike , "-2)\n"),
                              paste("rho_sd[2] <- pow(", rho_prior$slab , "-2)\n"),
                              paste("rho_sd[3] <- pow(", rho_prior$slab , "-2)\n"),
                              "p_rho[1] <- 1/3\n",
                              "p_rho[2] <- 1/3\n",
                              "p_rho[3] <- 1/3\n\n",
                              "T1[1] <- -10\n",
                              "T1[2] <- 0\n",
                              "T1[3] <- -10\n",
                              "T2[1] <- 10\n",
                              "T2[2] <- 10\n",
                              "T2[3] <- 0\n")

          } # end: three components (k = 3)

        # start: kuo & mallick (1998)
        } else if (mixture == "KM") {

          rho_prior <- paste("for(k in 1:cors){\n",
                             "k_rho[k] ~ dbern(0.5)\n",
                             paste("z[k] ~ dnorm(0, pow(", rho_prior$slab, ", -2))\n"),
                             "rho[k] <- tanh(z[k]) * k_rho[k] \n}\n")

          } # end: kuo & mallick (1998

      # start: user sd no mixture
      else {

        rho_prior <- paste("for(k in 1:cors){\n",
                           paste("z[k] ~ dnorm(0, pow(", rho_prior$slab, ", -2))\n"),
                           "rho[k] <- tanh(z[k]) }\n")

        } # end: user sd no mixture

      # start: rho priors (defaults)
      } else {

        # start: SSVS
        if (mixture == "SSVS") {

          # start: two components (k = 2)
          if (k == 2) {

            rho_prior <- paste("for(k in 1:cors){\n",
                               "k_rho[k] ~ dcat(pr_k[])\n",
                               "z[k] ~ dnorm(0, rho_sd[k_rho[k]])\n",
                               "rho[k] <- tanh(z[k])}\n",
                               "rho_sd[1] <- pow(0.01, -2)\n",
                               "rho_sd[2] <- pow(0.50, -2)\n",
                               "pr_k[1] <- 1/2\n",
                               "pr_k[2] <- 1/2\n")

            } # end: two components (k = 2)

          # start: three components (k = 3)
          if (k == 3) {

            rho_prior <-  paste("for(k in 1:cors){\n",
                                "k_rho[k] ~ dcat(pr_k[])\n",
                                "z[k] ~ dnorm(0,  rho_sd[k_rho[k]])T(T1[k_rho[k]], T2[k_rho[k]])\n",
                                "rho[k] <- tanh(z[k]) \n}\n\n",
                                "rho_sd[1] <- pow(0.01, -2)\n",
                                "rho_sd[2] <- pow(0.50, -2)\n",
                                "rho_sd[3] <- pow(0.50, -2)\n\n",
                                "pr_k[1] <- 1/3\n",
                                "pr_k[2] <- 1/3\n",
                                "pr_k[3] <- 1/3\n\n",
                                "T1[1] <- -10\n",
                                "T1[2] <- 0\n",
                                "T1[3] <- -10\n",
                                "T2[1] <- 10\n",
                                "T2[2] <- 10\n",
                                "T2[3] <- 0\n")

            } # end: three component (k = 3)

          # start: kuo & mallick (1998)
          } else if (mixture == "KM") {

            rho_prior <- paste("for(k in 1:cors){\n",
                               "k_rho[k] ~ dbern(0.5)\n",
                               "z[k] ~ dnorm(0, pow(0.5, -2))\n",
                               "rho[k] <- tanh(z[k]) * k_rho[k]}\n")

            # start: no mixture (default)
            } else {

              rho_prior <- paste("for(k in 1:cors){\n",
                           "z[k] ~ dnorm(0, pow(0.5, -2))\n",
                           "rho[k] <- tanh(z[k])\n}\n")
              } # end: no mixture (default)
        } # end: rho priors (defaults)

    # priors
    priors <-  paste(rho_prior, "\n",
                     beta_prior, "\n\n",
                     eta_prior, "\n\n",
                     tau_prior)

    # jags model
    mod <- paste(re_int_both,
                 priors, "}")

    # scale model matrix for the random int
    re_scl_mm <- model.matrix( ~ 1, data = data)

    # location model matrix for the random int
    re_loc_mm <- model.matrix( ~ 1, data = data)

    # model matrices
    model_matrices <- list(fe_loc_mm = fe_loc_mm,
                           fe_scl_mm = fe_scl_mm,
                           re_loc_mm = re_loc_mm,
                           re_scl_mm = re_scl_mm)

    # intialize model
    suppressWarnings(m_initialize <- jags.model(file = textConnection(mod),
                                               n.chains = chains,
                                               data = list(X_loc = fe_loc_mm,
                                                           y = y, X_scl = fe_scl_mm,
                                                           J = length(unique(data[,split_re_loc[1]])),
                                                           N = nrow(data), ID = data[,split_re_loc[1]],
                                                           cors = 1)))
    # update model with adaption phase
    update(m_initialize, n.iter = adapt)

    # estimate model
    suppressWarnings(m_samps <- coda.samples(model = m_initialize,thin = thin,
                                             variable.names = c("rho", "beta",
                                                                "eta", "u",
                                                                "Tau", "k_rho"),
                                             n.iter = iter))

    ret <- list(melsm_fit = m_samps,
                melsm_initialize = m_initialize,
                jags_model = mod,
                model_matrices = model_matrices,
                dat = data,
                iter = iter,
                cluster = split_re_loc[1],
                chains = chains,
                thin = thin,
                mixture = mixture,
                rho_test  = rho_test,
                k = k, priors = priors,
                call = match.call())

    } # end: both random loc and scl intercepts

  if (re_loc_type == "random_slope" &
      re_scl_type == "random_intercept") {

    # cluster
    clusters <- unique(data[,split_re_loc[2]])

    # check the random slope is valid (i.e., level 1)
    test_level_1 <- subset(data, data[,split_re_loc[2]]   == clusters[1] )
    test_level_1 <- length(unique(test_level_1[,split_re_loc[1]]))

    # random slope model matrix
    re_loc_mm <- model.matrix(as.formula(paste("~", split_re_loc[1])), data = data)

    # scale model matrix for the random int
    re_scl_mm <- model.matrix(~ 1, data = data)

    # stop if invalid random slope
    if(test_level_1 == 1){
      stop("random slope cannot be a level 2 predictor")
    }

    # check random slope is also a fixed effect
    if(!any(colnames(fe_loc_mm) == colnames(re_loc_mm)[2])){
      stop("random slope must included as a fixed effect")
    }

    # which FE predictor is the random slope
    ran_slp <- which(colnames(fe_loc_mm) == colnames(re_loc_mm)[2])

    # random intercept
    re_loc_int_form <- paste("beta[1] + u[ID[i], 1]")

    # random slope
    re_loc_slp_form <- paste("(beta[", ran_slp,
                             "]", " + u[ID[i],2])",
                             " * X_loc_", ran_slp, "[i]",
                             sep = ""
                             )

    # sequence for number of predictors
    temp <- 1:betas_loc

    # additional level two predictors
    if (length(temp) > 2) {

      # remove intercept and random slope, leaving remaining predictors
      level_2_loc <- temp[-c(1, ran_slp)]

      # level 2 location
      level_2_loc_form <- paste('beta[',
                                level_2_loc, "]",
                                " * X_loc_", level_2_loc,
                                "[i]", sep = "",
                                collapse = " + ")

      # likelihood
      likelihood <- paste("for(i in 1:N){\n",
                          "y[i] ~ dnorm(yhat[i],  1/exp(shat[i])^2) \n",
                          "yhat[i]<-", re_loc_int_form,  " + ",
                          re_loc_slp_form, " + ", level_2_loc_form,
                          "\n","shat[i] <- inprod(X_scl[i,], eta)  + u[ID[i], 3] \n",
                          "}\n", re_slp_loc)

      # only random intercept and slope
      } else {

        # likelihood
        likelihood <- paste("for(i in 1:N){\n",
                          "y[i] ~ dnorm(yhat[i],  1/exp(shat[i])^2) \n",
                          "yhat[i]<-", re_loc_int_form, " + ",
                          re_loc_slp_form,
                          "\n","shat[i] <- inprod(X_scl[i,], eta)  + u[ID[i], 3] \n",
                          "}\n", re_slp_loc)
        }

    # tau priors
    if (!unlist(any(prior_names %in% "tau"))) {

      # no user defined
      suppressWarnings(prior$tau <- NULL)

      }

    # tau prior
    tau_prior <- prior$tau
    tau_list <- list()

    # loop tau priors
    for(i in 1:3){

      # user defined
      if(any(names(tau_prior) ==  paste("tau_", i - 1,  sep = "") )){

        tau_list[[i]] <- paste("Tau[",i , ",", i ,  "]"," ~ ",
                               tau_prior[[paste("tau_", i - 1,  sep = "")]],
                               sep = "")
        # default
        } else {

          tau_list[[i]] <- paste("Tau[",i , ",", i ,  "]"," ~ ",
                                 "dt(0, pow(", sd_y, ",-2), 10)T(0,)",
                                 sep = "")
        }
      }

    # unlist
    tau_prior_temp <- paste(tau_list, collapse = "\n")

    # tau prior
    tau_prior <-  paste(tau_prior_temp, "\n",
                        "Tau[1,2] <- 0\n",
                        "Tau[2,1] <- 0\n",
                        "Tau[3,1] <- 0\n",
                        "Tau[1,3] <- 0\n",
                        "Tau[2,3] <- 0\n",
                        "Tau[3,2] <- 0\n",
                        collapse = "\n")

    # beta priors
    if (!any(prior_names %in% "beta")) {

      # no user defined
      prior$beta <- 1

    }

    # beta prior
    beta_prior <- prior$beta
    beta_list <- list()

    # user priors
    user_prior <- as.numeric(unlist(lapply(1:length(beta_prior),
                                           function(x) numextract(names(beta_prior)[x])))) + 1

    # check beta exists
    if (any(user_prior > betas_loc)) {

      stop("beta not found")

      }

    # loop betas
    for (i in 1:betas_loc) {

      # user defined
      if( any(names(beta_prior) ==  paste("beta_", i - 1,  sep = ""))) {

        beta_list[[i]] <- paste("beta[",i,"]"," ~ ",
                                 beta_prior[[paste("beta_", i - 1,  sep = "")]],
                                 "\n", sep = "")
        # defaults
        } else {

          # intercept
          if (i == 1) {

            beta_list[[i]] <- paste("beta[1] ~ ",
                                    "dnorm(", mean_y, ",0.001) \n",
                                    sep = "")

            # not intercept
            } else {

              beta_list[[i]] <- paste("beta[", i, "] ~ ",
                                      "dnorm(", 0, ",0.001)\n",
                                      sep = "")

            }
        }
      }

    # beta prior
    beta_prior <- paste(beta_list, collapse = "")

    # eta priors
    if (!any(prior_names %in% "eta")) {

      # no user defined
      prior$eta <- 1

      }

    # eta priors
    eta_prior <- prior$eta
    eta_list <- list()

    # user defined priors
    user_prior <- as.numeric(unlist(lapply(1:length(eta_prior),
                                           function(x) numextract(names(eta_prior)[x])))) + 1

    # check etas exist
    if(any(user_prior > etas_scl)){

      stop("eta not found")

      }

    # loop eta
    for(i in 1:etas_scl){

      # user defined
      if (any(names(eta_prior) ==  paste("eta_", i -1,  sep = ""))) {

        eta_list[[i]] <- paste("eta[",i,"]"," ~ ",
                               eta_prior[[paste("eta_", i - 1,  sep = "")]],
                               "\n", sep = "")

        # defaults
        } else {

          eta_list[[i]] <- paste("eta[", i, "] ~ ",
                                 "dnorm(", 0, ",0.001)\n",
                                 sep = "")

        }
      }

    # eta prior
    eta_prior <- paste(eta_list, collapse = "")


    # which RE correlations are tested
    # mean--variance
    if (rho_test == "muvar") {

      cors <- c(2,3)
      not_tested_cors <- which(1:3 %in% cors == FALSE)

      # all
      } else if (rho_test == "all"){

        cors <- 1:3
        not_tested_cors <- 0

        # error
        } else {

          stop("rho_test must be 'muvar' or 'all'")

        }

    # rho prior
    if (any(prior_names %in% "rho")) {

      # user defined
      rho_priors <- prior$rho

      } else {

        # default spike and slab
        rho_priors <- list(spike = 0.01, slab = 0.5)

        }

    # mixture: SSVS
    if (mixture == "SSVS") {

      # check spike and slab defined
      if (!all(c("spike", 'slab') %in% names(rho_priors))) {

        stop("spike and slab standard deviation must be specified")

        }

      if (k == 0 | is.null(k)) {

        stop("k (number of mixture components) must be specified)")

        }

      # two components
      if (k == 2) {

        rho_prior <- paste("for(k in cors){\n",
                           "k_rho[k] ~ dcat(pr_k[])\n",
                           "z[k] ~ dnorm(0, rho_sd[k_rho[k]])\n",
                           "rho[k] <- tanh(z[k])\n}\n",
                           paste("rho_sd[1] <- pow(", rho_priors$spike , ", -2)\n"),
                           paste("rho_sd[2] <- pow(", rho_priors$slab , ", -2)\n"),
                           "pr_k[1] <- 1/2\n",
                            "pr_k[2] <- 1/2\n")
        # not tested
        if(any(not_tested_cors != 0)){

          rho_prior <- paste("# rho(s) not tested \n",
                             "for(k in not_tested_cors){\n",
                             "z[k] ~ dnorm(0, pow(",   rho_priors$slab, ", -2)) \n",
                             "rho[k] <- tanh(z[k])\n}\n",
                             "# rho test \n",
                             rho_prior,
                             collapse = "\n\n")
          }
        }# end two component

      # three component
      if (k == 3) {

        rho_prior <-  paste("for(k in cors){\n",
                            "k_rho[k] ~ dcat(pr_k[])\n",
                            "z[k] ~ dnorm(0,  rho_sd[k_rho[k]])T(T1[k_rho[k]], T2[k_rho[k]])\n",
                            "rho[k] <- tanh(z[k]) \n}\n\n",
                            paste("rho_sd[1] <- pow(", rho_priors$spike, ",-2)\n"),
                            paste("rho_sd[2] <- pow(", rho_priors$slab, ",-2)\n"),
                            paste("rho_sd[3] <- pow(", rho_priors$slab, ",-2)\n"),
                            "pr_k[1] <- 1/3\n",
                            "pr_k[2] <- 1/3\n",
                            "pr_k[3] <- 1/3\n\n",
                            "T1[1] <- -10\n",
                            "T1[2] <- 0\n",
                            "T1[3] <- -10\n",
                            "T2[1] <- 10\n",
                            "T2[2] <- 10\n",
                            "T2[3] <- 0\n")

        # not tested
        if(any(not_tested_cors != 0)){

          rho_prior <- paste("# rho(s) not tested \n",
                             "for(k in not_tested_cors){\n",
                             "z[k] ~ dnorm(0, pow(", rho_priors$slab, ", -2)) \n",
                             "rho[k] <- tanh(z[k])\n}\n",
                             "# rho test \n",
                             rho_prior, collapse = "\n\n")
          }
      } # end three component

      # mixture KM
    } else if (mixture == "KM") {

      rho_prior <- paste("for(k in cors){\n",
                         "k_rho[k] ~ dbern(0.5)\n",
                         paste("z[k] ~ dnorm(0, pow(", rho_priors$slab, ", -2))\n"),
                        "rho[k] <- tanh(z[k]) * k_rho[k] \n}\n")

      # not tested
      if (any(not_tested_cors != 0)) {

        rho_prior <- paste("# rho(s) not tested \n",
                           "for(k in not_tested_cors){\n",
                           "z[k] ~ dnorm(0, pow(",   rho_priors$slab, ", -2)) \n",
                           "rho[k] <- tanh(z[k])\n}\n",
                           "# rho test \n",
                           rho_prior,
                           collapse = "\n\n")
        }

      # no mixture
      } else {

        rho_prior <- paste("# no cors tested\n",
                           "for(k in 1:3){\n",
                           paste("z[k] ~ dnorm(0, pow(", rho_priors$slab, ", -2))\n"),
                           "   rho[k] <- tanh(z[k])\n}\n", sep = "")
        }

    # model matrices
    model_matrices <- list(fe_loc_mm = fe_loc_mm,
                           fe_scl_mm = fe_scl_mm,
                           re_loc_mm = re_loc_mm,
                           re_scl_mm = re_scl_mm)

    # data for jags
    dat_list <- lapply(2:betas_loc, function(x) fe_loc_mm[,x])
    names(dat_list) <- paste("X_loc_", 2:betas_loc, sep = "")
    dat_list$y <- y
    dat_list$J <- length(unique(data[,split_re_loc[2]]))
    dat_list$N <- nrow(data)
    dat_list$ID <- data[,split_re_loc[2]]
    dat_list$X_scl <- fe_scl_mm
    dat_list$cors <- cors

    if(length(not_tested_cors) != 0){

      dat_list$not_tested_cors <- not_tested_cors

      }

    # jags model
    mod <- paste("model{\n", likelihood, "\n",
                 "# beta priors \n", beta_prior, "\n",
                 "# eta priors  \n",  eta_prior, "\n",
                 "# tau priors \n",
                 tau_prior, "\n\n", "# rho priors \n",
                 rho_prior, "\n}")

    # priors
    priors <-  paste(rho_prior, "\n",
                     beta_prior, "\n\n",
                     eta_prior, "\n\n",
                     tau_prior)

    # initialize model
    suppressWarnings(m_initialize <- jags.model(textConnection(mod),
                                                data = dat_list,
                                                n.chains = chains))
    # adaptive period
    update(m_initialize, n.iter = adapt)

    # fit model
    suppressWarnings(m_samps <- coda.samples(m_initialize,
                                             variable.names = c("beta", "eta", "u",
                                                                "rho", "Tau", "k_rho"),
                                             n.iter = iter))
    # returned object
    ret <- list(melsm_fit = m_samps,
                melsm_initialize = m_initialize,
                jags_model = mod,
                model_matrices = model_matrices,
                dat = data,
                cluster = split_re_loc[2],
                iter = iter,
                chains = chains,
                thin = thin,
                mixture = mixture,
                rho_test  = rho_test,
                k = k, priors = priors)



  } #end: random loc slopes and random intercepts scl


  if (re_loc_type == "random_intercept" &
      re_scl_type == "random_slope") {

    # cluster
    clusters <- unique(data[,split_re_scl[2]])

    # check the random slope is valid
    test_level_1 <- subset(data, data[,split_re_scl[2]]   == clusters[1])
    test_level_1 <- length(unique(test_level_1[,split_re_scl[1]]))

    if(test_level_1 == 1){
      stop("random slope cannot be a level 2 predictor")
    }

    # model matrix for the random slope
    re_scl_mm <- model.matrix(as.formula( paste("~", split_re_scl[1])), data = data)

    # random slope model matrix
    re_loc_mm <- model.matrix( ~ 1, data = data)

    # check random slope is also a fixed effect
    if(!any(colnames(fe_scl_mm) == colnames(re_scl_mm)[2])){
      stop("random slope must included as a fixed effect")
    }

    # which FE predictor is the random slope
    ran_slp <- which(colnames(fe_scl_mm) == colnames(re_scl_mm)[2])

    # random intercept
    re_scl_int_form <- paste("eta[1] + u[ID[i], 2]")

    # random slope
    re_scl_slp_form <- paste("(eta[", ran_slp,
                             "]", " + u[ID[i],3])",
                             " * X_scl_", ran_slp, "[i]",
                             sep = "")

    # sequence for number of predictors
    temp <- 1:etas_scl

    # level two predictors
    if(length(temp) > 2){

      #remove intercept and random slope, leaving remaining predictors
      level_2_scl <- temp[-c(1, ran_slp)]


      level_2_scl_form <- paste('eta[', level_2_scl, "]",
                                " * X_scl_", level_2_scl, "[i]",
                                sep = "", collapse = " + ")

      likelihood <- paste("for(i in 1:N){\n",
                          "y[i] ~ dnorm(yhat[i],  1/exp(shat[i])^2) \n",
                          "shat[i]<-", re_slc_int_form,  " + ",
                          re_scl_slp_form, " + ", level_2_scl_form,
                          "\n","yhat[i] <- inprod(X_loc[i,], beta)  + u[ID[i], 1] \n", "}\n",
                          re_slp_scl)

    } else{

      likelihood <- paste("for(i in 1:N){\n",
                          "y[i] ~ dnorm(yhat[i],  1/exp(shat[i])^2) \n",
                          "yhat[i] <- inprod(X_loc[i,], beta)  + u[ID[i], 1] \n",
                          "shat[i]<-", re_scl_int_form, " + ",
                          re_scl_slp_form, "\n}\n",
                          re_slp_scl)
    }

    # tau priors
    if (!unlist(any(prior_names %in% "tau"))) {

      # no user defined
      suppressWarnings(prior$tau <- NULL)

      }

    # tau priors
    tau_prior <- prior$tau
    tau_list <- list()

    # tau prior loop
    for(i in 1:3){

      # user defined
      if (any(names(tau_prior) ==  paste("tau_", i - 1,  sep = ""))) {

        tau_list[[i]] <- paste("Tau[", i , ",", i ,  "]"," ~ ",
                               tau_prior[[paste("tau_", i - 1,  sep = "")]],
                               sep = "")

        } else {

          # all defaults
          tau_list[[i]] <- paste("Tau[",i , ",", i ,  "]"," ~ ",
                                 "dt(0, pow(", sd_y, ",-2), 10)T(0,)",
                                 sep = "")
        }
      }

    # unlist tau priors
    tau_prior_temp <- paste(tau_list,
                            collapse = "\n",
                            sep = "")

    # tau priors
    tau_prior <-  paste(tau_prior_temp, "\n",
                        "Tau[1,2] <- 0\n",
                        "Tau[2,1] <- 0\n",
                        "Tau[3,1] <- 0\n",
                        "Tau[1,3] <- 0\n",
                        "Tau[2,3] <- 0\n",
                        "Tau[3,2] <- 0\n",
                        collapse = "\n",
                        sep = "")

    # beta priors
    if (!any(prior_names %in% "beta")) {

      # no user defined
      prior$beta <- 1

    }

    # beta priors
    beta_prior <- prior$beta
    beta_list <- list()

    # user priors numeric extract
    user_prior <- as.numeric(unlist(lapply(1:length(beta_prior),
                                           function(x) numextract(names(beta_prior)[x])))) + 1
    # check beta exists
    if(any(user_prior > betas_loc)){

      stop("beta not found")

      }

    # beta loop
    for(i in 1:betas_loc){

      # user defined
      if(any(names(beta_prior) ==  paste("beta_", i - 1,  sep = "") )){

        beta_list[[i]] <- paste("beta[", i ,"]"," ~ ",
                                beta_prior[[paste("beta_", i - 1,  sep = "")]],
                                "\n", sep = "")

        # defaults
        } else {

          # intercept
          if(i == 1){

            beta_list[[i]] <- paste("beta[1] ~ ",
                                  "dnorm(", mean_y, ",0.001) \n",
                                  sep = "")

            # non-intercepts
            } else {

              beta_list[[i]] <- paste("beta[", i , "] ~ ",
                                    "dnorm(", 0, ",0.001)\n",
                                    sep = "")
            }
        }
      }

    # beta prior
    beta_prior <- paste(beta_list, collapse = "")

    # eta priors
    if (!any(prior_names %in% "eta")) {

      # no user defined
      prior$eta <- 1

      }

    # eta priors
    eta_prior <- prior$eta
    eta_list <- list()

    # eta numeric extract
    user_prior <- as.numeric(unlist(lapply(1:length(eta_prior),
                                           function(x) numextract(names(eta_prior)[x])))) + 1

    # check eta exists
    if (any(user_prior > etas_scl)) {

      stop("eta not found")

      }

    # eta loop
    for(i in 1:etas_scl){

      # user defined
      if(any(names(eta_prior) ==  paste("eta_", i - 1,  sep = "") )){

        eta_list[[i]] <- paste("eta[",i,"]"," ~ ",
                               eta_prior[[paste("eta_", i - 1,  sep = "")]],
                               "\n", sep = "")
        # defaults
        } else {

          eta_list[[i]] <- paste("eta[", i, "] ~ ",
                                 "dnorm(", 0, ",0.001)\n",
                                 sep = "")
        }
      }

    # eta prior
    eta_prior <- paste(eta_list, collapse = "")

    # rho prior
    rho_test <- unlist(rho_test)
    cors <- NA

    # test mean--varinace relations
    if (rho_test == "muvar") {

      cors <- c(1,2)
      not_tested_cors <- which(1:3 %in% cors == FALSE)

      }

    # test all correlations
    if (rho_test == "all") {

      cors <- 1:3
      not_tested_cors <- 0

      }

    # start: build rho prior
    if(any(prior_names %in% "rho")){

      # user defined
      rho_priors <- prior$rho

      } else {

        # defaults
        rho_priors <- list(spike = 0.01, slab = 0.5)

        }

    # mixture SSVS
    if(mixture == "SSVS"){

      # check for spike and slab
      if (!all(c("spike", 'slab') %in% names(rho_priors))) {

        stop("spike and slab standard deviation must be specified")

        }

      # check k
      if(k == 0 | is.null(k)){

        stop("k (number of mixture components must be specified)")

        }

      # two components
      if (k == 2) {

        rho_prior <- paste("for(k in cors){\n",
                           "k_rho[k] ~ dcat(pr_k[])\n",
                           "z[k] ~ dnorm(0, rho_sd[k_rho[k]])\n",
                           "rho[k] <- tanh(z[k])\n}\n",
                           paste("rho_sd[1] <- pow(", rho_priors$spike , ", -2)\n"),
                           paste("rho_sd[2] <- pow(", rho_priors$slab , ", -2)\n"),
                           "pr_k[1] <- 1/2\n",
                           "pr_k[2] <- 1/2\n")

        # not tested
        if (any(not_tested_cors != 0)) {

          rho_prior <- paste("# rho(s) not tested \n",
                             "for(k in not_tested_cors){\n",
                             "z[k] ~ dnorm(0, pow(",   rho_priors$slab, ", -2)) \n",
                             "rho[k] <- tanh(z[k])\n}\n",
                             "# rho test \n",
                             rho_prior, collapse = "\n\n")
          }
        }# end: k = 2

      # three components
      if(k == 3){

        rho_prior <-  paste("for(k in cors){\n",
                            "k_rho[k] ~ dcat(pr_k[])\n",
                            "z[k] ~ dnorm(0,  rho_sd[k_rho[k]])T(T1[k_rho[k]], T2[k_rho[k]])\n",
                            "rho[k] <- tanh(z[k]) \n}\n\n",
                            paste("rho_sd[1] <- pow(", rho_priors$spike, ",-2)\n"),
                            paste("rho_sd[2] <- pow(", rho_priors$slab, ",-2)\n"),
                            paste("rho_sd[3] <- pow(", rho_priors$slab, ",-2)\n"),
                            "pr_k[1] <- 1/3\n",
                            "pr_k[2] <- 1/3\n",
                            "pr_k[3] <- 1/3\n\n",
                            "T1[1] <- -10\n",
                            "T1[2] <- 0\n",
                            "T1[3] <- -10\n",
                            "T2[1] <- 10\n",
                            "T2[2] <- 10\n",
                            "T2[3] <- 0\n")

        if (any(not_tested_cors != 0)) {

          rho_prior <- paste("# rho(s) not tested \n",
                             "for(k in not_tested_cors){\n",
                             "z[k] ~ dnorm(0, pow(", rho_priors$slab, ", -2)) \n",
                             "rho[k] <- tanh(z[k])\n}\n",
                             "# rho test \n",
                             rho_prior, collapse = "\n\n")
          }
        } # end: k = 3

      # mixture KM
      } else if (mixture == "KM") {

        rho_prior <- paste("for(k in cors){\n",
                           "k_rho[k] ~ dbern(0.5)\n",
                           paste("z[k] ~ dnorm(0, pow(", rho_priors$slab, ", -2))\n"),
                           "rho[k] <- tanh(z[k]) * k_rho[k] \n}\n")

        # not tested
        if (any(not_tested_cors != 0)) {

          rho_prior <- paste("# rho(s) not tested \n",
                           "for(k in not_tested_cors){\n",
                           "z[k] ~ dnorm(0, pow(",   rho_priors$slab, ", -2)) \n",
                           "rho[k] <- tanh(z[k])\n}\n",
                           "# rho test \n",
                           rho_prior, collapse = "\n\n")
        }

        } else {

          rho_prior <- paste("# no cors tested\n",
                             "for(k in 1:3){\n",
                             paste("   z[k] ~ dnorm(0, pow(", rho_priors$slab, ", -2))\n"),
                             "   rho[k] <- tanh(z[k])\n}\n", sep = "")
          }

    # model matrices
    model_matrices <- list(fe_loc_mm = fe_loc_mm,
                           fe_scl_mm = fe_scl_mm,
                           re_loc_mm = re_loc_mm,
                           re_scl_mm = re_scl_mm)

    # scale predictors
    dat_list <- lapply(2:etas_scl, function(x) fe_scl_mm[,x])

    # scale predictors names
    names(dat_list) <- paste("X_scl_", 2:etas_scl, sep = "")

    # outcome
    dat_list$y <- y
    dat_list$J <- length(unique(data[,split_re_loc[1]]))
    dat_list$N <- nrow(data)
    dat_list$ID <- data[,split_re_loc[1]]
    dat_list$X_loc <- fe_loc_mm
    dat_list$cors <- cors


    if (length(not_tested_cors) > 0) {

      dat_list$not_tested_cors <- not_tested_cors

      }

    mod <- paste("model{\n",
                 likelihood, "\n",
                 "# beta priors\n",
                 beta_prior, "\n",
                 "# eta priors\n",
                 eta_prior, "\n",
                 "# tau priors\n",
                 tau_prior, "\n\n",
                 "# rho priors\n",
                 rho_prior, "\n}")

    # priors
    priors <-  paste(rho_prior, "\n",
                     beta_prior, "\n\n",
                     eta_prior, "\n\n",
                     tau_prior)

    suppressWarnings(m_initialize <- jags.model(textConnection(mod),
                                            data = dat_list,
                                            n.chains = chains))

    update(m_initialize, n.iter = adapt)

    suppressWarnings(m_samps <- coda.samples(m_initialize,thin = thin,
                                             variable.names = c("beta", "eta",
                                                                "u", "rho",
                                                                "Tau", "k_rho"),
                            n.iter = iter))

    ret <- list(melsm_fit = m_samps,
                melsm_initialize = m_initialize,
                jags_model = mod,
                model_matrices = model_matrices,
                dat = data,
                cluster = split_re_loc[1],
                iter = iter,
                chains = chains,
                thin = thin,
                mixture = mixture,
                rho_test  = rho_test,
                k = k, priors = priors)
    } # end: scale random slope and location random intercept
#
  if (re_loc_type == "random_slope" &
      re_scl_type == "random_slope") {


    # cluster
    clusters <- unique(data[,split_re_scl[2]])

    # scale check to ensure the random slope is valid
    test_level_1_scl <- subset(data, data[,split_re_scl[2]]   == clusters[1] )
    test_level_1_scl <- length(unique(test_level_1_scl[,split_re_scl[1]]))

    # scale check to ensure the random slope is valid
    test_level_1_loc <- subset(data, data[,split_re_loc[2]]   == clusters[1] )
    test_level_1_loc <- length(unique(test_level_1_loc[,split_re_loc[1]]))

    # scale model matrix for the random slope
    re_scl_mm <- model.matrix(as.formula( paste("~", split_re_scl[1])), data = data)

    # location model matrix for the random slope
    re_loc_mm <- model.matrix(as.formula( paste("~", split_re_loc[1])), data = data)


    if(test_level_1_loc == 1){
      stop("location random slope cannot be a level 2 predictor")
    }

    if(test_level_1_scl == 1){
      stop("scale random slope cannot be a level 2 predictor")
    }

    # check random slope is also a fixed effect
    if(!any(colnames(fe_scl_mm) == colnames(re_scl_mm)[2])){
      stop("scale random slope must included as a fixed effect")
    }

    if(!any(colnames(fe_loc_mm) ==  colnames(re_loc_mm)[2])){
      stop("location random slope must included as a fixed effect")
    }

    model_matrices <- list(fe_loc_mm = fe_loc_mm,
                           fe_scl_mm = fe_scl_mm,
                           re_loc_mm = re_loc_mm,
                           re_scl_mm = re_scl_mm)


    # scale: which FE predictor is the random slope
    ran_slp_scl <- which(colnames(fe_scl_mm) == colnames(re_scl_mm)[2])

    # location: which FE predictor is the random slope
    ran_slp_loc <- which(colnames(fe_loc_mm) ==  colnames(re_loc_mm)[2])


    # scale: random intercepts
    re_scl_int_form <- paste("eta[1] + u[ID[i], 3]")

    # location: random intercepts
    re_loc_int_form <- paste("beta[1] + u[ID[i], 1]")

    # scale: random slope form
    re_scl_slp_form <- paste("(eta[", ran_slp_scl,
                             "]", " + u[ID[i], 4])",
                             " * X_scl_", ran_slp_scl, "[i]",
                             sep = "")


    # location: random slope form
    re_loc_slp_form <- paste("(beta[", ran_slp_loc,
                             "]", " + u[ID[i], 2])",
                             " * X_loc_", ran_slp_loc, "[i]",
                             sep = "")



    # sequence for number of predictors
    temp_scl <- 1:etas_scl
    temp_loc <- 1:betas_loc

    if(length(temp_scl) > 2){

      #remove intercept and random slope, leaving remaining predictors
      level_2_scl <- temp_scl[-c(1, ran_slp_scl)]


      level_2_scl_form <- paste('eta[', level_2_scl, "]",
                                " * X_scl_", level_2_scl, "[i]",
                                sep = "", collapse = " + ")

      scl_form <- paste("shat[i] <- ", re_scl_int_form,
                        " + ", re_scl_slp_form,
                        " + ", level_2_scl_form,
                        sep = "")

    } else {

      scl_form <- paste("shat[i] <- ", re_scl_int_form,
                        " + ", re_scl_slp_form,
                        sep = "")


    }



    if(length(temp_loc) > 2){

      #remove intercept and random slope, leaving remaining predictors
      level_2_loc <- temp_loc[-c(1, ran_slp_loc)]


      level_2_loc_form <- paste('beta[', level_2_loc, "]",
                                " * X_loc_", level_2_loc, "[i]",
                                sep = "", collapse = " + ")

      loc_form <- paste("yhat[i] <- ", re_loc_int_form,
                        " + ", re_loc_slp_form,
                        " + ", level_2_loc_form,
                        sep = "")

    } else {

      loc_form <- paste("yhat[i] <- ", re_loc_int_form,
                        " + ", re_loc_slp_form,
                        sep = "")

    }


    likelihood <- paste("for(i in 1:N){\n",
                        "y[i] ~ dnorm(yhat[i], 1/exp(shat[i])^2)\n",
                        loc_form, "\n",
                        scl_form, "\n",
                        "}\n\n", Omega_re_slp_both)

    # start: tau priors
    if(!unlist(any(prior_names %in% "tau"))){
      suppressWarnings(prior$tau <- NULL)
    }


    tau_prior <- prior$tau
    tau_list <- list()

    for(i in 1:4){
      if(any(names(tau_prior) ==  paste("tau_", i - 1,  sep = "") )){
        tau_list[[i]] <- paste("Tau[",i , ",", i ,  "]"," ~ ",
                               tau_prior[[paste("tau_", i - 1,  sep = "")]],
                               sep = "")
      } else {
        tau_list[[i]] <- paste("Tau[",i , ",", i ,  "]"," ~ ",
                               "dt(0, pow(", mean_y, ",-2), 10)T(0,)", sep = "")
      }
    }
    tau_prior_temp <- paste(tau_list,
                            collapse = "\n",
                            sep = "")
    tau_prior <-      paste(tau_prior_temp, "\n",
                            "Tau[1,2] <- 0\n",
                            "Tau[2,1] <- 0\n",
                            "Tau[3,1] <- 0\n",
                            "Tau[1,3] <- 0\n",
                            "Tau[1,4] <- 0\n",
                            "Tau[4,1] <- 0\n",
                            "Tau[2,3] <- 0\n",
                            "Tau[3,2] <- 0\n",
                            "Tau[2,4] <- 0\n",
                            "Tau[4,2] <- 0\n",
                            "Tau[3,4] <- 0\n",
                            "Tau[4,3] <- 0\n",
                            collapse = "\n",
                            sep = "")
    # end: tau priors
    #-----


    #-----
    # start: beta priors
    if(!any(prior_names %in% "beta")){
      prior$beta <- 1
    }
    beta_prior <- prior$beta
    beta_list <- list()

    user_prior <- as.numeric(unlist(lapply(
      1:length(beta_prior),
      function(x) numextract(names(beta_prior)[x])
    ))) + 1

    if(any(user_prior > betas_loc)){
      stop("beta not found")
    }


    for(i in 1:betas_loc){

      if(any(names(beta_prior) ==  paste("beta_", i - 1,  sep = "") )){
        beta_list[[i]] <- paste("beta[",i,"]"," ~ ",
                                beta_prior[[paste("beta_", i - 1,  sep = "")]],
                                "\n", sep = "")
      } else {
        if(i == 1){
          beta_list[[i]] <- paste("beta[1] ~ ", "dnorm(", mean_y, ",0.001) \n", sep = "")
        } else{
          beta_list[[i]] <- paste("beta[", i, "] ~ ", "dnorm(", 0, ",0.001)\n", sep = "")
        }
      }
    }

    beta_prior <- paste(beta_list, collapse = "")
    # end: beta prior
    #-----


    #-----
    # start: eta priors
    if (!any(prior_names %in% "eta")) {
      prior$eta <- 1
    }

    eta_prior <- prior$eta
    eta_list <- list()

    user_prior <- as.numeric(unlist(lapply(
      1:length(eta_prior),
      function(x) numextract(names(eta_prior)[x])
    ))) + 1

    if(any(user_prior > etas_scl)){
      stop("eta not found")
    }
    for(i in 1:etas_scl){
      if(any(names(eta_prior) ==  paste("eta_", i -1,  sep = "") )){
        eta_list[[i]] <- paste("eta[",i,"]"," ~ ",
                               eta_prior[[paste("eta_", i - 1,  sep = "")]],
                               "\n", sep = "")
      } else {
        eta_list[[i]] <- paste("eta[", i, "] ~ ", "dnorm(", 0, ",0.001)\n", sep = "")
      }
    }
    eta_prior <- paste(eta_list, collapse = "")
    # end: eta prior
    #-----


    rho_test <- unlist(rho_test)
    cors <- NA

    if(rho_test == "muvar"){
      cors <- c(2,3,4,5)
      not_tested_cors <- which(1:6 %in% cors == FALSE)
    }
    if (rho_test == "all"){
      cors <- 1:6
      not_tested_cors <- 0
    }

    # start: build rho prior
    if(any(prior_names %in% "rho")){

      rho_priors <- prior$rho
    } else {
      rho_priors <- list(spike = 0.01, slab = 0.5)
    }


    if(mixture == "SSVS"){

      if(!all(c("spike", 'slab') %in% names(rho_priors))) {
        stop("spike and slab standard deviation must be specified")
      }

      if(k == 0 | is.null(k)){
        stop("k (number of mixture components must be specified)")
      }

      if (k == 2) {# start: two components (k = 2)

        rho_prior <- paste(
          "for(k in cors){\n",
          "k_rho[k] ~ dcat(pr_k[])\n",
          "z[k] ~ dnorm(0, rho_sd[k_rho[k]])\n",
          "rho[k] <- tanh(z[k])\n}\n",
          paste("rho_sd[1] <- pow(", rho_priors$spike , ", -2)\n"),
          paste("rho_sd[2] <- pow(", rho_priors$slab , ", -2)\n"),
          "pr_k[1] <- 1/2\n",
          "pr_k[2] <- 1/2\n"
        )

        if(any(not_tested_cors != 0)){
          rho_prior <- paste("# rho(s) not tested \n",
                             "for(k in not_tested_cors){\n",
                             "z[k] ~ dnorm(0, pow(",   rho_priors$slab, ", -2)) \n",
                             "rho[k] <- tanh(z[k])\n}\n",
                             "# rho test \n",
                             rho_prior, collapse = "\n\n")
        }


      }# end: k = 2

      # start: k = 3
      if(k == 3){

        rho_prior <-  paste(
          "for(k in cors){\n",
          "k_rho[k] ~ dcat(pr_k[])\n",
          "z[k] ~ dnorm(0,  rho_sd[k_rho[k]])T(T1[k_rho[k]], T2[k_rho[k]])\n",
          "rho[k] <- tanh(z[k]) \n}\n\n",
          paste("rho_sd[1] <- pow(", rho_priors$spike, ",-2)\n"),
          paste("rho_sd[2] <- pow(", rho_priors$slab, ",-2)\n"),
          paste("rho_sd[3] <- pow(", rho_priors$slab, ",-2)\n"),
          "pr_k[1] <- 1/3\n",
          "pr_k[2] <- 1/3\n",
          "pr_k[3] <- 1/3\n\n",
          "T1[1] <- -10\n",
          "T1[2] <- 0\n",
          "T1[3] <- -10\n",
          "T2[1] <- 10\n",
          "T2[2] <- 10\n",
          "T2[3] <- 0\n"
        )

        if(any(not_tested_cors != 0)){
          rho_prior <- paste("# rho(s) not tested \n",
                             "for(k in not_tested_cors){\n",
                             "z[k] ~ dnorm(0, pow(", rho_priors$slab, ", -2)) \n",
                             "rho[k] <- tanh(z[k])\n}\n",
                             "# rho test \n",
                             rho_prior, collapse = "\n\n")
        }
      } # end: k = 3
      # end: SSVS
    } else if(mixture == "KM"){# start: kuo & mallick (1998)

      rho_prior <- paste(
        "for(k in cors){\n",
        "k_rho[k] ~ dbern(0.5)\n",
        paste("z[k] ~ dnorm(0, pow(", rho_priors$slab, ", -2))\n"),
        "rho[k] <- tanh(z[k]) * k_rho[k] \n}\n"
      )
      if(any(not_tested_cors != 0)){
        rho_prior <- paste("# rho(s) not tested \n",
                           "for(k in not_tested_cors){\n",
                           "z[k] ~ dnorm(0, pow(",   rho_priors$slab, ", -2)) \n",
                           "rho[k] <- tanh(z[k])\n}\n",
                           "# rho test \n",
                           rho_prior, collapse = "\n\n")
      }
      # end: kuo and mallick (1998)
    }
    else{
      rho_prior <- paste(
        "# no cors tested\n",
        "for(k in 1:6){\n",
        paste("   z[k] ~ dnorm(0, pow(", rho_priors$slab, ", -2))\n"),
        "   rho[k] <- tanh(z[k])\n}\n", sep = ""
      )
    }



    mod <- paste("model{\n",
                 likelihood, "\n",
                 "# beta priors\n",
                 beta_prior, "\n",
                 "# eta priors\n",
                 eta_prior, "\n",
                 "# tau priors\n",
                 tau_prior, "\n\n",
                 "# rho priors\n",
                 rho_prior, "\n}\n")

    dat_list_scl <- lapply(2:etas_scl, function(x) fe_scl_mm[,x])
    dat_list_loc <- lapply(2:betas_loc, function(x) fe_loc_mm[,x])
    names(dat_list_scl) <- paste("X_scl_", 2:etas_scl, sep = "")
    names(dat_list_loc) <- paste("X_loc_", 2:betas_loc, sep = "")
    dat_list <- c(dat_list_loc, dat_list_scl)

    dat_list$y <- y
    dat_list$J <- length(unique(data[,split_re_loc[2]]))
    dat_list$N <- nrow(data)
    dat_list$ID <- data[,split_re_loc[2]]
    dat_list$cors <- cors

    if (length(not_tested_cors) > 0) {

      dat_list$not_tested_cors <- not_tested_cors

      }

    suppressWarnings(m_initialize <- jags.model(textConnection(mod),
                                               data = dat_list,
                                               n.chains = chains))

    update(m_initialize, n.iter = adapt)

    suppressWarnings(m_samps <- coda.samples(m_initialize,
                                              thin = thin,
                            variable.names = c("beta", "eta", "u",
                                               "rho", "Tau", "k_rho"),
                            n.iter = iter))
     ret <- list(melsm_fit = m_samps,
                 melsm_initialize = m_initialize,
                 jags_model = mod,
                 model_matrices = model_matrices,
                 dat = data,
                 cluster = split_re_loc[2],
                 iter = iter,
                 chains = chains,
                 thin = thin,
                 mixture = mixture,
                 rho_test  = rho_test,
                 k = k, call = match.call())
     }
  returned_object <- ret
  class(returned_object) <- "melsm"
  return(returned_object)
}

#' @title S3 \code{melsm} method
#' @inherit melsm.default
#'
#' @description S3 \code{melsm} method
#' @seealso \code{\link{melsm.default}}
#' @export
melsm <- function(fixed_location,
                  random_location,
                  fixed_scale,
                  random_scale,
                  prior = NULL,
                  mixture = NULL,
                  k = 2,
                  rho_test = "muvar",
                  adapt = 1000,
                  chains = 4,
                  iter = 5000,
                  thin = 1,
                  data, ...) {
  UseMethod("melsm")
}



#' Summary Method for \code{melsm} Objects
#'
#' @param object \code{melsm} object
#' @param cred credible interval
#' @param ... currently ignored
#'
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
#' summary(fit)
#'
#'}
summary.melsm <- function(object, cred = 0.95, ...){
  samps <- do.call(rbind.data.frame, object$melsm_fit)
  cat(paste("hypMuVar: Bayesian Hypothesis Testing",
            "\n          of Mean--Variance Relations\n"))
  cat("Model: MELSM\n")
  cat("Chains:", object$chains, "\n")
  if(isTRUE(object$bma)){
    cat("Samples:", object$bma_samps, "\n")
  } else {
  cat("Samples:",  (object$iter * object$chains) / object$thin , "\n")
  }
  cat("Mixture:", object$mixture, "\n")
  if(object$mixture == "none"){
    cat("Rho Test:", "none\n")
  } else {
    cat("Rho Test:", object$rho_test, "\n")
  }
  cat("Credible Interval:", cred, "\n")
  cat("Bayesian Model Averaged:", "(Top", object$top, "Models)", "\n")
  cat("----\n")
  cat("Call:\n")
  print(object$call)
  cat("----\n\n")
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


#' @title Print Method for \code{melsm} objects
#'
#' @param x \code{melsm} object
#' @param ... currently ignored
#' @export
print.melsm <- function(x, ...){
  summary(x)
  # object <- x
  # samps <- do.call(rbind.data.frame, object$melsm_fit)
  # cat(paste("hypMuVar: Bayesian Hypothesis Testing",
  #           "\n          of Mean--Variance Relations\n"))
  # cat("Model: MELSM\n")
  # cat("Chains:", object$chains, "\n")
  # cat("Samples:",  (object$iter * object$chains) / object$thin , "\n")
  # cat("Mixture:", object$mixture, "\n")
  # if(object$mixture == "none"){
  #   cat("Rho Test:", "none\n")
  # } else {
  #   cat("Rho Test:", object$rho_test, "\n")
  # }
  # cat("----\n")
}


