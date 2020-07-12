
#'
#'
#' @import Rfast
#' @import superheat
#' @import tidyverse
#' @export

log_like_hamming = function(res, thresh = 1.00, risks_method = "linear",
                            LipschitzPlot = FALSE, scale_factor = 1.0,
                            shift_factor = 0.0, weightsPlot = FALSE)
{
  require(Rfast) ## needed for rowMaxs() and colMaxs()
  require(superheat) ## needed for heatmap
  log_lik           <- res$log_lik ## S x N matrix, including weighting
  N                 <- ncol(log_lik)
  S                 <- nrow(log_lik)

  log_ratio         <- matrix(0,S,N)
  log_ratio_theta   <- matrix(0,S,1)
  pos               <- rep(TRUE,N)

  for( s in 1:S )
  {
    log_like_xs  <- sum(log_lik[s,]) ## full data
    for(i in 1:N)
    {
      pos_i               <- pos
      pos_i[i]            <- FALSE
      log_like_xsi        <- sum(log_lik[s,pos_i])
      log_ratio[s,i]      <- abs(log_like_xs - log_like_xsi) ## L in 3.1(b)
    } ## leave-one-out neighborhood
  } ## end loop over S MCMC iterations
  log_ratio_theta         <- rowMaxs(log_ratio, value = TRUE) ## S x 1 (by row) maxima of log-ratio over records for *each* MCMC draw
  L                       <- quantile(log_ratio_theta,thresh)
  S_excl                  <- which(log_ratio_theta > L) ## Excluded MCMC iterations based on thres log_ratio_theta

  ## remove rows from log_ratio for excluded MCMC observations so that weights tie to L
  if( !is.null(S_excl) && length(S_excl) > 0 )
  {
    logthresh_ratio         <- log_ratio[-S_excl,]
  }else{
    logthresh_ratio         <- log_ratio
  }

  ## Compute by record risks \in (0,1) - only include MCMC draws used to compose L
  Lipschitz_by_record       <- colMaxs(logthresh_ratio,value=TRUE) ## 1 x N (by col) maxima of log-ratio over records for *each* record
  if( risks_method == "linear")
  {
    f_linres  <- function(x){(x-min(x))/(max(x)-min(x))}
    risks     <- f_linres( Lipschitz_by_record ) ## \in (0,1)
  }else{
    if( risks_method == "softmax")
    {
      f_softmax <- function(x){(1.0/(1.0+exp(-1.69897*(x-mean(x))/sd(x))))}
      risks     <- f_softmax( Lipschitz_by_record ) ## \in (0,1)
    }else{ ## weight_method == "logistic"
      f_logistic <- function(x,factor=0.5){f_x = (1/(1+exp(-factor*x)))*2-1}
      risks      <- f_logistic( Lipschitz_by_record ) ## \in (0,1)
    }
  } ## end condition on which weight method to apply to convert L_i from R^+ to (0,1)

  if (LipschitzPlot) {
    require(tidyverse)
    dat_lr  <-  Lipschitz_by_record %>% as_tibble()
    names(dat_lr) <- c("Lipschitz")
    ## violin plot
    print(
      ggplot(dat_lr, aes(y=Lipschitz, x = "")) +
        geom_violin(trim=TRUE, alpha = 0.3) +
        theme_bw(base_size = 15) +
        ylab("Record-level Lipschitz bounds") + xlab("")
    )
  }

  weights = scale_factor * (1 - risks) + shift_factor
  weights[weights <= 0] = 0
  weights[weights >= 1] = 1

  if (weightsPlot) {
    require(tidyverse)
    dat_lr  <-  weights %>% as_tibble()
    names(dat_lr) <- c("Weights")
    ## violin plot
    print(
      ggplot(dat_lr, aes(y=Weights, x= "")) +
        geom_violin(trim=TRUE, alpha = 0.3) +
        theme_bw(base_size = 15) +
        ylab("Record-level Weights") + xlab("")
    )
  }

  return(list(log_ratio = log_ratio, logthresh_ratio = logthresh_ratio, log_ratio_theta = log_ratio_theta,
              Lipschitz_by_record = Lipschitz_by_record, S_excl = S_excl, L = L, risks = risks, weights = weights))

}


#'
#' @import rstanarm
#' @import gridExtra
#' @import tidyverse
#' @import reshape2
#' @export

syn = function(origonal_data, model_string = "outcome ~ 1", chains = 1, iterations = 1000,
               scale_factor = 1.0, shift_factor = 0.0, c = 0.95, thresh = 1.00,
               risks_method = "linear", m = 20, thin = 5, plots = FALSE) {

  require(rstanarm)
  ###########################################################################
  ######################## unweighted, fit, and DP risks ####################
  ###########################################################################
  fit_notw <- stan_glm(
    model_string,
    data = origonal_data,
    family = poisson(link="log"),
    prior = normal(0, 2, autoscale = FALSE),
    refresh = 0,
    chains = chains, iter = iterations
  )

  log_lik_notw <- log_lik(fit_notw)

  res_notw <- list(log_lik = log_lik_notw)

  loglike_notw <- log_like_hamming(res = res_notw, thresh = 1.00,
                                         risks_method = "linear",
                                         scale_factor = scale_factor,
                                         shift_factor = shift_factor,
                                         LipschitzPlot = FALSE, weightsPlot = FALSE)

  risks_notw        <- loglike_notw$risks
  Lipschitz_notw    <- loglike_notw$Lipschitz_by_record

  ###########################################################################
  ################## LW weights, fit, and DP risks ##########################
  ###########################################################################

  weights_LW <- loglike_notw$weights

  assign("weights_LW", weights_LW, sys.nframe())

  fit_LW <- stan_glm(
    model_string,
    data = origonal_data,
    family = poisson(link="log"),
    weights = weights_LW,
    prior = normal(0, 2, autoscale = FALSE),
    refresh = 0,
    chains = chains, iter = iterations
  )

  rm(list = c("weights_LW"), envir = globalenv())

  log_lik_LW <- log_lik(fit_LW)

  res_LW <- list(log_lik = log_lik_LW)

  loglike_LW <- log_like_hamming(res = res_LW, thresh = 1.00,
                                       risks_method = "linear",
                                       scale_factor = scale_factor,
                                       shift_factor = shift_factor,
                                       LipschitzPlot = FALSE, weightsPlot = FALSE)
  risks_LW       <- loglike_LW$risks
  Lipschitz_LW    <- loglike_LW$Lipschitz_by_record

  L_LW_all        <- loglike_LW$log_ratio_theta
  L_LW            <- loglike_LW$L
  ###########################################################################
  ################## LW weights reweight, fit, DP risks, and synthesize #####
  ###########################################################################
  weights_LW_final <- weights_LW * L_LW / L_LW_all * c
  weights_LW_final[weights_LW_final > 1] <- 1
  weights_LW_final[weights_LW_final < 0] <- 0

  assign("weights_LW_final", weights_LW_final, globalenv())

  fit_LW_final <- stan_glm(
    model_string,
    data = origonal_data,
    family = poisson(link="log"),
    weights = weights_LW_final,
    prior = normal(0, 2, autoscale = FALSE),
    refresh = 0,
    chains = chains, iter = iterations
  )

  rm(list = c("weights_LW_final"), envir = globalenv())

  log_lik_LW_final <- log_lik(fit_LW_final)

  res_LW_final <- list(log_lik = log_lik_LW_final)

  loglike_LW_final <- log_like_hamming(res = res_LW_final, thresh = 1.00,
                                             risks_method = "linear",
                                             scale_factor = scale_factor,
                                             shift_factor = shift_factor,
                                             LipschitzPlot = FALSE, weightsPlot = FALSE)
  risks_LW_final        <- loglike_LW_final$risks
  Lipschitz_LW_final    <- loglike_LW_final$Lipschitz_by_record
  #### synthesis ####
  N = length(origonal_data[,1])
  draws <- as.data.frame(fit_LW_final)
  start <- 500 / 2
  syndata <- vector("list", m)
  for (i in 1:m){
    draws_exp <- exp(draws[start + thin * (i - 1), ])
    syndata[[i]] <- rpois(N, draws_exp)
  }

  ###########################################################################
  ################## Plots ##################################################
  ###########################################################################

  cbPalette_withoutdata <- c("#999999",  "#E69F00", "#56B4E9",  "#009E73", "#CC79A7", "#D55E00",  "#F0E442")

  if (plots) {
    defaultW <- getOption("warn")
    options(warn = -1)

    require(tidyverse)
    Lbounds = cbind(Lipschitz_notw,
                    Lipschitz_LW,
                    Lipschitz_LW_final) %>% as_tibble()
    names(Lbounds) = c("Unweighted", "LW", "LW_final")
    Lbounds_long = reshape2::melt(Lbounds, id.vars = NULL)
    plot_L = ggplot(Lbounds_long, aes(x = variable, y = value, fill = variable, color = variable)) +
      geom_violin(trim=TRUE, alpha = 0.3) +
      scale_colour_manual(values = cbPalette_withoutdata) + scale_fill_manual(values = cbPalette_withoutdata) +
      theme_bw(base_size = 15)   +
      theme(legend.position = "none")  +
      ylab("Lipschitz Bounds") + xlab("")
    print(plot_L)

    df1 <- data.frame(weights_LW, Lipschitz_LW)
    names(df1) <- c("weights", "Lipchitz")
    plot_weightsL_1 <- ggplot(df1, aes(x = weights, y = Lipchitz)) +
      geom_point() +
      theme_bw(base_size = 15) + ggtitle("LW") +
      geom_hline(yintercept=max(Lipschitz_LW),
                 linetype="dashed", color = "red", size=2) +
      xlim(0, 1)

    df2 <- data.frame(weights_LW_final, Lipschitz_LW_final)
    names(df2) <- c("weights", "Lipchitz")
    plot_weightsL_2 <- ggplot(df2, aes(x = weights, y = Lipchitz)) +
      geom_point() +
      theme_bw(base_size = 15) + ggtitle("LW_final") +
      geom_hline(yintercept=max(Lipschitz_LW_final),
                 linetype="dashed", color = "red", size=2) +
      xlim(0, 1)

    require(gridExtra)
    grid.arrange(plot_weightsL_1, plot_weightsL_2, ncol = 2)
    options(warn = defaultW)
  }

  return(syndata)
}
