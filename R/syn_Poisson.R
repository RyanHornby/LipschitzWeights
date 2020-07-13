#'
#' @import rstanarm
#' @import gridExtra
#' @import tidyverse
#' @import reshape2
#' @export

syn_Poisson = function(origonal_data, model_string = "outcome ~ 1", chains = 1, iterations = 1000,
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
  
  syn1 = syn_part_1(fit_notw, scale_factor, shift_factor)
  
  # log_lik_notw <- log_lik(fit_notw)
  # 
  # res_notw <- list(log_lik = log_lik_notw)
  # 
  # loglike_notw <- log_like_hamming(res = res_notw, thresh = 1.00,
  #                                  risks_method = "linear",
  #                                  scale_factor = scale_factor,
  #                                  shift_factor = shift_factor,
  #                                  LipschitzPlot = FALSE, weightsPlot = FALSE)
  # 
  # risks_notw        <- loglike_notw$risks
  # Lipschitz_notw    <- loglike_notw$Lipschitz_by_record
  
  ###########################################################################
  ################## LW weights, fit, and DP risks ##########################
  ###########################################################################
  
  weights_LW <- syn1[[2]]
  
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
  
  syn2 = syn_part_2(fit_LW, scale_factor, shift_factor, weights_LW, c)
  
  # log_lik_LW <- log_lik(fit_LW)
  # 
  # res_LW <- list(log_lik = log_lik_LW)
  # 
  # loglike_LW <- log_like_hamming(res = res_LW, thresh = 1.00,
  #                                risks_method = "linear",
  #                                scale_factor = scale_factor,
  #                                shift_factor = shift_factor,
  #                                LipschitzPlot = FALSE, weightsPlot = FALSE)
  # risks_LW       <- loglike_LW$risks
  # Lipschitz_LW    <- loglike_LW$Lipschitz_by_record
  # 
  # L_LW_all        <- loglike_LW$log_ratio_theta
  # L_LW            <- loglike_LW$L
  ###########################################################################
  ################## LW weights reweight, fit, DP risks, and synthesize #####
  ###########################################################################
  # weights_LW_final <- weights_LW * L_LW / L_LW_all * c
  # weights_LW_final[weights_LW_final > 1] <- 1
  # weights_LW_final[weights_LW_final < 0] <- 0
  
  weights_LW_final = syn2[[2]]
  
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
  
  syn3 = syn_part_3(fit_LW_final, scale_factor, shift_factor)
  
  # log_lik_LW_final <- log_lik(fit_LW_final)
  # 
  # res_LW_final <- list(log_lik = log_lik_LW_final)
  # 
  # loglike_LW_final <- log_like_hamming(res = res_LW_final, thresh = 1.00,
  #                                      risks_method = "linear",
  #                                      scale_factor = scale_factor,
  #                                      shift_factor = shift_factor,
  #                                      LipschitzPlot = FALSE, weightsPlot = FALSE)
  # risks_LW_final        <- loglike_LW_final$risks
  # Lipschitz_LW_final    <- loglike_LW_final$Lipschitz_by_record
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
  
  
  
  if (plots) {
    syn_plots(syn1[[1]], syn2[[1]], syn3[[1]], syn1[[2]], syn2[[2]])
    # defaultW <- getOption("warn")
    # options(warn = -1)
    # 
    # require(tidyverse)
    # Lbounds = cbind(syn1[[1]],
    #                 syn2[[1]],
    #                 syn3[[1]]) %>% as_tibble()
    # names(Lbounds) = c("Unweighted", "LW", "LW_final")
    # Lbounds_long = reshape2::melt(Lbounds, id.vars = NULL)
    # plot_L = ggplot(Lbounds_long, aes(x = variable, y = value, fill = variable, color = variable)) +
    #   geom_violin(trim=TRUE, alpha = 0.3) +
    #   scale_colour_manual(values = cbPalette_withoutdata) + scale_fill_manual(values = cbPalette_withoutdata) +
    #   theme_bw(base_size = 15)   +
    #   theme(legend.position = "none")  +
    #   ylab("Lipschitz Bounds") + xlab("")
    # print(plot_L)
    # 
    # df1 <- data.frame(weights_LW, syn2[[1]])
    # names(df1) <- c("weights", "Lipchitz")
    # plot_weightsL_1 <- ggplot(df1, aes(x = weights, y = Lipchitz)) +
    #   geom_point() +
    #   theme_bw(base_size = 15) + ggtitle("LW") +
    #   geom_hline(yintercept=max(syn2[[1]]),
    #              linetype="dashed", color = "red", size=2) +
    #   xlim(0, 1)
    # 
    # df2 <- data.frame(weights_LW_final, syn3[[1]])
    # names(df2) <- c("weights", "Lipchitz")
    # plot_weightsL_2 <- ggplot(df2, aes(x = weights, y = Lipchitz)) +
    #   geom_point() +
    #   theme_bw(base_size = 15) + ggtitle("LW_final") +
    #   geom_hline(yintercept=max(syn3[[1]]),
    #              linetype="dashed", color = "red", size=2) +
    #   xlim(0, 1)
    # 
    # require(gridExtra)
    # grid.arrange(plot_weightsL_1, plot_weightsL_2, ncol = 2)
    # options(warn = defaultW)
  }
  
  return(syndata)
}
