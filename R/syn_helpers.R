
#'
#'
#'

syn_part_1 = function (fit_notw, scale_factor, shift_factor) {
  log_lik_notw <- log_lik(fit_notw)
  
  res_notw <- list(log_lik = log_lik_notw)
  
  loglike_notw <- log_like_hamming(res = res_notw, thresh = 1.00,
                                   risks_method = "linear",
                                   scale_factor = scale_factor,
                                   shift_factor = shift_factor,
                                   LipschitzPlot = FALSE, weightsPlot = FALSE)
  
  return(list(loglike_notw$Lipschitz_by_record, loglike_notw$weights))
}

#'
#'
#'

syn_part_2 = function (fit_LW, scale_factor, shift_factor, weights_LW, c) {
  log_lik_LW <- log_lik(fit_LW)
  
  res_LW <- list(log_lik = log_lik_LW)
  
  loglike_LW <- log_like_hamming(res = res_LW, thresh = 1.00,
                                 risks_method = "linear",
                                 scale_factor = scale_factor,
                                 shift_factor = shift_factor,
                                 LipschitzPlot = FALSE, weightsPlot = FALSE)
  
  L_LW_all        <- loglike_LW$log_ratio_theta
  L_LW            <- loglike_LW$L
  
  weights_LW_final <- weights_LW * L_LW / L_LW_all * c
  weights_LW_final[weights_LW_final > 1] <- 1
  weights_LW_final[weights_LW_final < 0] <- 0
  
  return(list(loglike_LW$Lipschitz_by_record, weights_LW_final))
}


#'
#'
#'

syn_part_3 = function (fit_LW_final, scale_factor, shift_factor) {
  log_lik_LW_final <- log_lik(fit_LW_final)
  
  res_LW_final <- list(log_lik = log_lik_LW_final)
  
  loglike_LW_final <- log_like_hamming(res = res_LW_final, thresh = 1.00,
                                       risks_method = "linear",
                                       scale_factor = scale_factor,
                                       shift_factor = shift_factor,
                                       LipschitzPlot = FALSE, weightsPlot = FALSE)

  Lipschitz_LW_final    <- loglike_LW_final$Lipschitz_by_record
  
  return(list(loglike_LW_final$Lipschitz_by_record))
}

#'
#'
#'

syn_plots = function (Lipschitz_notw, Lipschitz_LW, Lipschitz_LW_final, weights_LW, weights_LW_final) {
  defaultW <- getOption("warn")
  options(warn = -1)
  
  cbPalette_withoutdata <- c("#999999",  "#E69F00", "#56B4E9",  "#009E73", "#CC79A7", "#D55E00",  "#F0E442")
  
  require(tidyverse)
  Lbounds = cbind(Lipschitz_notw,
                  Lipschitz_LW,
                  Lipschitz_LW_final) %>% as_tibble()
  names(Lbounds) = c("Unweighted", "LW", "LW_final")
  Lbounds_long = reshape2::melt(Lbounds)
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

