#'
#' @import rstanarm
#' @import gridExtra
#' @import tidyverse
#' @import reshape2
#' @export

syn_NegativeBinomial = function(origonal_data, model_string = "outcome ~ 1", chains = 1, iterations = 1000,
                      scale_factor = 1.0, shift_factor = 0.0, c = 0.95, thresh = 1.00,
                      risks_method = "linear", m = 20, thin = 5, plots = FALSE) {

  require(rstanarm)

  fit_notw <- stan_glm.nb(
    model_string,
    data = origonal_data,
    link = "log",
    prior_aux = exponential(1.5),
    refresh = 0,
    chains = chains, iter = iterations
  )

  syn1 = syn_part_1(fit_notw, scale_factor, shift_factor)

  weights_LW <- syn1[[2]]

  assign("weights_LW", weights_LW, sys.nframe())

  fit_LW <- stan_glm.nb(
    model_string,
    data = origonal_data,
    link = "log",
    prior_aux = exponential(1.5),
    weights = weights_LW,
    refresh = 0,
    chains = chains, iter = iterations
  )

  rm(list = c("weights_LW"), envir = globalenv())

  syn2 = syn_part_2(fit_LW, scale_factor, shift_factor, weights_LW, c)

  weights_LW_final = syn2[[2]]

  assign("weights_LW_final", weights_LW_final, globalenv())

  fit_LW_final <- stan_glm.nb(
    model_string,
    data = origonal_data,
    link = "log",
    prior_aux = exponential(1.5),
    weights = weights_LW_final,
    refresh = 0,
    chains = chains, iter = iterations
  )

  rm(list = c("weights_LW_final"), envir = globalenv())

  syn3 = syn_part_3(fit_LW_final, scale_factor, shift_factor)

  #### synthesis ####
  N = length(origonal_data[,1])
  draws <- as.data.frame(fit_LW_final)
  #start <- 500 / 2
  start = length(draws[,1]) - thin  * (m - 1)
  syndata <- vector("list", m)
  for (i in 1:m){
    draws_exp <- exp(draws[start + thin * (i - 1), 1])
    trials = draws[start + thin * (i - 1), 2]
    prob = trails/draws_exp
    syndata[[i]] <- rnbinom(N, trials, prob)
  }


  if (plots) {
    syn_plots(syn1[[1]], syn2[[1]], syn3[[1]], syn1[[2]], syn2[[2]])
  }

  return(syndata)
}
