
#'
#'
#' @import Rfast
#' @import superheat
#' @export

log_like_hamming = function(res, thresh = 0.975, risks_method = "linear",
                            riskPlot = FALSE, scale_factor = 1.0,
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
  log_ratio_data       <- colMaxs(logthresh_ratio,value=TRUE) ## 1 x N (by col) maxima of log-ratio over records for *each* record
  if( risks_method == "linear")
  {
    f_linres  <- function(x){(x-min(x))/(max(x)-min(x))}
    risks     <- f_linres( log_ratio_data ) ## \in (0,1)
  }else{
    if( risks_method == "softmax")
    {
      f_softmax <- function(x){(1.0/(1.0+exp(-1.69897*(x-mean(x))/sd(x))))}
      risks     <- f_softmax( log_ratio_data ) ## \in (0,1)
    }else{ ## weight_method == "logistic"
      f_logistic <- function(x,factor=0.5){f_x = (1/(1+exp(-factor*x)))*2-1}
      risks      <- f_logistic( log_ratio_data ) ## \in (0,1)
    }
  } ## end condition on which weight method to apply to convert L_i from R^+ to (0,1)

  if (riskPlot) {
    dat_lr  <-  risks %>% as_tibble()
    names(dat_lr) <- c("Risks")
    ## density plot
    print(
      ggplot(dat_lr, aes(Risks)) +
        geom_density(alpha=0.1, adjust = 1.5) + theme_bw() +
        xlab("Record level Risks in (0,1)") + ylab("Density")
    )
  }

  weights = scale_factor * (1 - risks) + shift_factor
  weights[weights <= 0] = 0
  weights[weights >= 1] = 1

  if (weightsPlot) {
    dat_lr  <-  weights %>% as_tibble()
    names(dat_lr) <- c("Weights")
    ## density plot
    print(
      ggplot(dat_lr, aes(Weights)) +
        geom_density(alpha=0.1, adjust = 1.5) + theme_bw() +
        xlab("Record level Weights") + ylab("Density")
    )
  }

  return(list(log_ratio = log_ratio, logthresh_ratio = logthresh_ratio, log_ratio_theta = log_ratio_theta,
              log_ratio_data = log_ratio_data, S_excl = S_excl, L = L, risks = risks, weights = weights))

} ## end function log_like_hamming() to compute Lipschitz bound under Hamming-1 distance
