## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ATbounds)

## ----nsw-eval-options, include=FALSE------------------------------------------
knitr::opts_chunk$set(eval = FALSE, purl = FALSE)

## ----packaged-data-eval-options, include=FALSE, eval=TRUE, purl=TRUE----------
knitr::opts_chunk$set(eval = TRUE, purl = TRUE)

## -----------------------------------------------------------------------------
  Y <- RHC[,"survival"]
  D <- RHC[,"RHC"]
  X <- as.matrix(RHC[,-c(1,2)])

## -----------------------------------------------------------------------------
  # Logit estimation of propensity score
  glm_ps <- stats::glm(D~X,family=binomial("logit"))
  ps <- glm_ps$fitted.values
  ps_treated <- ps[D==1]
  ps_control <- ps[D==0]
  # Plotting histograms of propensity scores
  df <- data.frame(cbind(D,ps))
  colnames(df)<-c("RHC","PS")
  df$RHC <- as.factor(df$RHC)
  levels(df$RHC) <- c("No RHC (Control)", "RHC (Treated)")
  
  ggplot2::ggplot(df, ggplot2::aes(x=PS, color=RHC, fill=RHC)) +
    ggplot2::geom_histogram(breaks=seq(0,1,0.1),alpha=0.5,position="identity")

## -----------------------------------------------------------------------------
  # ATT normalized estimation
  y1_att <- mean(D*Y)/mean(D)
  att_wgt <- ps/(1-ps)
  y0_att_num <- mean((1-D)*att_wgt*Y)
  y0_att_den <- mean((1-D)*att_wgt)
  y0_att <- y0_att_num/y0_att_den
  att_ps <- y1_att - y0_att
  print(att_ps)

## -----------------------------------------------------------------------------
  rps <- rep(mean(D),length(D))  

## -----------------------------------------------------------------------------
  att_rps <- mean(D*Y)/mean(D) - mean((1-D)*Y)/mean(1-D)
  print(att_rps)

## -----------------------------------------------------------------------------
   Xunique <- mgcv::uniquecombs(X) # A matrix of unique rows from X
   print(c("no. of unique rows:", nrow(Xunique))) 
   print(c("sample size       :", nrow(X)))   

## -----------------------------------------------------------------------------
  summary(attbounds(Y, D, X, rps))

## ----eval=FALSE---------------------------------------------------------------
#   # Bounding ATT: sensitivity analysis
#   # This chunk is not evaluated by default because the analysis is time-intensive.
#   nhc_set <- c(5, 10, 20)
#   results_att <- {}
# 
#   for (hc in nhc_set){
#     nhc <- ceiling(length(Y)/hc)
# 
#     for (q in c(1,2,3,4)){
#       res <- attbounds(Y, D, X, rps, Q = q, n_hc = nhc)
#       results_att <- rbind(results_att,c(hc,q,res$est_lb,res$est_ub,res$ci_lb,res$ci_ub))
#     }
#   }
#   colnames(results_att) = c("L","Q","LB","UB","CI-LB","CI-UB")
#   print(results_att, digits = 3)

## -----------------------------------------------------------------------------
  summary(atebounds(Y, D, X, rps, Q = 1))

## -----------------------------------------------------------------------------
  summary(atebounds(Y, D, X, rps, Q = 2))

## -----------------------------------------------------------------------------
  summary(atebounds(Y, D, X, rps, Q = 3))

## -----------------------------------------------------------------------------
  summary(atebounds(Y, D, X, rps, Q = 4))

## -----------------------------------------------------------------------------
  Y <- EFM[,"cesarean"]
  D <- EFM[,"monitor"]
  X <- as.matrix(EFM[,c("arrest", "breech", "nullipar", "year")])
  year <- EFM[,"year"] 

## -----------------------------------------------------------------------------
  ate_rps <- mean(D*Y)/mean(D) - mean((1-D)*Y)/mean(1-D)
  print(ate_rps)

## -----------------------------------------------------------------------------
  rps <- rep(mean(D),length(D))  
  print(rps[1])

## -----------------------------------------------------------------------------
  summary(atebounds(Y, D, X, rps, Q = 1, x_discrete = TRUE))

## -----------------------------------------------------------------------------
  summary(atebounds(Y, D, X, rps, Q = 2, x_discrete = TRUE))

## -----------------------------------------------------------------------------
  summary(atebounds(Y, D, X, rps, Q = 3, x_discrete = TRUE))

## -----------------------------------------------------------------------------
  summary(atebounds(Y, D, X, rps, Q = 5, x_discrete = TRUE))
  summary(atebounds(Y, D, X, rps, Q = 10, x_discrete = TRUE))
  summary(atebounds(Y, D, X, rps, Q = 20, x_discrete = TRUE))
  summary(atebounds(Y, D, X, rps, Q = 50, x_discrete = TRUE))
  summary(atebounds(Y, D, X, rps, Q = 100, x_discrete = TRUE))

