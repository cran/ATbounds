## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ATbounds)

## -----------------------------------------------------------------------------
  nsw_treated <- read.table("http://users.nber.org/~rdehejia/data/nsw_treated.txt")
  colnames(nsw_treated) <- c("treat","age","edu","black","hispanic",
                           "married","nodegree","RE75","RE78")

  nsw_control <- read.table("http://users.nber.org/~rdehejia/data/nsw_control.txt")
  colnames(nsw_control) <- c("treat","age","edu","black","hispanic",
                           "married","nodegree","RE75","RE78")


## -----------------------------------------------------------------------------
  nsw <- rbind(nsw_treated,nsw_control)
  attach(nsw)
  D <- treat  
  Y <- (RE78 > 0) 

## -----------------------------------------------------------------------------
  rps <- rep(mean(D),length(D))  

## -----------------------------------------------------------------------------
  ate_nsw <- mean(D*Y)/mean(D)-mean((1-D)*Y)/mean(1-D)
  print(ate_nsw)

## -----------------------------------------------------------------------------
  model <- lm(Y ~ D)
  summary(model)
  confint(model)

## -----------------------------------------------------------------------------
  detach(nsw)
  nswre_treated <- read.table("http://users.nber.org/~rdehejia/data/nswre74_treated.txt")
  colnames(nswre_treated) <- c("treat","age","edu","black","hispanic",
                           "married","nodegree","RE74","RE75","RE78")

  nswre_control <- read.table("http://users.nber.org/~rdehejia/data/nswre74_control.txt")
  colnames(nswre_control) <- c("treat","age","edu","black","hispanic",
                           "married","nodegree","RE74","RE75","RE78")
  nswre <- rbind(nswre_treated,nswre_control)
  attach(nswre)
  D <- treat  
  Y <- (RE78 > 0) 
  X <- cbind(age,edu,black,hispanic,married,nodegree,RE74/1000,RE75/1000)

## -----------------------------------------------------------------------------
  rps <- rep(mean(D),length(D))  

## -----------------------------------------------------------------------------
  ate_nswre <- mean(D*Y)/mean(D)-mean((1-D)*Y)/mean(1-D)
  print(ate_nswre)

## -----------------------------------------------------------------------------
  model <- lm(Y ~ D)
  summary(model)
  confint(model)

## -----------------------------------------------------------------------------
  bns_nsw <- atebounds(Y, D, X, rps)

## -----------------------------------------------------------------------------
  summary(bns_nsw)

## -----------------------------------------------------------------------------
summary(atebounds(Y, D, X, rps, Q = 2))

## -----------------------------------------------------------------------------
summary(atebounds(Y, D, X, rps, Q = 4))

## -----------------------------------------------------------------------------
print(ate_nswre)

## -----------------------------------------------------------------------------
summary(atebounds(Y, D, X, rps))
summary(atebounds(Y, D, X, rps, n_hc = ceiling(length(Y)/5)))
summary(atebounds(Y, D, X, rps, n_hc = ceiling(length(Y)/20)))

## -----------------------------------------------------------------------------
  print(bns_nsw)

## -----------------------------------------------------------------------------
  bns_nsw_att <- attbounds(Y, D, X, rps)
  summary(bns_nsw_att)

## -----------------------------------------------------------------------------
  summary(attbounds(Y, D, X, rps, Q = 2))

## -----------------------------------------------------------------------------
  summary(attbounds(Y, D, X, rps, Q = 4))

## -----------------------------------------------------------------------------
summary(attbounds(Y, D, X, rps))
summary(attbounds(Y, D, X, rps, n_hc = ceiling(length(Y)/5)))
summary(attbounds(Y, D, X, rps, n_hc = ceiling(length(Y)/20)))

## -----------------------------------------------------------------------------
  psid2_control <- read.table("http://users.nber.org/~rdehejia/data/psid2_controls.txt")
  colnames(psid2_control) <- c("treat","age","edu","black","hispanic",
                           "married","nodegree","RE74","RE75","RE78")
  psid <- rbind(nswre_treated,psid2_control)
  detach(nswre)
  attach(psid)
  D <- treat  
  Y <- (RE78 > 0) 
  X <- cbind(age,edu,black,hispanic,married,nodegree,RE74/1000,RE75/1000)

## -----------------------------------------------------------------------------
  rps_sp <- rep(mean(D),length(D))  
  bns_psid <- atebounds(Y, D, X, rps_sp)
  summary(bns_psid)

## -----------------------------------------------------------------------------
  summary(atebounds(Y, D, X, rps_sp, Q=1))

## -----------------------------------------------------------------------------
  summary(attbounds(Y, D, X, rps_sp))
  detach(psid)

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

## ---- eval=FALSE--------------------------------------------------------------
#    # Bounding  ATT: sensitivity analysis
#    # not run to save time
#    nhc_set <- c(5, 10, 20)
#    results_att <- {}
#  
#    for (hc in nhc_set){
#      nhc <- ceiling(length(Y)/hc)
#  
#      for (q in c(1,2,3,4)){
#        res <- attbounds(Y, D, X, rps, Q = q, n_hc = nhc)
#        results_att <- rbind(results_att,c(hc,q,res$est_lb,res$est_ub,res$ci_lb,res$ci_ub))
#      }
#    }
#    colnames(results_att) = c("L","Q","LB","UB","CI-LB","CI-UB")
#    print(results_att, digits = 3)

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

